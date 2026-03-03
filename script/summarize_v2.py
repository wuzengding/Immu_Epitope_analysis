import os
import shutil
import argparse
import pandas as pd
import numpy as np
import xlsxwriter
from Bio import SeqIO
from difflib import SequenceMatcher

# ==========================================
# 辅助函数：提取信息与富集区滑窗聚类
# ==========================================

def extract_mut_info(identity):
    parts = identity.split('|')
    gene = parts[2] if len(parts) > 2 else "-"
    hgvs = parts[4] if len(parts) > 4 and parts[4] != "-" else (parts[3] if len(parts) > 3 else "-")
    uniprot = parts[5] if len(parts) > 5 else "-"
    return gene, hgvs, uniprot

def find_mismatch(s1, s2):
    for i, (a, b) in enumerate(zip(s1, s2)):
        if a != b: return i
    return 0

def get_optimized_candidates(full_seq, epitopes, max_len=33, min_gap=10):
    """加权富集区滑窗聚类算法"""
    if not epitopes: return[]
    seq_len = len(full_seq)
    weights = np.zeros(seq_len)
    for ep in epitopes:
        w = 1.0 / (float(ep['rank']) + 0.01)
        weights[ep['start']:ep['end']] += w
    
    active = np.where(weights > 0)[0]
    if len(active) == 0: return []

    clusters =[]
    curr_c =[active[0]]
    for i in range(1, len(active)):
        if active[i] - active[i-1] <= min_gap:
            curr_c.append(active[i])
        else:
            clusters.append((min(curr_c), max(curr_c)))
            curr_c = [active[i]]
    clusters.append((min(curr_c), max(curr_c)))

    results =[]
    for c_start, c_end in clusters:
        c_len = c_end - c_start + 1
        if c_len <= max_len:
            pad = (max_len - c_len) // 2
            s = max(0, c_start - pad)
            e = min(seq_len, s + max_len)
            if e == seq_len: s = max(0, e - max_len)
            results.append(full_seq[s:e])
        else:
            best_s, max_w = c_start, 0
            for j in range(c_start, c_end - max_len + 2):
                win_w = np.sum(weights[j:j+max_len])
                if win_w > max_w:
                    max_w, best_s = win_w, j
            results.append(full_seq[best_s : best_s + max_len])
    return list(dict.fromkeys(results))

def filter_data(df, aff_col, rank_col, aff_limit, rank_limit, bind_levels, tm_drop_identities, whitelist=None):
    df = df.copy()
    bind_list =[x.strip() for x in bind_levels.split(',')]
    if whitelist is None: whitelist = set()
    
    def check(row):
        try:
            ident = str(row['Identity']).strip()
            if ident in tm_drop_identities:
                return "FAIL"
            if str(row.get('TransMemb', '')).strip() == 'TMhelix': 
                return "FAIL"
            if ident in whitelist: 
                return "PASS"
            if str(row['BindLevel']) not in bind_list: return "FAIL"
            if float(row[aff_col]) > aff_limit: return "FAIL"
            if float(row[rank_col]) > rank_limit: return "FAIL"
            return "PASS"
        except: return "FAIL"
    return df[df.apply(check, axis=1) == "PASS"]

def prepare_data(df_pass, var_fasta, ref_fasta, other_df_pass, is_mhci, gene_order=None):
    rows =[]
    df = df_pass.copy()
    df['gene_name'], df['hgvs_p'], df['uniprot_id'] = zip(*df['Identity'].apply(extract_mut_info))
    
    if gene_order:
        gene_dict = {g: i for i, g in enumerate(gene_order)}
        df['sort_idx'] = df['gene_name'].apply(lambda x: gene_dict.get(x, 999))
        df = df.sort_values(['sort_idx', 'Identity'])
    else:
        def get_custom_sort_key(identity):
            ident_upper = identity.upper()
            if 'LSTV' in ident_upper: type_priority = 0
            elif 'FUSION' in ident_upper: type_priority = 1
            else: type_priority = 2 
            ref_id = identity.replace('Var|', 'Ref|', 1).strip()
            ref_s = ref_fasta.get(ref_id, "")
            novelty_priority = 0 if "AAAAA" in ref_s else 1
            return (type_priority, novelty_priority, identity)

        df['sort_key'] = df['Identity'].apply(get_custom_sort_key)
        temp_sort_df = pd.DataFrame(df['sort_key'].tolist(), index=df.index)
        df['p_type'] = temp_sort_df[0]
        df['p_novel'] = temp_sort_df[1]
        df = df.sort_values(['p_type', 'p_novel', 'Identity'])

    for identity, group in df.groupby('Identity', sort=False):
        identity_clean = identity.strip()
        cleaned_bindlevel = group['BindLevel'].astype(str).str.strip()
        mhc_counts_sb = group[cleaned_bindlevel == 'SB']['MHC'].nunique()
        mhc_counts_sbwb = group[cleaned_bindlevel.isin(['SB', 'WB'])]['MHC'].nunique()

        var_s = var_fasta.get(identity_clean, var_fasta.get(identity_clean.split('|')[0], ""))
        ref_id = identity_clean.replace('Var|', 'Ref|', 1)
        ref_s = ref_fasta.get(ref_id, ref_fasta.get(ref_id.split('|')[0], ""))
        
        is_fusion_or_lstv = any(kw in identity_clean.upper() for kw in ['FUSION', 'LSTV'])
        
        # 判断是否在另一种 MHC 中存在
        in_other_bool = False
        if not other_df_pass.empty and identity in set(other_df_pass['Identity']):
            in_other_bool = True

        candidate_seqs =[]
        

        if var_s:
            hgvs_p_val = group.iloc[0]['hgvs_p']
            is_fs = "AAAAA" in ref_s if ref_s else False
            
            if (not is_mhci) or (is_mhci and in_other_bool):
                candidate_seqs =[var_s]
            elif is_fusion_or_lstv or is_fs or hgvs_p_val == '-':
                # [修复逻辑]：如果它是融合、移码(AAAAA)或者HGVS为空(-)，均判定为非点突变，直接保留原序列，不截断！
                candidate_seqs = [var_s]
            else:
                # [修复逻辑]：只有确认为“点突变”时，才以突变为中心，左右各保留12aa进行截短
                mut_pt = find_mismatch(ref_s, var_s) if ref_s and var_s else 0
                start_idx = max(0, mut_pt - 12)
                end_idx = min(len(var_s), mut_pt + 13)
                candidate_seqs = [var_s[start_idx : end_idx]]

        for _, row in group.iterrows():
            rows.append({
                'Gene_name': row['gene_name'], 'Uniprot_id': row['uniprot_id'], 'HGVS_p': row['hgvs_p'],
                'MHC': row['MHC'], 'Wildtype_peptide': row['Wildtype_peptide'], 'Mut_Peptide': row['Peptide'],
                'Core': row['Core'], 'BindingLevel': row['BindLevel'],
                'MHC_counts_SB': mhc_counts_sb, 'MHC_counts_SBWB': mhc_counts_sbwb,
                'Visual_Meta': (var_s, row['Peptide'], row['BindLevel']),
                'In_Other_Flag': "YES" if in_other_bool else "NO",  # 存入内部通用变量
                'Candidate_Peptides': "\n".join(candidate_seqs),
                'Identity': identity, 'is_empty': False
            })
        rows.append({'is_empty': True})
    return rows

def write_xlsx_with_visual(data, filename, is_mhci):
    workbook = xlsxwriter.Workbook(filename)
    worksheet = workbook.add_worksheet('Neoantigens')
    
    header_fmt = workbook.add_format({'bold': True, 'bg_color': '#D9EAD3', 'border': 1})
    gray_fmt = workbook.add_format({'bg_color': '#D3D3D3'})
    text_wrap = workbook.add_format({'text_wrap': True, 'valign': 'vcenter'})
    sb_font = workbook.add_format({'color': '#00B050', 'bold': True}) 
    wb_font = workbook.add_format({'color': '#FF0000', 'bold': True}) 
    
    # 修改点 1: 列名的动态变化
    other_mhc_header = 'In_MHCII' if is_mhci else 'In_MHCI'
    
    headers =['No.', 'Gene_name', 'Uniprot_id', 'HGVS_p', 'MHC', 'MHC_counts_SB', 'MHC_counts_SBWB', 'Wildtype_peptide', 'Mut_Peptide', 'Core', 
               'BindingLevel', 'Epitope_Position_Visual', other_mhc_header, '候选疫苗多肽序列']
    
    for col, h in enumerate(headers):
        worksheet.write(0, col, h, header_fmt)
    
    row_idx = 1
    no_idx = 1
    current_identity = None
    merge_start_row = 1

    for item in data:
        if item.get('is_empty'):
            worksheet.set_row(row_idx, 15, gray_fmt)
            if current_identity is not None and row_idx > merge_start_row:
                if row_idx - merge_start_row > 1:
                    worksheet.merge_range(merge_start_row, 13, row_idx-1, 13, merge_val, text_wrap)
                else:
                    worksheet.write(merge_start_row, 13, merge_val, text_wrap)
            current_identity = None
            row_idx += 1
            continue

        worksheet.write(row_idx, 0, no_idx)
        worksheet.write(row_idx, 1, item['Gene_name'])
        worksheet.write(row_idx, 2, item['Uniprot_id'])
        worksheet.write(row_idx, 3, item['HGVS_p'])
        worksheet.write(row_idx, 4, item['MHC'])
        worksheet.write(row_idx, 5, item['MHC_counts_SB'])
        worksheet.write(row_idx, 6, item['MHC_counts_SBWB'])
        worksheet.write(row_idx, 7, item['Wildtype_peptide'])
        worksheet.write(row_idx, 8, item['Mut_Peptide'])
        worksheet.write(row_idx, 9, item['Core'])
        worksheet.write(row_idx, 10, item['BindingLevel'])
        worksheet.write(row_idx, 12, item['In_Other_Flag']) # 取通用变量写入对应的动态列名
            
        full, pep, lvl = item['Visual_Meta']
        p_start = full.find(pep)
        if p_start != -1:
            font = sb_font if lvl == 'SB' else wb_font
            rich_args =[]
            if p_start > 0: rich_args.append(full[:p_start])
            rich_args.extend([font, pep])
            if p_start + len(pep) < len(full): rich_args.append(full[p_start+len(pep):])
            
            # [修复逻辑]：如果序列全是表位本身，rich_args只有2个元素，此时改用普通写入带格式的文本
            if len(rich_args) == 2:
                worksheet.write_string(row_idx, 11, rich_args[1], rich_args[0])
            else:
                worksheet.write_rich_string(row_idx, 11, *rich_args)
        else:
            worksheet.write(row_idx, 11, full)

        if current_identity is None:
            current_identity = item['Identity']
            merge_start_row = row_idx
            merge_val = item['Candidate_Peptides']
        
        row_idx += 1
        no_idx += 1

    worksheet.freeze_panes(1, 0)
    worksheet.autofilter(0, 0, row_idx-1, len(headers)-1)
    worksheet.set_column(1, 1, 30)
    worksheet.set_column(5, 6, 18)
    worksheet.set_column(11, 11, 50)
    worksheet.set_column(13, 13, 50)
    workbook.close()

# ==========================================
# 核心类
# ==========================================

class NetMHCSummarizer:
    def __init__(self, data_dir, output_dir, prefix, sample_name, mhc_genotype):
        self.data_dir = data_dir
        self.output_dir = output_dir
        self.prefix = prefix
        self.sample_name = sample_name
        self.mhc_genotype = mhc_genotype
        self.tmr_data = None
        self.ref_netmhcpan_df = None
        self.ref_netmhciipan_df = None
        self.netmhcpan_df = None
        self.netmhciipan_df = None
        self.homo_netmhcpan_df = None
        self.homo_netmhcpan_faa_file = None
        self.homo_netmhciipan_df = None
        self.homo_netmhciipan_faa_file = None
        self.id_transformation_dict = {}
        self.all_format_input = True
        self.var_sequences = {}
        self.ref_sequences = {}

    def load_source_sequences(self):
        var_faa_path = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_var_seq.faa')
        ref_faa_path = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_ref_seq.faa')

        if os.path.exists(var_faa_path):
            for record in SeqIO.parse(var_faa_path, "fasta"):
                self.var_sequences[record.id.strip()] = str(record.seq).strip()
        
        if os.path.exists(ref_faa_path):
            for record in SeqIO.parse(ref_faa_path, "fasta"):
                self.ref_sequences[record.id.strip()] = str(record.seq).strip()

    def get_wildtype_peptide_from_source(self, identity, peptide, pos):
        var_seq = self.var_sequences.get(identity)
        ref_identity = identity.replace('Var', 'Ref')
        ref_seq = self.ref_sequences.get(ref_identity)

        if not var_seq or not ref_seq: return '-'
        search_start = max(0, pos - 5)
        
        try: exact_pos = var_seq.index(peptide, search_start)
        except ValueError:
            try: exact_pos = var_seq.index(peptide)
            except ValueError: return '-'
        
        peptide_len = len(peptide)
        if exact_pos + peptide_len > len(ref_seq): return '-'
            
        wildtype_peptide = ref_seq[exact_pos : exact_pos + peptide_len]
        if wildtype_peptide == peptide: return '-'
        return wildtype_peptide

    def get_transmembrane_region(self, peptide_id, pos, peptide_length):
        midle_pos = pos + (peptide_length // 2)
        tmr_region = "-"
        if not self.tmr_data: return tmr_region
            
        in_record_block = False
        for line in self.tmr_data:
            if line.startswith('##'): continue
            if line.startswith(f'# {peptide_id} '):
                in_record_block = True
                continue
            if line.startswith('//'):
                in_record_block = False
                continue
            
            if in_record_block:
                if line.startswith('#'): continue
                fields = line.strip().split()
                if len(fields) >= 4:
                    try:
                        seq_id = fields[0]
                        if seq_id != peptide_id: continue 
                        region_type = fields[1]
                        tmr_start, tmr_end = int(fields[2]), int(fields[3])
                        
                        if tmr_start <= midle_pos + 1 <= tmr_end:
                            tmr_region = region_type 
                            break 
                    except (ValueError, IndexError):
                        continue
        return tmr_region

    def get_homologous_peptide(self, identity, peptide, homo_df, homo_faa_list):
        identity = identity.replace('Var', 'Ref')
        homo_peptides = homo_df.loc[homo_df['Identity'] == identity, 'Peptide'].tolist()
        if len(homo_peptides) == 0: return 'N', 'N', 'N'
        similarity_scores =[SequenceMatcher(None, peptide, homo_peptide).ratio() for homo_peptide in homo_peptides]
        max_similarity_idx = similarity_scores.index(max(similarity_scores))
        homo_peptide = homo_peptides[max_similarity_idx]
        for pep_id, homo_seq, homo_id in homo_faa_list:
            if identity == pep_id and homo_peptide == homo_seq:
                return 'Y', homo_peptide, homo_id
        return 'N', 'N', 'N'

    def get_affinity(self, wildtype_peptide, homo_peptide, identity, mhc_genotype, ref_df, keyword):
        identity = identity.replace('Var', 'Ref')
        competitor_peptide = wildtype_peptide if wildtype_peptide != '-' else homo_peptide
        
        if competitor_peptide not in ['-', 'N']:
            affinity_values = ref_df.loc[
                (ref_df['Peptide'] == competitor_peptide) & 
                (ref_df['Identity'] == identity) & 
                (ref_df['MHC'] == mhc_genotype), keyword
            ].values
            if len(affinity_values) > 0: return affinity_values[0]
        return '-'

    def load_id_transformation(self):
        id_transformation_file = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_id_transformation.txt')
        if not os.path.exists(id_transformation_file): return

        with open(id_transformation_file, 'r') as f:
            for line in f.readlines():
                if '->' in line:
                    parts = line.strip().split('->')
                    origin = parts[0].replace('>', '').strip()
                    trans = parts[1].split()[0].replace('>', '').strip()
                    
                    ident_items = origin.split('|')
                    if len(ident_items) >= 4 and ident_items[0] == 'Var': 
                        self.id_transformation_dict[trans] =[origin, ident_items[2], ident_items[1], ident_items[3]]
                    else:
                        self.all_format_input = False
                        self.id_transformation_dict[trans] =[origin]
        
    def load_tmr_data(self):
        tmr_file = os.path.join(self.data_dir, '04.TransMembrane.DeepTMHMM', 'TMRs.gff3')
        if not os.path.exists(tmr_file): 
            self.tmr_data = None
            return
        with open(tmr_file, 'r') as f:
            self.tmr_data = f.readlines()

    def load_homo_faa(self, homo_netmhc_faa_file):
        homo_netmhc_faa_list =[]
        if not os.path.exists(homo_netmhc_faa_file): return homo_netmhc_faa_list
        for record in SeqIO.parse(homo_netmhc_faa_file, 'fasta'):
            homo_netmhc_faa_list.append([record.id, str(record.seq), record.description.split("_")[-1]])
        return homo_netmhc_faa_list

    def load_dataframes(self):
        if self.mhc_genotype in ['mhci', 'all']:
            csv_path = os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCpan.csv')
            if os.path.exists(csv_path): self.netmhcpan_df = pd.read_csv(csv_path)

            ref_path = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCpan.csv')
            if os.path.exists(ref_path): self.ref_netmhcpan_df = pd.read_csv(ref_path, low_memory=False)

            homo_csv = os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCpan.csv')
            if os.path.exists(homo_csv): self.homo_netmhcpan_df = pd.read_csv(homo_csv)
            self.homo_netmhcpan_faa_list = self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCpan_homologous.faa'))

        if self.mhc_genotype in['mhcii', 'all']:
            csv_path = os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCIIpan.csv')
            if os.path.exists(csv_path): self.netmhciipan_df = pd.read_csv(csv_path)
            
            ref_path = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCIIpan.csv')
            if os.path.exists(ref_path): self.ref_netmhciipan_df = pd.read_csv(ref_path, low_memory=False)
            
            homo_csv = os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCIIpan.csv')
            if os.path.exists(homo_csv): self.homo_netmhciipan_df = pd.read_csv(homo_csv)
            self.homo_netmhciipan_faa_list = self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCIIpan_homologous.faa'))

    def summarize(self, args):
        self.load_id_transformation()
        self.load_tmr_data()
        self.load_dataframes()
        self.load_source_sequences() 

        def process_df(df, ref_df, homo_df, homo_list, is_mhcii=False):
            if df is None: return None
            df['TransMemb'] = df.apply(lambda row: self.get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide'])), axis=1)
            df['InCutmerRate'] = '-'
            df['InCutmerRegion'] = '-'
            df['Wildtype_peptide'] = df.apply(lambda row: self.get_wildtype_peptide_from_source(row['Identity'], row['Peptide'], row['Pos']), axis=1)
            df['HomoExsit'] = 'N'
            df['Homo_id'] = '-'
            df['Homo_peptide'] = '-'
            
            if homo_df is not None:
                homo_df['Identity'] = homo_df['Identity'].str.split('_').str[0]
                results = df.loc[df['Wildtype_peptide'] == '-',['Peptide', 'Identity']].apply(
                    lambda row: self.get_homologous_peptide(row['Identity'], row['Peptide'], homo_df, homo_list), axis=1, result_type='expand')
                if not results.empty:
                    results.columns =['HomoExsit', 'Homo_peptide', 'Homo_id']
                    df.update(results)

            if ref_df is not None:
                kw = 'Affinity(nM)' if is_mhcii else 'Aff(nM)'
                df['Aff(nM)_competitor'] = df.apply(lambda row: self.get_affinity(row['Wildtype_peptide'], row['Homo_peptide'], row['Identity'], row["MHC"], ref_df, kw), axis=1)
                df['Aff(nM)_competitor/Aff(nM)'] = df.apply(
                    lambda row: round(float(row['Aff(nM)_competitor']) / float(row[kw]), 2) if pd.notna(row['Aff(nM)_competitor']) and str(row['Aff(nM)_competitor']) != '-' and pd.notna(row[kw]) else '-', axis=1)
            return df

        self.netmhcpan_df = process_df(self.netmhcpan_df, self.ref_netmhcpan_df, self.homo_netmhcpan_df, self.homo_netmhcpan_faa_list)
        self.netmhciipan_df = process_df(self.netmhciipan_df, self.ref_netmhciipan_df, self.homo_netmhciipan_df, self.homo_netmhciipan_faa_list, is_mhcii=True)

        def restore_identity(df):
            if df is None: return None
            df['Identity_raw'] = df['Identity']
            if self.all_format_input:
                df['Gene_name'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x.strip(), ['-','-','-','-'])[1])
                df['Transcript_id'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x.strip(), ['-','-','-','-'])[2])
                df['HGVS_p'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x.strip(), ['-','-','-','-'])[3])
            
            df['Identity'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x.strip(), [x.strip()])[0])
            df.drop(columns=['Identity_raw'], inplace=True)
            
            cols = list(df.columns)
            if self.all_format_input:
                new_cols =['Gene_name', 'Transcript_id', 'HGVS_p', 'Identity']
                ordered_cols = new_cols +[c for c in cols if c not in new_cols]
                return df[ordered_cols]
            return df
        
        self.netmhcpan_df = restore_identity(self.netmhcpan_df)
        self.netmhciipan_df = restore_identity(self.netmhciipan_df)
        
        if self.netmhcpan_df is not None: self.netmhcpan_df = self.netmhcpan_df[self.netmhcpan_df['Wildtype_peptide'] != '-']
        if self.netmhciipan_df is not None: self.netmhciipan_df = self.netmhciipan_df[self.netmhciipan_df['Wildtype_peptide'] != '-']
            
        if self.netmhcpan_df is not None:
            self.netmhcpan_df.drop(columns=[c for c in['Of', 'Gp', 'Gl', 'Ip', 'Il', 'Pos'] if c in self.netmhcpan_df.columns], inplace=True)
        if self.netmhciipan_df is not None:
            self.netmhciipan_df.drop(columns=[c for c in['Of', 'Exp_Bind', 'Pos', 'Core_Rel'] if c in self.netmhciipan_df.columns], inplace=True)

        if os.path.exists(self.output_dir): shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

        if self.netmhcpan_df is not None and not self.netmhcpan_df.empty:
            self.netmhcpan_df.to_csv(os.path.join(self.output_dir, f'{self.prefix}_netMHCpan_deliverable.csv'), index=False)
        if self.netmhciipan_df is not None and not self.netmhciipan_df.empty:
            self.netmhciipan_df.to_csv(os.path.join(self.output_dir, f'{self.prefix}_netMHCIIpan_deliverable.csv'), index=False)

        self.generate_excel_deliverables(args)


    def generate_excel_deliverables(self, args):
        whitelist = set()
        if args.include_list:
            if os.path.exists(args.include_list):
                with open(args.include_list, 'r') as f: whitelist = set(line.strip() for line in f if line.strip())
            else:
                items = [x.strip() for x in args.include_list.split(',') if x.strip()]
                if items: whitelist = set(items)

        gene_order = None
        if args.genelist and os.path.exists(args.genelist):
            gene_order = pd.read_csv(args.genelist, sep='\t')['gene_name'].tolist()

        var_fasta = {}
        ref_fasta = {}
        for short_id, seq in self.var_sequences.items():
            full_id = self.id_transformation_dict.get(short_id, [short_id])[0]
            var_fasta[full_id.strip()] = seq.strip()
            var_fasta[short_id.strip()] = seq.strip() 
            
        for short_id, seq in self.ref_sequences.items():
            var_short_id = short_id.replace('Ref', 'Var')
            var_full_id = self.id_transformation_dict.get(var_short_id, [var_short_id])[0]
            ref_full_id = var_full_id.replace('Var', 'Ref', 1)
            ref_fasta[ref_full_id.strip()] = seq.strip()
            ref_fasta[short_id.strip()] = seq.strip()

        tm_drop_identities = set()
        all_identities = set()
        if self.netmhcpan_df is not None and not self.netmhcpan_df.empty:
            all_identities.update(self.netmhcpan_df['Identity'].dropna().unique())
        if self.netmhciipan_df is not None and not self.netmhciipan_df.empty:
            all_identities.update(self.netmhciipan_df['Identity'].dropna().unique())
            
        for ident in all_identities:
            ident_clean = ident.strip()
            var_s = var_fasta.get(ident_clean, var_fasta.get(ident_clean.split('|')[0], ""))
            ref_id = ident_clean.replace('Var|', 'Ref|', 1)
            ref_s = ref_fasta.get(ref_id, ref_fasta.get(ref_id.split('|')[0], ""))
            
            is_point_mutation = (len(var_s) <= 37) and ("AAAAA" not in ref_s)
            has_tm = False
            if self.netmhcpan_df is not None and not self.netmhcpan_df.empty:
                tm_hits = self.netmhcpan_df[(self.netmhcpan_df['Identity'] == ident) & (self.netmhcpan_df['TransMemb'] == 'TMhelix')]
                if not tm_hits.empty: has_tm = True
                
            if not has_tm and self.netmhciipan_df is not None and not self.netmhciipan_df.empty:
                tm_hits = self.netmhciipan_df[(self.netmhciipan_df['Identity'] == ident) & (self.netmhciipan_df['TransMemb'] == 'TMhelix')]
                if not tm_hits.empty: has_tm = True
                
            if is_point_mutation and has_tm:
                tm_drop_identities.add(ident)

        m1_pass = pd.DataFrame()
        m2_pass = pd.DataFrame()

        if self.netmhcpan_df is not None and not self.netmhcpan_df.empty:
            m1_pass = filter_data(self.netmhcpan_df, 'Aff(nM)', '%Rank_EL', args.mhc1_aff, args.mhc1_rank, args.bind_levels, tm_drop_identities, whitelist)
        
        if self.netmhciipan_df is not None and not self.netmhciipan_df.empty:
            m2_pass = filter_data(self.netmhciipan_df, 'Affinity(nM)', '%Rank_EL', args.mhc2_aff, args.mhc2_rank, args.bind_levels, tm_drop_identities, whitelist)

        if m1_pass.empty and 'Identity' not in m1_pass.columns: m1_pass['Identity'] =[]
        if m2_pass.empty and 'Identity' not in m2_pass.columns: m2_pass['Identity'] =[]

        # 修改点 3: 显式传入 is_mhci 的状态到处理与写入函数
        if not m1_pass.empty:
            res1 = prepare_data(m1_pass, var_fasta, ref_fasta, m2_pass, is_mhci=True, gene_order=gene_order)
            write_xlsx_with_visual(res1, os.path.join(self.output_dir, f"{self.prefix}_netMHCI_deliverable.xlsx"), is_mhci=True)
            
        if not m2_pass.empty:
            res2 = prepare_data(m2_pass, var_fasta, ref_fasta, m1_pass, is_mhci=False, gene_order=gene_order)
            write_xlsx_with_visual(res2, os.path.join(self.output_dir, f"{self.prefix}_netMHCII_deliverable.xlsx"), is_mhci=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize netMHC results and generate formatted deliverables.')
    parser.add_argument('-d', '--data_dir', type=str, required=True, help='Path to the data directory')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to the output directory')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix of the output result')
    parser.add_argument('-n', '--sample_name', type=str, required=True, help='Sample name of this sample')
    parser.add_argument('-m', '--mhc_genotype', choices=['mhci', 'mhcii', 'all'], required=True, help='MHC genotypes')
    
    parser.add_argument('--mhc1_aff', type=float, default=500.0, help='MHC-I Affinity filter limit')
    parser.add_argument('--mhc1_rank', type=float, default=2.0, help='MHC-I Rank filter limit')
    parser.add_argument('--mhc2_aff', type=float, default=1000.0, help='MHC-II Affinity filter limit')
    parser.add_argument('--mhc2_rank', type=float, default=10.0, help='MHC-II Rank filter limit')
    parser.add_argument('--bind_levels', type=str, default='SB,WB', help='Bind levels to keep (e.g., SB,WB)')
    parser.add_argument('--genelist', required=False, help='Gene list file for sorting (Optional)')
    parser.add_argument('--include_list', required=False, help='File path OR comma-separated string of Identities to force include')
    
    args = parser.parse_args()
    summarizer = NetMHCSummarizer(args.data_dir, args.output_dir, args.prefix, args.sample_name, args.mhc_genotype)
    summarizer.summarize(args)
