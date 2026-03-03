import os
import re
import shutil
import argparse
import pandas as pd
from Bio import SeqIO
from difflib import SequenceMatcher

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
        """Loads the full variant and reference protein sequences from FASTA files."""
        var_faa_path = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_var_seq.faa')
        ref_faa_path = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_ref_seq.faa')

        if os.path.exists(var_faa_path):
            for record in SeqIO.parse(var_faa_path, "fasta"):
                self.var_sequences[record.id] = str(record.seq)
        
        if os.path.exists(ref_faa_path):
            for record in SeqIO.parse(ref_faa_path, "fasta"):
                self.ref_sequences[record.id] = str(record.seq)

    def get_wildtype_peptide_from_source(self, identity, peptide, pos):
        """
        Directly extracts the corresponding wild-type peptide from the full
        reference protein sequence based on the neoantigen's position.
        """
        var_seq = self.var_sequences.get(identity)
        ref_identity = identity.replace('Var', 'Ref')
        ref_seq = self.ref_sequences.get(ref_identity)

        if not var_seq or not ref_seq:
            return '-'

        search_start = max(0, pos - 5)
        
        try:
            exact_pos = var_seq.index(peptide, search_start)
        except ValueError:
            try:
                exact_pos = var_seq.index(peptide)
            except ValueError:
                return '-'
        
        peptide_len = len(peptide)
        if exact_pos + peptide_len > len(ref_seq):
            return '-'
            
        wildtype_peptide = ref_seq[exact_pos : exact_pos + peptide_len]
        
        if wildtype_peptide == peptide:
            return '-'

        return wildtype_peptide

    def get_transmembrane_region(self, peptide_id, pos, peptide_length):
        """
        Determines if a peptide is within a transmembrane region based on DeepTMHMM GFF3 output.
        """
        midle_pos = pos + (peptide_length // 2)
        tmr_region = "-"
        if not self.tmr_data: 
            return tmr_region
            
        in_record_block = False
        for line in self.tmr_data:
            if line.startswith('##'): continue
            if line.startswith(f'# {peptide_id}'):
                in_record_block = True
                continue
            if line.startswith('//'):
                in_record_block = False
                continue
            
            if in_record_block:
                # GFF3 feature lines should not start with '#', but we check for safety
                if line.startswith('#'): continue
                
                fields = line.strip().split('\t')
                
                # *** FIX STARTS HERE ***
                # Add a defensive check to ensure the line has enough columns for a GFF feature
                if len(fields) >= 5:
                    try:
                        # GFF3 is 1-based, columns are seqid, source, type, start, end
                        # fields[0] is seqid, fields[2] is type, fields[3] is start, fields[4] is end
                        seq_id = fields[0]
                        if seq_id != peptide_id: continue # Ensure we are still in the correct protein block

                        tmr_start, tmr_end = int(fields[3]), int(fields[4])
                        
                        # GFF3 coordinates are 1-based and inclusive. pos is 0-based.
                        # So we check if (pos+1) is within [start, end]
                        if tmr_start <= midle_pos + 1 <= tmr_end:
                            tmr_region = fields[2] # The 'type' of region (e.g., TRANSMEM)
                            break # Found the region, no need to check further for this peptide
                    except (ValueError, IndexError):
                        # If conversion to int fails or another unexpected index error occurs, skip this malformed line
                        continue
                # *** FIX ENDS HERE ***
        return tmr_region

    def get_homologous_peptide(self, identity, peptide, homo_df, homo_faa_list):
        identity = identity.replace('Var', 'Ref')
        homo_peptides = homo_df.loc[homo_df['Identity'] == identity, 'Peptide'].tolist()
        if len(homo_peptides) == 0:
            return 'N', 'N', 'N'
        similarity_scores = [SequenceMatcher(None, peptide, homo_peptide).ratio() for homo_peptide in homo_peptides]
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
            if len(affinity_values) > 0:
                return affinity_values[0]
        
        return '-'

    def load_id_transformation(self):
        id_transformation_file = os.path.join(self.data_dir, '01.protein_sequence', f'{self.sample_name}_id_transformation.txt')
        if not os.path.exists(id_transformation_file): return

        with open(id_transformation_file, 'r') as f:
            for line in f.readlines():
                if line.startswith(">Var") and '->' in line:
                    parts = line.strip().split('->')
                    identity_origin = parts[0].lstrip('>')
                    identity_transformated = parts[1].split()[0].lstrip('>')
                    
                    identity_origin_item = identity_origin.split('|')
                    if len(identity_origin_item) == 5  and identity_origin_item[0] == 'Var' : 
                        transcript_id = identity_origin_item[1]
                        gene_name = identity_origin_item[2]
                        hgvs_p = identity_origin_item[3]
                        self.id_transformation_dict[identity_transformated] = [identity_origin, gene_name, transcript_id, hgvs_p]
                    else:
                        self.all_format_input = False
                        self.id_transformation_dict[identity_transformated] = [identity_origin]
        
    def load_tmr_data(self):
        tmr_file = os.path.join(self.data_dir, '04.TransMembrane.DeepTMHMM', 'TMRs.gff3')
        if not os.path.exists(tmr_file): 
            self.tmr_data = None
            return
        with open(tmr_file, 'r') as f:
            self.tmr_data = f.readlines()

    def load_homo_faa(self, homo_netmhc_faa_file):
        homo_netmhc_faa_list = []
        if not os.path.exists(homo_netmhc_faa_file): return homo_netmhc_faa_list
        for record in SeqIO.parse(homo_netmhc_faa_file, 'fasta'):
            pep_id = record.id
            pep_description = record.description.split("_")[-1]
            pep_seq = str(record.seq)
            homo_netmhc_faa_list.append([pep_id, pep_seq, pep_description])
        return homo_netmhc_faa_list

    def load_dataframes(self):
        if self.mhc_genotype in ['mhci', 'all']:
            netmhcpan_csv = os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCpan.csv')
            if os.path.exists(netmhcpan_csv): self.netmhcpan_df = pd.read_csv(netmhcpan_csv)

            ref_netmhcpan_file = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCpan.csv')
            if os.path.exists(ref_netmhcpan_file): self.ref_netmhcpan_df = pd.read_csv(ref_netmhcpan_file, low_memory=False) # FIX: Add low_memory=False

            homo_csv = os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCpan.csv')
            if os.path.exists(homo_csv): self.homo_netmhcpan_df = pd.read_csv(homo_csv)
            self.homo_netmhcpan_faa_list = self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCpan_homologous.faa'))

        if self.mhc_genotype in ['mhcii', 'all']:
            netmhciipan_csv = os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCIIpan.csv')
            if os.path.exists(netmhciipan_csv): self.netmhciipan_df = pd.read_csv(netmhciipan_csv)
            
            ref_netmhciipan_file = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCIIpan.csv')
            if os.path.exists(ref_netmhciipan_file): self.ref_netmhciipan_df = pd.read_csv(ref_netmhciipan_file, low_memory=False) # FIX: Add low_memory=False
            
            homo_csv = os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCIIpan.csv')
            if os.path.exists(homo_csv): self.homo_netmhciipan_df = pd.read_csv(homo_csv)
            self.homo_netmhciipan_faa_list = self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCIIpan_homologous.faa'))

    def summarize(self):
        self.load_id_transformation()
        self.load_tmr_data()
        self.load_dataframes()
        self.load_source_sequences() 

        if self.mhc_genotype in ['mhci', 'all'] and self.netmhcpan_df is not None:
            self.netmhcpan_df['TransMemb'] = self.netmhcpan_df.apply(lambda row: self.get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide'])), axis=1)
            self.netmhcpan_df['InCutmerRate'] = '-'
            self.netmhcpan_df['InCutmerRegion'] = '-'
            self.netmhcpan_df['Wildtype_peptide'] = self.netmhcpan_df.apply(
                lambda row: self.get_wildtype_peptide_from_source(row['Identity'], row['Peptide'], row['Pos']), axis=1
            )
            self.netmhcpan_df['HomoExsit'] = 'N'
            self.netmhcpan_df['Homo_id'] = '-'
            self.netmhcpan_df['Homo_peptide'] = '-'
            
            if self.homo_netmhcpan_df is not None:
                self.homo_netmhcpan_df['Identity'] = self.homo_netmhcpan_df['Identity'].str.split('_').str[0]
                results = self.netmhcpan_df.loc[self.netmhcpan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(
                    lambda row: self.get_homologous_peptide(row['Identity'], row['Peptide'], self.homo_netmhcpan_df, self.homo_netmhcpan_faa_list), 
                    axis=1, result_type='expand')
                if not results.empty:
                    results.columns = ['HomoExsit', 'Homo_peptide', 'Homo_id']
                    self.netmhcpan_df.update(results)

            if self.ref_netmhcpan_df is not None:
                self.netmhcpan_df['Aff(nM)_competitor'] = self.netmhcpan_df.apply(lambda row: self.get_affinity(row['Wildtype_peptide'], row['Homo_peptide'], row['Identity'], row["MHC"], self.ref_netmhcpan_df, 'Aff(nM)'), axis=1)
                self.netmhcpan_df['Aff(nM)_competitor/Aff(nM)'] = self.netmhcpan_df.apply(
                    lambda row: round(float(row['Aff(nM)_competitor']) / float(row['Aff(nM)']), 2) if pd.notna(row['Aff(nM)_competitor']) and str(row['Aff(nM)_competitor']) != '-' and pd.notna(row['Aff(nM)']) else '-', axis=1
                )
        
        if self.mhc_genotype in ['mhcii', 'all'] and self.netmhciipan_df is not None:
            self.netmhciipan_df['TransMemb'] = self.netmhciipan_df.apply(lambda row: self.get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide'])), axis=1)
            self.netmhciipan_df['InCutmerRate'] = '-'
            self.netmhciipan_df['InCutmerRegion'] = '-'
            self.netmhciipan_df['Wildtype_peptide'] = self.netmhciipan_df.apply(
                lambda row: self.get_wildtype_peptide_from_source(row['Identity'], row['Peptide'], row['Pos']), axis=1
            )
            self.netmhciipan_df['HomoExsit'] = 'N'
            self.netmhciipan_df['Homo_id'] = '-'
            self.netmhciipan_df['Homo_peptide'] = '-'
            
            if self.homo_netmhciipan_df is not None:
                self.homo_netmhciipan_df['Identity'] = self.homo_netmhciipan_df['Identity'].str.split('_').str[0]
                results = self.netmhciipan_df.loc[self.netmhciipan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(
                    lambda row: self.get_homologous_peptide(row['Identity'], row['Peptide'], self.homo_netmhciipan_df, self.homo_netmhciipan_faa_list), 
                    axis=1, result_type='expand'
                )
                if not results.empty:
                    results.columns = ['HomoExsit', 'Homo_peptide', 'Homo_id']
                    self.netmhciipan_df.update(results)
            
            if self.ref_netmhciipan_df is not None:
                self.netmhciipan_df['Aff(nM)_competitor'] = self.netmhciipan_df.apply(lambda row: self.get_affinity(row['Wildtype_peptide'], row['Homo_peptide'], row['Identity'], row["MHC"], self.ref_netmhciipan_df, 'Affinity(nM)'), axis=1)
                self.netmhciipan_df['Aff(nM)_competitor/Affinity(nM)'] = self.netmhciipan_df.apply(
                    lambda row: round(float(row['Aff(nM)_competitor']) / float(row['Affinity(nM)']), 2) if pd.notna(row['Aff(nM)_competitor']) and str(row['Aff(nM)_competitor']) != '-' and pd.notna(row['Affinity(nM)']) else '-', axis=1
                )

        def restore_identity(df):
            if df is None: return None
            df['Identity_raw'] = df['Identity']
            
            if self.all_format_input:
                df['Gene_name'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x, ['-','-','-','-'])[1])
                df['Transcript_id'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x, ['-','-','-','-'])[2])
                df['HGVS_p'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x, ['-','-','-','-'])[3])
            
            df['Identity'] = df['Identity_raw'].apply(lambda x: self.id_transformation_dict.get(x, [x])[0])
            
            df.drop(columns=['Identity_raw'], inplace=True)
            
            cols = list(df.columns)
            if self.all_format_input:
                new_cols = ['Gene_name', 'Transcript_id', 'HGVS_p', 'Identity']
                ordered_cols = new_cols + [c for c in cols if c not in new_cols]
                return df[ordered_cols]
            return df
        
        self.netmhcpan_df = restore_identity(self.netmhcpan_df)
        self.netmhciipan_df = restore_identity(self.netmhciipan_df)
        
        if self.netmhcpan_df is not None:
            self.netmhcpan_df = self.netmhcpan_df[self.netmhcpan_df['Wildtype_peptide'] != '-']
        if self.netmhciipan_df is not None:
            self.netmhciipan_df = self.netmhciipan_df[self.netmhciipan_df['Wildtype_peptide'] != '-']
            
        if self.netmhcpan_df is not None:
            netmhcpan_columns_to_drop = ['Of', 'Gp', 'Gl', 'Ip', 'Il', 'Pos']
            self.netmhcpan_df.drop(columns=[col for col in netmhcpan_columns_to_drop if col in self.netmhcpan_df.columns], inplace=True)
        if self.netmhciipan_df is not None:
            netmhciipan_columns_to_drop = ['Of', 'Exp_Bind', 'Pos', 'Core_Rel']
            self.netmhciipan_df.drop(columns=[col for col in netmhciipan_columns_to_drop if col in self.netmhciipan_df.columns], inplace=True)

        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

        if self.mhc_genotype in ['mhci', 'all'] and self.netmhcpan_df is not None and not self.netmhcpan_df.empty:
            self.netmhcpan_df.to_csv(os.path.join(self.output_dir, f'{self.prefix}_netMHCpan_deliverable.csv'), index=False)
        if self.mhc_genotype in ['mhcii', 'all'] and self.netmhciipan_df is not None and not self.netmhciipan_df.empty:
            self.netmhciipan_df.to_csv(os.path.join(self.output_dir, f'{self.prefix}_netMHCIIpan_deliverable.csv'), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize netMHC results.')
    parser.add_argument('-d', '--data_dir', type=str, required=True, help='Path to the data directory')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to the output directory')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='Prefix of the output result')
    parser.add_argument('-n', '--sample_name', type=str, required=True, help='Sample name of this sample')
    parser.add_argument('-m', '--mhc_genotype', choices=['mhci', 'mhcii', 'all'], required=True, help='MHC genotypes')
    args = parser.parse_args()

    summarizer = NetMHCSummarizer(args.data_dir, args.output_dir, args.prefix, args.sample_name, args.mhc_genotype)
    summarizer.summarize()