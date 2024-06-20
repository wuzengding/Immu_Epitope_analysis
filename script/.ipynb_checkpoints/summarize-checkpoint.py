import os
import shutil
import argparse
import pandas as pd
from difflib import SequenceMatcher

def get_transmembrane_region(peptide_id, pos, peptide_length, tmr_data):
    """
    Determine the transmembrane region of a peptide based on the provided TMR data.
    """
    midle_pos = pos + (peptide_length//2)

    tmr_region="-"
    for i, line in enumerate(tmr_data):
        if line.startswith(f'# {peptide_id}'):
            for j in range(i+1, len(tmr_data)):
                if tmr_data[j].startswith('//'):
                    break
                if tmr_data[j].startswith("Var"):  ##绑定分析peptide的identiy必须是Var开头的
                    #print(line)
                    fields = tmr_data[j].strip().split('\t')
                    tmr_start, tmr_end = int(fields[2]), int(fields[3])
                    if  tmr_start <= midle_pos < tmr_end:
                        tmr_region = fields[1]
    return tmr_region

def get_wildtype_peptide(identity, peptide, ref_df):
    """
    Find the most similar peptide for the given identity in the reference file.
    """
    ref_peptides = ref_df.loc[ref_df['Identity'] == identity, 'Peptide'].tolist()

    if not ref_peptides:
        return '-'

    similarity_scores = [SequenceMatcher(None, peptide, ref_peptide).ratio() for ref_peptide in ref_peptides]
    max_similarity_idx = similarity_scores.index(max(similarity_scores))
    wildtype_peptide = ref_peptides[max_similarity_idx]

    return wildtype_peptide

def get_homologous_peptide(identity, peptide, homo_df, faa_file):
    """
    Find the most similar homologous peptide for the given identity in the homologous file.
    """
    homo_peptides = homo_df.loc[homo_df['Identity'] == identity, 'Peptide'].tolist()

    if not homo_peptides:
        return '-', '-', '-'

    similarity_scores = [SequenceMatcher(None, peptide, homo_peptide).ratio() for homo_peptide in homo_peptides]
    max_similarity_idx = similarity_scores.index(max(similarity_scores))
    homo_peptide = homo_peptides[max_similarity_idx]
    homo_id = None

    with open(faa_file, 'r') as f:
        for line in f:
            if line.startswith('>') and homo_peptide in line:
                homo_id = line.strip('>').split()[1]
                break

    if homo_id is None:
        return 'Y', homo_peptide, '-'
    else:
        return 'Y', homo_peptide, homo_id

def get_affinity(wildtype_peptide, identity, ref_df):
    """
    Get the affinity (nM) of the wildtype peptide from the reference file.
    """
    affinity = ref_df.loc[(ref_df['Peptide'] == wildtype_peptide) & (ref_df['Identity'] == identity), 'Aff(nM)'].values

    if len(affinity) > 0:
        return affinity[0]
    else:
        return '-'

def summarize(data_dir, output_dir, prefix, sample_name):
    netmhcpan_file = os.path.join(data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{sample_name}_netMHCpan.csv')
    netmhciipan_file = os.path.join(data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{sample_name}_netMHCIIpan.csv')
    ref_netmhcpan_file = os.path.join(data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{sample_name}_netMHCpan.csv')
    ref_netmhciipan_file = os.path.join(data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{sample_name}_netMHCIIpan.csv')
    homo_netmhcpan_file = os.path.join(data_dir, '03.homologous', 'parsed', f'{sample_name}_homologous_netMHCpan.csv')
    homo_netmhciipan_file = os.path.join(data_dir, '03.homologous', 'parsed', f'{sample_name}_homologous_netMHCIIpan.csv')
    homo_netmhcpan_faa_file = os.path.join(data_dir, '03.homologous', f'{sample_name}_netMHCpan_homologous.faa')
    homo_netmhciipan_faa_file = os.path.join(data_dir, '03.homologous', f'{sample_name}_netMHCIIpan_homologous.faa')
    tmr_file = os.path.join(data_dir, '04.TransMembrane.DeepTMHMM', 'TMRs.gff3')

    netmhcpan_df = pd.read_csv(netmhcpan_file)
    netmhciipan_df = pd.read_csv(netmhciipan_file)
    ref_netmhcpan_df = pd.read_csv(ref_netmhcpan_file)
    ref_netmhciipan_df = pd.read_csv(ref_netmhciipan_file)
    homo_netmhcpan_df = pd.read_csv(homo_netmhcpan_file)
    homo_netmhciipan_df = pd.read_csv(homo_netmhciipan_file)

    with open(tmr_file, 'r') as f:
        tmr_data = f.readlines()

    netmhcpan_df['TransMemb'] = netmhcpan_df.apply(lambda row: get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide']), tmr_data), axis=1)
    netmhciipan_df['TransMemb'] = netmhciipan_df.apply(lambda row: get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide']), tmr_data), axis=1)

    netmhcpan_df['InCutmerRate'] = '-'
    netmhcpan_df['InCutmerRegion'] = '-'
    netmhciipan_df['InCutmerRate'] = '-'
    netmhciipan_df['InCutmerRegion'] = '-'

    netmhcpan_df['Wildtype_peptide'] = netmhcpan_df.apply(lambda row: get_wildtype_peptide(row['Identity'], row['Peptide'], ref_netmhcpan_df), axis=1)
    netmhciipan_df['Wildtype_peptide'] = netmhciipan_df.apply(lambda row: get_wildtype_peptide(row['Identity'], row['Peptide'], ref_netmhciipan_df), axis=1)

    netmhcpan_df['HomoExsit'] = '-'
    netmhcpan_df['Homo_id'] = '-'
    netmhcpan_df['Homo_peptide'] = '-'
    netmhciipan_df['HomoExsit'] = '-'
    netmhciipan_df['Homo_id'] = '-'
    netmhciipan_df['Homo_peptide'] = '-'

    netmhcpan_df.loc[netmhcpan_df['Wildtype_peptide'] == '-', ['HomoExsit', 'Homo_peptide', 'Homo_id']] = netmhcpan_df.loc[netmhcpan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(lambda row: get_homologous_peptide(row['Identity'], row['Peptide'], homo_netmhcpan_df, homo_netmhcpan_faa_file), axis=1, result_type='expand')
    netmhciipan_df.loc[netmhciipan_df['Wildtype_peptide'] == '-', ['HomoExsit', 'Homo_peptide', 'Homo_id']] = netmhciipan_df.loc[netmhciipan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(lambda row: get_homologous_peptide(row['Identity'], row['Peptide'], homo_netmhciipan_df, homo_netmhciipan_faa_file), axis=1, result_type='expand')

    netmhcpan_df['Aff(nM)_competitor'] = netmhcpan_df.apply(lambda row: get_affinity(row['Wildtype_peptide'], row['Identity'], ref_netmhcpan_df) if row['Wildtype_peptide'] != '-' else '-', axis=1)
    netmhciipan_df['Aff(nM)_competitor'] = netmhciipan_df.apply(lambda row: get_affinity(row['Wildtype_peptide'], row['Identity'], ref_netmhciipan_df) if row['Wildtype_peptide'] != '-' else '-', axis=1)

    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    else:
        os.mkdir(output_dir)
    netmhcpan_df.to_csv(os.path.join(data_dir, f'{prefix}_netMHCpan_deliverable.csv'), index=False)
    netmhciipan_df.to_csv(os.path.join(data_dir, f'{prefix}_netMHCIIpan_deliverable.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize netMHC results.')
    parser.add_argument('-d','--data_dir', type=str, help='Path to the data directory')
    parser.add_argument('-o','--output_dir', type=str, help='Path to the output directory')
    parser.add_argument('-p','--prefix', type=str, help='Prefix of the ouput result')
    parser.add_argument('-n','--sample_name', type=str, help='Sample name of this sample ')
    args = parser.parse_args()

    summarize(args.data_dir, args.output_dir, args.prefix, args.sample_name)