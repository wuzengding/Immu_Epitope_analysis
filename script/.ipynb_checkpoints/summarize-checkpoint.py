import os
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


    def get_transmembrane_region(self, peptide_id, pos, peptide_length):
        midle_pos = pos + (peptide_length // 2)
        tmr_region = "-"
        for i, line in enumerate(self.tmr_data):
            if line.startswith(f'# {peptide_id}'):
                for j in range(i+1, len(self.tmr_data)):
                    if self.tmr_data[j].startswith('//'):
                        break
                    if self.tmr_data[j].startswith("Var"):
                        fields = self.tmr_data[j].strip().split('\t')
                        tmr_start, tmr_end = int(fields[2]), int(fields[3])
                        if tmr_start <= midle_pos < tmr_end:
                            tmr_region = fields[1]
        return tmr_region

    def get_wildtype_peptide(self, identity, peptide, ref_df):
        similarity_cutoff = 0.75
        length_diff_cutoff = 4
        """
        similarity_score of ("SKQTTEVVE","SKQTTEVVE") = 1.0
                            ("SKQTTEVVE","SKKTTEVVE") = 0.888
                            ("SKQTTEVV","SKKTTEVV")   = 0.875
                            ("SKQATEVLNEY","SKQATEVLNEV") = 0.9
                            ("SKQTTEVVE","QSKQEVVEA") = 0.777
                            ("ISGAENTRSY","SLLLSGVVAY") = 0.4
        recommed: similarity_cutoff = 0.75;
                  length_diff_cutoff = 4
        """

        # Replace 'Var' with 'Ref' in the identity
        identity = identity.replace('Var', 'Ref')
    
        # Filter ref_df to get the list of ref_peptides with the same identity
        ref_peptides = ref_df.loc[ref_df['Identity'] == identity, 'Peptide'].tolist()
        
        # Return '-' if there are no matching ref_peptides
        if not ref_peptides:
            return '-'
    
        # Function to get peptides with a specified length difference
        def get_similar_length_peptides(ref_peptides, peptide_length, length_diff):
            return [ref_peptide for ref_peptide in ref_peptides if (abs(len(ref_peptide) - peptide_length) == length_diff)]
    
        peptide_length = len(peptide)
        length_diff = 0
    
        while True:
            # Get ref_peptides with the current length difference
            similar_length_peptides = get_similar_length_peptides(ref_peptides, peptide_length, length_diff)
        
            # If we found peptides with the current length difference, calculate similarity scores
            if similar_length_peptides:
                similarity_scores = [SequenceMatcher(None, peptide, ref_peptide).ratio() for ref_peptide in similar_length_peptides]

                if max(similarity_scores) >= similarity_cutoff:
                    max_similarity_idx = similarity_scores.index(max(similarity_scores))
                    similar_length_peptides[max_similarity_idx]
                    return similar_length_peptides[max_similarity_idx]
        
            # If no peptides found, increase the length difference and try again
            length_diff += 1
        
            # If length_diff exceeds a reasonable limit, return '-'
            if (length_diff > peptide_length + max(len(r) for r in ref_peptides)) or (length_diff >= length_diff_cutoff):
                return '-'

    def get_homologous_peptide(self, identity, peptide, homo_df, homo_faa_list):
        
        # Replace 'Var' with 'Ref' in the identity
        homo_id_list = []
        homo_seq_list = []
        
        identity = identity.replace('Var', 'Ref')
        homo_peptides = homo_df.loc[homo_df['Identity'] == identity, 'Peptide'].tolist()
        if len(homo_peptides) == 0:
            return 'N', 'N', 'N'

        similarity_scores = [SequenceMatcher(None, peptide, homo_peptide).ratio() for homo_peptide in homo_peptides]
        max_similarity_idx = similarity_scores.index(max(similarity_scores))
        homo_peptide = homo_peptides[max_similarity_idx]
    
        for pep_id, homo_seq, homo_id  in homo_faa_list:
            if identity == pep_id  and homo_peptide == pep_seq:
                homo_id_list.append(homo_id)
                homo_seq_list.append(homo_seq)
                
        print("homo_id_list, homo_seq_list",homo_id_list, homo_seq_list)
        if  len(homo_id_list) == 0:
            return 'N','N','N'
        else:
            return 'Y', homo_peptide, homo_id

    def get_affinity(self, wildtype_peptide, homo_peptide,  identity, mhc_genetype, ref_df, keyword):
        """
        keyword='Affinity(nM)' if MHCii, else 'Aff(nM)'
        """
        # Replace 'Var' with 'Ref' in the identity
        identity = identity.replace('Var', 'Ref')
        
        if wildtype_peptide != '-':
            affinity = ref_df.loc[(ref_df['Peptide'] == wildtype_peptide) & (ref_df['Identity'] == identity) & (ref_df['MHC'] == mhc_genetype), keyword].values
            return affinity[0]
        else:
            if homo_peptide != 'N':
                affinity = ref_df.loc[(ref_df['Peptide'] == homo_peptide) & (ref_df['Identity'] == identity) & (ref_df['MHC'] == mhc_genetype), keyword].values
                return affinity[0]
            else:
                return '-'
    
    def load_tmr_data(self):
        tmr_file = os.path.join(self.data_dir, '04.TransMembrane.DeepTMHMM', 'TMRs.gff3')
        with open(tmr_file, 'r') as f:
            self.tmr_data = f.readlines()
            
    def load_homo_faa(self, homo_netmhc_faa_file):
        self.homo_netmhc_faa_list = []
        for record in SeqIO.parse(homo_netmhc_faa_file, 'fasta'):
            pep_id = record.id
            pep_description = record.description.split("_")[-1]
            pep_seq = record.seq
            self.homo_netmhc_faa_list.append([pep_id, pep_seq, pep_description])
        
    def load_dataframes(self):
        if self.mhc_genotype in ['mhci', 'all']:
            self.netmhcpan_df = pd.read_csv(os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCpan.csv'))
            ref_netmhcpan_file = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCpan.csv')
            if os.path.exists(ref_netmhcpan_file):
                self.ref_netmhcpan_df = pd.read_csv(ref_netmhcpan_file)
            self.homo_netmhcpan_df = pd.read_csv(os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCpan.csv'))
            self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCpan_homologous.faa'))

        if self.mhc_genotype in ['mhcii', 'all']:
            self.netmhciipan_df = pd.read_csv(os.path.join(self.data_dir, '02.protein_antigen_prediction_var', 'parsed', f'{self.sample_name}_netMHCIIpan.csv'))
            ref_netmhciipan_file = os.path.join(self.data_dir, '02.protein_antigen_prediction_ref', 'parsed', f'{self.sample_name}_netMHCIIpan.csv')
            if os.path.exists(ref_netmhciipan_file):
                self.ref_netmhciipan_df = pd.read_csv(ref_netmhciipan_file)
            self.homo_netmhciipan_df = pd.read_csv(os.path.join(self.data_dir, '03.homologous', 'parsed', f'{self.sample_name}_homologous_netMHCIIpan.csv'))
            self.load_homo_faa(os.path.join(self.data_dir, '03.homologous', f'{self.sample_name}_netMHCIIpan_homologous.faa'))

    def summarize(self):
        self.load_tmr_data()
        self.load_dataframes()
        if self.mhc_genotype in ['mhci', 'all']:
            self.netmhcpan_df['TransMemb'] = self.netmhcpan_df.apply(lambda row: self.get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide'])), axis=1)
            self.netmhcpan_df['InCutmerRate'] = '-'
            self.netmhcpan_df['InCutmerRegion'] = '-'
            self.netmhcpan_df['Wildtype_peptide'] = self.netmhcpan_df.apply(lambda row: self.get_wildtype_peptide(row['Identity'], row['Peptide'], self.ref_netmhcpan_df), axis=1)
            self.netmhcpan_df['HomoExsit'] = '-'
            self.netmhcpan_df['Homo_id'] = '-'
            self.netmhcpan_df['Homo_peptide'] = '-'
            self.homo_netmhcpan_df['Identity'] = self.homo_netmhcpan_df['Identity'].str.split('_').str[0]
            results = self.netmhcpan_df.loc[self.netmhcpan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(
                lambda row: self.get_homologous_peptide(row['Identity'], row['Peptide'], self.homo_netmhcpan_df, self.homo_netmhc_faa_list), 
                axis=1, result_type='expand')

            for idx, (homo_exist, homo_peptide, homo_id) in results.iterrows():
                self.netmhcpan_df.loc[idx, 'HomoExsit'] = homo_exist
                self.netmhcpan_df.loc[idx, 'Homo_peptide'] = homo_peptide
                self.netmhcpan_df.loc[idx, 'Homo_id'] = homo_id
            self.netmhcpan_df['Aff(nM)_competitor'] = self.netmhcpan_df.apply(lambda row: self.get_affinity(row['Wildtype_peptide'], row['Homo_peptide'], row['Identity'], row["MHC"], self.ref_netmhcpan_df, 'Aff(nM)'), axis=1)

        if self.mhc_genotype in ['mhcii', 'all']:
            self.netmhciipan_df['TransMemb'] = self.netmhciipan_df.apply(lambda row: self.get_transmembrane_region(row['Identity'], row['Pos'], len(row['Peptide'])), axis=1)
            self.netmhciipan_df['InCutmerRate'] = '-'
            self.netmhciipan_df['InCutmerRegion'] = '-'
            self.netmhciipan_df['Wildtype_peptide'] = self.netmhciipan_df.apply(lambda row: self.get_wildtype_peptide(row['Identity'], row['Peptide'], self.ref_netmhciipan_df), axis=1)
            self.netmhciipan_df['HomoExsit'] = '-'
            self.netmhciipan_df['Homo_id'] = '-'
            self.netmhciipan_df['Homo_peptide'] = '-'
            self.homo_netmhciipan_df['Identity'] = self.homo_netmhciipan_df['Identity'].str.split('_').str[0]
            results = self.netmhciipan_df.loc[self.netmhciipan_df['Wildtype_peptide'] == '-', ['Peptide', 'Identity']].apply(
                lambda row: self.get_homologous_peptide(row['Identity'], row['Peptide'], self.homo_netmhciipan_df, self.homo_netmhc_faa_list), 
                axis=1, result_type='expand'
            )
            for idx, (homo_exist, homo_peptide, homo_id) in results.iterrows():
                self.netmhciipan_df.loc[idx, 'HomoExsit'] = homo_exist
                self.netmhciipan_df.loc[idx, 'Homo_peptide'] = homo_peptide
                self.netmhciipan_df.loc[idx, 'Homo_id'] = homo_id
            
            self.netmhciipan_df['Aff(nM)_competitor'] = self.netmhciipan_df.apply(lambda row: self.get_affinity(row['Wildtype_peptide'], row['Homo_peptide'], row['Identity'], row["MHC"], self.ref_netmhciipan_df, 'Affinity(nM)'), axis=1)


        # 删除 Peptide 与 Wildtype_peptide 完全相同的行
        self.netmhcpan_df = self.netmhcpan_df[self.netmhcpan_df['Peptide'] != self.netmhcpan_df['Wildtype_peptide']]
        self.netmhciipan_df = self.netmhciipan_df[self.netmhciipan_df['Peptide'] != self.netmhciipan_df['Wildtype_peptide']]

        # 删除不需要的列
        netmhcpan_columns_to_drop = ['Of', 'Gp', 'Gl', 'Ip', 'Il']
        netmhciipan_columns_to_drop = ['Of', 'Exp_Bind']
        self.netmhcpan_df.drop(columns=[col for col in netmhcpan_columns_to_drop if col in self.netmhcpan_df.columns], inplace=True)
        self.netmhciipan_df.drop(columns=[col for col in netmhciipan_columns_to_drop if col in self.netmhciipan_df.columns], inplace=True)


        # 计算 Aff(nM)_competitor/Aff(nM) 的比值
        self.netmhcpan_df['Aff(nM)_competitor/Aff(nM)'] = self.netmhcpan_df.apply(
            lambda row: round(float(row['Aff(nM)_competitor']) / float(row['Aff(nM)']), 2) if row['Aff(nM)_competitor'] != '-' and row['Aff(nM)'] != '-' else '-', axis=1
        )
        self.netmhciipan_df['Aff(nM)_competitor/Affinity(nM)'] = self.netmhciipan_df.apply(
            lambda row: round(float(row['Aff(nM)_competitor']) / float(row['Affinity(nM)']), 2) if row['Aff(nM)_competitor'] != '-' and row['Affinity(nM)'] != '-' else '-', axis=1
        )

        
        output_deliverable_dir = os.path.join(self.output_dir, 'Deliverable')
        if os.path.exists(output_deliverable_dir):
            shutil.rmtree(output_deliverable_dir)
        os.mkdir(output_deliverable_dir)

        if self.mhc_genotype in ['mhci', 'all']:
            self.netmhcpan_df.to_csv(os.path.join(output_deliverable_dir, f'{self.prefix}_netMHCpan_deliverable.csv'), index=False)
        if self.mhc_genotype in ['mhcii', 'all']:
            self.netmhciipan_df.to_csv(os.path.join(output_deliverable_dir, f'{self.prefix}_netMHCIIpan_deliverable.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Summarize netMHC results.')
    parser.add_argument('-d', '--data_dir', type=str, help='Path to the data directory')
    parser.add_argument('-o', '--output_dir', type=str, help='Path to the output directory')
    parser.add_argument('-p', '--prefix', type=str, help='Prefix of the output result')
    parser.add_argument('-n', '--sample_name', type=str, help='Sample name of this sample')
    parser.add_argument('-m', '--mhc_genotype', choices=['mhci', 'mhcii', 'all'], help='MHC genotypes')
    args = parser.parse_args()

    summarizer = NetMHCSummarizer(args.data_dir, args.output_dir, args.prefix, args.sample_name, args.mhc_genotype)
    summarizer.summarize()
