import argparse
import pandas as pd
from Bio import SeqIO
import os

class BlastpParser:
    def __init__(self, blastp_file,  pident_cutoff=0.8):
        self.blastp_file = blastp_file
        #self.fasta_file = fasta_file
        #self.output_file = output_file
        self.pident_cutoff = pident_cutoff

    def read_blastp_results(self):
        self.blastp_df = pd.read_csv(self.blastp_file, sep='\s+', header=None,
                                     names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
                                            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

    def filter_results(self):
        # 先过滤pident大于cutoff的行
        pident_filtered_df = self.blastp_df[self.blastp_df['pident'] >= self.pident_cutoff]

        max_length_df = pident_filtered_df.loc[pident_filtered_df.groupby('qseqid')['length'].idxmax()]
        max_pident_df = max_length_df.loc[max_length_df.groupby('qseqid')['pident'].idxmax()]
        final_df = max_pident_df.loc[max_pident_df.groupby('qseqid')['evalue'].idxmin()]
        return final_df

    def read_fasta_sequences(self, fasta_file):
        protein_sequences = {record.id.split("|")[1]:str(record.seq) for record in SeqIO.parse(fasta_file, 'fasta')}
        return protein_sequences
        
    def extract_peptides(self, df, sequences, output_dir, prefix ):
        outfile = os.path.join(output_dir, prefix+"_homologous.faa")
        with open(outfile, 'w') as out_file:
            for idx, row in df.iterrows():
                subject_seq = sequences[row['sseqid']]
                start = row['sstart'] - 1
                end = row['send']
                peptide_seq = subject_seq[start:end]
                seq_id_parts =row["qseqid"].split("_")
                seq_id = seq_id_parts[0]+ " " + "_".join([seq_id_parts[1], seq_id_parts[2], row["sseqid"]])
                out_file.write(f'>{seq_id}\n{peptide_seq}\n')

    def run(self, fasta_file, output_dir, prefix):
        self.read_blastp_results()
        final_df = self.filter_results()
        protein_sequences = self.read_fasta_sequences(fasta_file)
        self.extract_peptides(final_df, protein_sequences, output_dir, prefix)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract peptides from BLASTP results and FASTA file.')
    parser.add_argument('-b','--blastp_file', help='Path to BLASTP results file')
    parser.add_argument('-f','--protein_fasta_file', help='Path to protein db FASTA file')
    parser.add_argument('-o','--output_dir', help='Path to output FASTA file')
    parser.add_argument('-p','--prefix', help="Prefix of output FASTA file")
    parser.add_argument('-c','--pident_cutoff', type=float, default=0.8, help='Minimum pident cutoff value (default: 0.8)')
    args = parser.parse_args()

    parser = BlastpParser(args.blastp_file,  args.pident_cutoff)
    parser.run(args.protein_fasta_file, args.output_dir, args.prefix)