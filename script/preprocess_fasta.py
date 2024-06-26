import argparse
from datetime import datetime
import os

class FastaProcessor:
    def __init__(self, fasta_file, outpath, prefix):
        self.fasta_file = fasta_file
        self.outpath = outpath
        self.prefix = prefix
        self.today = datetime.today().strftime('%y%m%d')

    def read_fasta(self):
        """
        读取FASTA文件并返回Ref和Var序列列表
        """
        with open(self.fasta_file, 'r') as file:
            lines = file.readlines()
        self.ref_sequences = []
        self.var_sequences = []
        for line in lines:
            if line.startswith('>Ref'):
                self.ref_sequences.append((line.strip(), ''))
            elif line.startswith('>Var'):
                self.var_sequences.append((line.strip(), ''))
            elif line.startswith(">"):
                line = ">Ref" + line[1:] 
                self.var_sequences.append((line.strip(), ''))
            else:
                if self.ref_sequences and self.ref_sequences[-1][1] == '':
                    self.ref_sequences[-1] = (self.ref_sequences[-1][0], line.strip())
                elif self.var_sequences and self.var_sequences[-1][1] == '':
                    self.var_sequences[-1] = (self.var_sequences[-1][0], line.strip())

    def _read_fasta(self):
        """
        读取FASTA文件并返回Ref和Var序列列表
        """
        with open(self.fasta_file, 'r') as file:
            lines = file.readlines()
        self.ref_sequences = []
        self.var_sequences = []
        multi_ref = 0
        for line in lines:
            if line.startswith('>Ref'):
                if multi_ref == 0:
                    self.ref_sequences.append([(line.strip(), '')])
                    multi_ref += 1
                elif multi_ref >= 1:
                    self.ref_sequences[-1].append((line.strip(), ''))
            elif line.startswith('>Var'):
                self.var_sequences.append([(line.strip(), '')])
                multi_ref = 0
            elif line.startswith(">"):
                line = ">Ref" + line[1:] 
                self.var_sequences.append([(line.strip(), '')])
            else:
                if self.ref_sequences and self.ref_sequences[-1][-1][1] == '':
                    self.ref_sequences[-1][-1] = (self.ref_sequences[-1][-1][0], line.strip())
                elif self.var_sequences and self.var_sequences[-1][-1][1] == '':
                    self.var_sequences[-1][-1] = (self.var_sequences[-1][-1][0], line.strip())

                    
    def transform_ids(self, sequences, type_):
        """
        转换序列ID并生成新的ID对应关系
        """
        transformed_sequences = []
        id_transformations = []
        for count, (old_id, seq) in enumerate(sequences, start=1):
            new_id = f">{type_}{self.today}{count:04d}"
            transformed_sequences.append((new_id, seq))
            id_transformations.append(f"{old_id} -> {new_id}")
        return transformed_sequences, id_transformations

    def _transform_ids(self, sequences, type_):
        """
        转换序列ID并生成新的ID对应关系
        """
        transformed_sequences = []
        id_transformations = []
        for count, old_id_seq_list in enumerate(sequences, start=1):
            for (old_id, seq) in  old_id_seq_list:
                new_id = f">{type_}{self.today}{count:04d}"
                transformed_sequences.append((new_id, seq))
                id_transformations.append(f"{old_id} -> {new_id}")
        return transformed_sequences, id_transformations
        
    def write_fasta(self, sequences, filename):
        """
        将序列写入FASTA文件
        """
        outfile = os.path.join(self.outpath, filename)
        with open(outfile, 'w') as file:
            for seq in sequences:
                file.write(f"{seq[0]}\n{seq[1]}\n")

    def write_transformations(self, transformations):
        """
        将ID转换关系写入文件
        """
        outfile = os.path.join(self.outpath, f"{self.prefix}_id_transformation.txt")
        with open(outfile, 'w') as file:
            for transformation in transformations:
                file.write(f"{transformation}\n")

    def process(self):
        """
        主函数：处理FASTA文件并输出结果
        """
        # 读取FASTA文件
        self._read_fasta()
        # 转换Ref序列的ID
        transformed_ref_sequences, ref_id_transformations = self._transform_ids(self.ref_sequences, 'Ref')
        # 转换Var序列的ID
        transformed_var_sequences, var_id_transformations = self._transform_ids(self.var_sequences, 'Var')
       
        '''
        transformed_sequences = transformed_ref_sequences + transformed_var_sequences
        self.write_fasta(transformed_sequences, f"{self.prefix}_ref_var_seq.faa")
        '''
        # 写入Ref序列到输出文件
        self.write_fasta(transformed_ref_sequences, f"{self.prefix}_ref_seq.faa")
        # 写入Var序列到输出文件
        self.write_fasta(transformed_var_sequences, f"{self.prefix}_var_seq.faa")
        # 写入ID转换关系到输出文件
        self.write_transformations(ref_id_transformations + var_id_transformations)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process FASTA files')
    parser.add_argument('-f', '--fasta-file', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--outpath', required=True, help='Output directory')
    parser.add_argument('-p', '--prefix', required=True, help='Output file prefix')
    args = parser.parse_args()
    processor = FastaProcessor(args.fasta_file, args.outpath, args.prefix)
    processor.process()