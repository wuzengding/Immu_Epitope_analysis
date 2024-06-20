import re
import pandas as pd
import argparse
import os

class NetMHCParser:
    def __init__(self, tool_type, input_file, BindLevel="SBWB"):
        self.tool_type = tool_type
        self.input_file = input_file
        self.BindLevel = BindLevel
        
        if tool_type == "netMHCpan":
            self.header_keywords = ["Pos", "MHC", "Peptide", "Core", "Of", "Gp",\
                                    "Gl", "Ip", "Il", "Icore", "Identity", "Score_EL",\
                                    "%Rank_EL", "Score_BA", "%Rank_BA", "Aff(nM)", "BindLevel"]
        elif tool_type == "netMHCIIpan":
            self.header_keywords = ["Pos", "MHC", "Peptide", "Of", "Core", "Core_Rel", \
                                    "Identity", "Score_EL", "%Rank_EL", "Exp_Bind",\
                                    "Score_BA", "Affinity(nM)", "%Rank_BA", "BindLevel"]
        else:
            raise ValueError("Unknown tool type")

    def escape_keywords(self, keywords):
        keywords = [re.escape(keyword) for keyword in keywords]
        return keywords

    def parse_netmhc_output(self):
        data = []
        parsing = False
        header_detected = False
        end_sign = 0
        header_keywords = self.escape_keywords(self.header_keywords)
        start_pattern = re.compile(r"^\s*" + r"\s+".join(header_keywords) + r"\s*$")
        end_pattern = re.compile(r"^-+(-+)+$")  ## -+: 这个部分匹配行开头(^)的一个或多个(+)破折号(-)。
        
        with open(self.input_file,"r") as netmhc_f:
            for line in netmhc_f.readlines():
                if start_pattern.match(line):
                    parsing = True
                    header_detected = True
                    end_sign = 0
                    continue
                
                if end_pattern.match(line):
                    if end_sign == 0:
                        end_sign += 1
                    elif end_sign == 1:
                        end_sign = 0
                        parsing = False
                    continue
                
                if parsing:
                    line = re.sub(r'<=', '', line)
                    line = re.sub(r'<=', '', line)
                    line_item = line.split()
                    line_item_num = len(line_item)
                    header_item_num = len(self.header_keywords)
                    
                    if line_item_num < header_item_num:
                        line_item += [""] * (header_item_num - line_item_num)
                        data.append(line_item)
                    elif line_item_num == header_item_num:
                        data.append(line_item)
                    else:
                        print("error line:", line)
        return data

    def to_dataframe(self):
        data = self.parse_netmhc_output()
        df = pd.DataFrame(data, columns=self.header_keywords)
        return df
        
    def filter(self):
        df = self.to_dataframe()
        if self.BindLevel == "all":
            # 不进行任何过滤，返回原始的 DataFrame
            df_filtered =  df
        elif self.BindLevel == "SBWB" or self.BindLevel == "WBSB":
            df_filtered = df[(df["BindLevel"] == "SB") | (df["BindLevel"] == "WB")]
        elif self.BindLevel == "SB":
            df_filtered = df[df["BindLevel"] == "SB"]
        elif self.BindLevel == "WB":
            df_filtered = df[df["BindLevel"] == "WB"]
        else:
            # 如果输入的 BindLevel 不在预期的选项中，可以选择抛出异常或者返回原始的 DataFrame
            df_filtered = df
            raise ValueError("Invalid BindLevel value. Expected 'SB,WB', 'SB', 'WB', or 'all'.")
        #print(df_filtered)
        return df_filtered
        
    def ensure_path_exists(self, path):
        if not os.path.exists(path):
            os.makedirs(path)
            #print(f"Path '{path}' created.")
        else:
            #print(f"Path '{path}' already exists.")
            pass
            
    def to_csv(self, output_dir, prefix):
        df = self.to_dataframe()
        df = self.filter()
        self.ensure_path_exists(output_dir)
        df.to_csv(os.path.join(output_dir, prefix+"_"+self.tool_type+".csv"), index=False)
        
    def to_fasta(self, output_dir, prefix):
        df = self.to_dataframe()
        df = self.filter()
        fasta_lines = []
        for _, row in df.iterrows():
            identity = row['Identity']
            pos = int(row['Pos'])
            peptide = row['Peptide']
            length = len(peptide)
            fasta_header = f">{identity}_{pos}_{pos + length}"
            fasta_lines.append(fasta_header)
            fasta_lines.append(peptide)
        
        fasta_content = "\n".join(fasta_lines)
        self.ensure_path_exists(output_dir)
        with open(os.path.join(output_dir, prefix+"_"+ self.tool_type+ ".faa"), 'w') as fasta_file:
            fasta_file.write(fasta_content)
        

# Example usage
file_content_netmhcpan = """
---------------------------------------------------------------------------------------------------------------------------
 Pos         MHC        Peptide      Core Of Gp Gl Ip Il        Icore        Identity  Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
---------------------------------------------------------------------------------------------------------------------------
   7 HLA-A*11:01      SVESKMSNK SVESKMSNK  0  0  0  0  0    SVESKMSNK PB_1036_2_FASTK 0.7558220    0.127 0.617933    0.371    62.42 <=SB
   4 HLA-A*11:01       GSVSVESK GSV-SVESK  0  0  0  3  1     GSVSVESK PB_1036_2_FASTK 0.0213460    3.830 0.196121    5.920  5989.68
---------------------------------------------------------------------------------------------------------------------------
"""

file_content_netmhciipan = """
--------------------------------------------------------------------------------------------------------------------------------------------
 Pos           MHC              Peptide   Of        Core  Core_Rel        Identity      Score_EL %Rank_EL Exp_Bind      Score_BA  Affinity(nM) %Rank_BA  BindLevel
--------------------------------------------------------------------------------------------------------------------------------------------
   3     DRB1_0301      FGSVSVESKMSNKAG    3   VSVESKMSN     1.000 PB_1036_2_FASTK      0.182126     5.96       NA      0.347062       1169.87    11.16   <=WB     
  12     DRB1_0301       MSNKAGSFFWNLRQ    5   GSFFWNLRQ     0.427 PB_1036_2_FASTK      0.000110    95.00       NA      0.111220      15008.96    82.19       
--------------------------------------------------------------------------------------------------------------------------------------------

"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse netMHCpan/netMHCIIpan results.')
    parser.add_argument('-t','--tool_type', choices=['netMHCpan', 'netMHCIIpan'], required=True, help='Type of the netMHC tool')
    parser.add_argument('-i','--input_file', required=True, help='Path to the input file containing netMHC results')
    parser.add_argument('-o','--output_dir', required=True, help='Directory to save the output CSV file')
    parser.add_argument('-y','--output_type', choices=["csvfasta","csv","fasta"],default="csv",help="output file type")
    parser.add_argument('-p','--prefix', required=True, help='Prefix for the output CSV file name')
    parser.add_argument('-b','--bindlevel', choices=['SBWB', 'SB', 'WB', 'all'], default='SBWB', help="Filter with bindlevel")

    args = parser.parse_args()

    df = NetMHCParser(args.tool_type, args.input_file, args.bindlevel)
    if  args.output_type == "csv":
        df.to_csv(args.output_dir, args.prefix)
    elif args.output_type == "fasta":
        df.to_fasta(args.output_dir, args.prefix)
    elif args.output_type == "csvfasta":
        df.to_csv(args.output_dir, args.prefix)
        df.to_fasta(args.output_dir, args.prefix)

