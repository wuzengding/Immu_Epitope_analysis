import pandas as pd
import argparse
import os

def process_mhc_file_sakit(file_path, sample_id):
    # Load the Excel file
    df = pd.read_excel(file_path)

    # Filter the row based on the sample ID
    row = df[df['样本姓名'].str.contains(sample_id)]

    if row.empty:
        return f"Sample ID {sample_id} not found."

    # Extract relevant columns
    a1, a2, b1, b2, c1, c2 = [mhc.replace('*','') for mhc in row[['A1*', 'A2*', 'B1*', 'B2*', 'C1*', 'C2*']].values[0]]
    dqa1_1, dqa1_2, dqb1_1, dqb1_2 = [mhc.replace('*','').replace(':','') for mhc in row[['DQA1-1*', 'DQA1-2*', 'DQB1-1*', 'DQB1-2*']].values[0]]
    dpa1_1, dpa1_2, dpb1_1, dpb1_2 = [mhc.replace('*','').replace(':','') for mhc in row[['DPA1-1*', 'DPA1-2*', 'DPB1-1*', 'DPB1-2*']].values[0]]
    drb1_1, drb1_2, drb345_1, drb345_2 = [mhc.replace('*','_').replace(':','') for mhc in  row[['DRB1-1*', 'DRB1-2*', 'DRB345*-1', 'DRB345*-2']].values[0]]

    # Format the strings
    string_a = [f'HLA-{a1[:-3]}', f'HLA-{a2[:-3]}', f'HLA-{b1[:-3]}',f'HLA-{b2[:-3]}', f'HLA-{c1[:-3]}', f'HLA-{c2[:-3]}']
    string_b = [f'HLA-{dqa1_1[:-2]}-{dqb1_1[:-2]}', f'HLA-{dqa1_1[:-2]}-{dqb1_2[:-2]}', f'HLA-{dqa1_2[:-2]}-{dqb1_1[:-2]}', f'HLA-{dqa1_2[:-2]}-{dqb1_2[:-2]}']
    string_c = [f'HLA-{dpa1_1[:-2]}-{dpb1_1[:-2]}', f'HLA-{dpa1_1[:-2]}-{dpb1_2[:-2]}', f'HLA-{dpa1_2[:-2]}-{dpb1_1[:-2]}', f'HLA-{dpa1_2[:-2]}-{dpb1_2[:-2]}']
    string_d = [f'{drb1_1[:-2]}', f'{drb1_2[:-2]}', f'{drb345_1[:-2]}',f'{drb345_2[:-2]}']

    # Combine the strings
    string_mhc = string_a + string_b + string_c + string_d
    string_mhc_final = []
    for mhc_type in string_mhc:
        if mhc_type not in string_mhc_final:
            string_mhc_final.append(mhc_type)

    return ','.join(string_mhc_final)

def clean_allele(allele_text, gene_prefix):
    """Clean allele text by removing gene prefix and limiting to 2-field notation"""
    # Remove gene prefix if present
    if gene_prefix in allele_text:
        allele_text = allele_text.replace(gene_prefix, '')
    
    # Remove * if present
    allele_text = allele_text.replace('*', '')
    
    # Split by colon and take only first two fields
    parts = allele_text.split(':')
    if len(parts) >= 2:
        return f"{parts[0]}:{parts[1]}"
    return allele_text

def process_mhc_file_scanno2(file_path, sample_id, output_dir='.'):
    # Load the Excel file
    df = pd.read_excel(file_path)

    # Filter the row based on the sample ID
    row = df[df['样本姓名'].str.contains(sample_id)]

    if row.empty:
        return f"Sample ID {sample_id} not found."

    # Extract relevant columns
    mhc_i_data = []
    mhc_ii_data = []
    
    # Class I (A, B, C)
    a1_value = row['A1*'].values[0]
    a2_value = row['A2*'].values[0]
    
    if pd.notna(a1_value) and a1_value.strip():
        mhc_i_data.append(f"custom\tHLA-A*{clean_allele(a1_value, 'A')}")
    
    if pd.notna(a2_value) and a2_value.strip():
        mhc_i_data.append(f"custom\tHLA-A*{clean_allele(a2_value, 'A')}")
    
    b1_value = row['B1*'].values[0]
    b2_value = row['B2*'].values[0]
    
    if pd.notna(b1_value) and b1_value.strip():
        mhc_i_data.append(f"custom\tHLA-B*{clean_allele(b1_value, 'B')}")
    
    if pd.notna(b2_value) and b2_value.strip():
        mhc_i_data.append(f"custom\tHLA-B*{clean_allele(b2_value, 'B')}")
    
    c1_value = row['C1*'].values[0]
    c2_value = row['C2*'].values[0]
    
    if pd.notna(c1_value) and c1_value.strip():
        mhc_i_data.append(f"custom\tHLA-C*{clean_allele(c1_value, 'C')}")
    
    if pd.notna(c2_value) and c2_value.strip():
        mhc_i_data.append(f"custom\tHLA-C*{clean_allele(c2_value, 'C')}")
    
    # Class II
    dqa1_1_value = row['DQA1-1*'].values[0]
    dqa1_2_value = row['DQA1-2*'].values[0]
    
    if pd.notna(dqa1_1_value) and dqa1_1_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DQA1*{clean_allele(dqa1_1_value, 'DQA1')}")
    
    if pd.notna(dqa1_2_value) and dqa1_2_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DQA1*{clean_allele(dqa1_2_value, 'DQA1')}")
    
    dqb1_1_value = row['DQB1-1*'].values[0]
    dqb1_2_value = row['DQB1-2*'].values[0]
    
    if pd.notna(dqb1_1_value) and dqb1_1_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DQB1*{clean_allele(dqb1_1_value, 'DQB1')}")
    
    if pd.notna(dqb1_2_value) and dqb1_2_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DQB1*{clean_allele(dqb1_2_value, 'DQB1')}")
    
    dpa1_1_value = row['DPA1-1*'].values[0]
    dpa1_2_value = row['DPA1-2*'].values[0]
    
    if pd.notna(dpa1_1_value) and dpa1_1_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DPA1*{clean_allele(dpa1_1_value, 'DPA1')}")
    
    if pd.notna(dpa1_2_value) and dpa1_2_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DPA1*{clean_allele(dpa1_2_value, 'DPA1')}")
    
    dpb1_1_value = row['DPB1-1*'].values[0]
    dpb1_2_value = row['DPB1-2*'].values[0]
    
    if pd.notna(dpb1_1_value) and dpb1_1_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DPB1*{clean_allele(dpb1_1_value, 'DPB1')}")
    
    if pd.notna(dpb1_2_value) and dpb1_2_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DPB1*{clean_allele(dpb1_2_value, 'DPB1')}")
    
    drb1_1_value = row['DRB1-1*'].values[0]
    drb1_2_value = row['DRB1-2*'].values[0]
    
    if pd.notna(drb1_1_value) and drb1_1_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DRB1*{clean_allele(drb1_1_value, 'DRB1')}")
    
    if pd.notna(drb1_2_value) and drb1_2_value.strip():
        mhc_ii_data.append(f"custom\tHLA-DRB1*{clean_allele(drb1_2_value, 'DRB1')}")
    
    # DRB345 handling
    drb345_1 = row['DRB345*-1'].values[0]
    drb345_2 = row['DRB345*-2'].values[0]
    
    if pd.notna(drb345_1) and drb345_1.strip():
        if 'DRB3' in drb345_1:
            mhc_ii_data.append(f"custom\tHLA-DRB3*{clean_allele(drb345_1.split('*')[1], '')}")
        elif 'DRB4' in drb345_1:
            mhc_ii_data.append(f"custom\tHLA-DRB4*{clean_allele(drb345_1.split('*')[1], '')}")
        elif 'DRB5' in drb345_1:
            mhc_ii_data.append(f"custom\tHLA-DRB5*{clean_allele(drb345_1.split('*')[1], '')}")
    
    if pd.notna(drb345_2) and drb345_2.strip():
        if 'DRB3' in drb345_2:
            mhc_ii_data.append(f"custom\tHLA-DRB3*{clean_allele(drb345_2.split('*')[1], '')}")
        elif 'DRB4' in drb345_2:
            mhc_ii_data.append(f"custom\tHLA-DRB4*{clean_allele(drb345_2.split('*')[1], '')}")
        elif 'DRB5' in drb345_2:
            mhc_ii_data.append(f"custom\tHLA-DRB5*{clean_allele(drb345_2.split('*')[1], '')}")
    
    # Remove duplicates while maintaining order
    mhc_i_data = list(dict.fromkeys(mhc_i_data))
    mhc_ii_data = list(dict.fromkeys(mhc_ii_data))
    
    # Write to files
    mhc_i_file = os.path.join(output_dir, f"{sample_id}_MHC-I.txt")
    mhc_ii_file = os.path.join(output_dir, f"{sample_id}_MHC-II.txt")
    
    with open(mhc_i_file, 'w') as f:
        f.write('\n'.join(mhc_i_data))
    
    with open(mhc_ii_file, 'w') as f:
        f.write('\n'.join(mhc_ii_data))
    
    return f"Files created: {mhc_i_file} and {mhc_ii_file}"

def show_file_content(file_path):
    """Function to display the content of the Excel file."""
    df = pd.read_excel(file_path)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process MHC genotyping data from an Excel file.')
    parser.add_argument('-f','--file_path', type=str, help='Path to the MHC genetype file', required=True)
    parser.add_argument('-s','--sample_id', type=str, help='Sample ID to filter', required=False)
    parser.add_argument('-w', '--show', action='store_true', help='Show the content of the file')
    parser.add_argument('-m', '--mode', type=str, choices=['SAKit', 'Scanno2'], 
                        default='SAKit', help='Output mode (SAKit or Scanno2)')
    parser.add_argument('-o', '--output_dir', type=str, default='.', 
                        help='Output directory for Scanno2 mode files')

    args = parser.parse_args()

    if args.show:
        # Show the content of the file if the --show flag is used
        file_content = show_file_content(args.file_path)
        print(file_content)
    elif args.sample_id:
        # Process the file based on the selected mode
        if args.mode == 'SAKit':
            result = process_mhc_file_sakit(args.file_path, args.sample_id)
            print(result)
        elif args.mode == 'Scanno2':
            result = process_mhc_file_scanno2(args.file_path, args.sample_id, args.output_dir)
            print(result)
    else:
        print("Please provide either --sample_id to process or --show to display the file.")

