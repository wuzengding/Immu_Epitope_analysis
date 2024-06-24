import pandas as pd
import argparse

def process_mhc_file(file_path, sample_id):
    # Load the Excel file
    df = pd.read_table(file_path)

    # Filter the row based on the sample ID
    row = df[df['样本编号'].str.contains(sample_id)]

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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process MHC genotyping data from an Excel file.')
    parser.add_argument('-f','--file_path', type=str, help='Path to the MHC genetype file')
    parser.add_argument('-s','--sample_id', type=str, help='Sample ID to filter')

    args = parser.parse_args()

    result = process_mhc_file(args.file_path, args.sample_id)
    print(result)
