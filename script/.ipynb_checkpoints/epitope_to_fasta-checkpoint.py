import os
import argparse
import re
from collections import defaultdict

def parse_mhc_species_mapping(filename):
    mapping = {}
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('|')
            if len(parts) < 2:
                continue
            species = parts[0]
            mhc_genotype = parts[1].strip()
            mapping[mhc_genotype] = species
    return mapping

# 表头模式1
MHCII_header = re.compile(r'Pos\s+MHC\s+Peptide\s+Of\s+Core\s+Core_Rel\s+Identity\s+Score_EL\s+%Rank_EL\s+Exp_Bind\s+Score_BA\s+Affinity\(nM\)\s+%Rank_BA\s+BindLevel')
# 表头模式2
MHCI_header= re.compile(r'Pos\s+MHC\s+Peptide\s+Core\s+Of\s+Gp\s+Gl\s+Ip\s+Il\s+Icore\s+Identity\s+Score_EL\s+%Rank_EL\s+Score_BA\s+%Rank_BA\s+Aff\(nM\)\s+BindLevel')

def determine_header_type(header_line):
    if MHCI_header.match(header_line):
        return "MHCI"
    elif MHCII_header.match(header_line):
        return "MHCII"
    else:
        return "Unknown"

def parse_netmhcpan_results(input_file, mhc_species_mapping):
    species_dict = defaultdict(list)
    
    with open(input_file, 'r') as infile:
        header_type = None
    
        for line in infile:
            print(line)
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-'):
                continue
            if header_type is None:
                header_type = determine_header_type(line)
                continue

            if header_type == "MHCI":
                parts = line.split()
                if  len(parts)  == 16 or len(parts)  == 17:
                    if line.startswith(" Pos"):
                        continue
                    allele = parts[1]
                    peptide = parts[2]
                    pos_start = int(parts[0])
                    pos_end = pos_start + len(peptide)
                    Identity = parts[10]
                
            elif header_type == "MHCII":
                parts = line.split()
                if len(parts) == 13 or len(parts) == 14:
                    if line.startswith(" Pos"):
                        continue
                    allele = parts[1]
                    peptide = parts[2]
                    pos_start = int(parts[0])
                    pos_end = pos_start + len(peptide)
                    Identity = parts[6]
            else:
                print("Unknown header type")
                return species_dict

            print("allele", allele, 'peptide', peptide)
            
            species = determine_species(allele, mhc_species_mapping)
            if species:
                species_dict[species].append((Identity, pos_start, pos_end, peptide))
            else:
                print(f"No match found for allele: {allele}")    
    return species_dict

def determine_species(allele, mhc_species_mapping):
    patterns = [
        allele,                       # Direct match
        allele.replace('*', ''),      # Remove asterisk
        re.sub(r'\d{2}$', '', allele) # Remove trailing two digits
    ]

    for pattern in patterns:
        for mhc_genotype_pattern, species in mhc_species_mapping.items():
            regex_pattern = re.escape(mhc_genotype_pattern).replace(r'\*', '*')
            if re.match(regex_pattern, pattern):
                return species
    return None

def write_to_fasta(species_dict, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for species, entries in species_dict.items():
        output_file = os.path.join(output_dir, f"{species}.faa")
        with open(output_file, 'w') as outfile:
            for idx, (allele, peptide) in enumerate(entries):
                outfile.write(f">{Identity}_{pos_start}_{pos_end}\n")
                outfile.write(f"{peptide}\n")

def main():
    parser = argparse.ArgumentParser(description='Process NetMHCpan results and categorize by MHC species.')
    parser.add_argument('-i', '--input_file', metavar='INPUT_FILE', type=str,
                        help='input file containing NetMHCpan results')
    parser.add_argument('--mhci_genetype_list', metavar='MHCI_GENETYPE_LIST', type=str,
                        help='file containing MHCI genotypes and species mappings')
    parser.add_argument('--mhcii_genetype_list', metavar='MHCII_GENETYPE_LIST', type=str,
                        help='file containing MHCII genotypes and species mappings')
    parser.add_argument('-o', '--output_dir', metavar='OUTPUT_DIR', type=str, default='output_fasta_files',
                        help='output directory for generated FASTA files (default: output_fasta_files)')
    
    args = parser.parse_args()
    
    mhci_species_mapping = parse_mhc_species_mapping(args.mhci_genetype_list)
    mhcii_species_mapping = parse_mhc_species_mapping(args.mhcii_genetype_list)
    
    mhc_species_mapping = {**mhci_species_mapping, **mhcii_species_mapping}
    
    species_dict = parse_netmhcpan_results(args.input_file, mhc_species_mapping)
    
    write_to_fasta(species_dict, args.output_dir)
    
    print(f"FASTA files have been generated in {args.output_dir}.")

if __name__ == "__main__":
    main()
