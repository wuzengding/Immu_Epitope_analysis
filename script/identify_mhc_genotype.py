import argparse
import json

def load_mhc_file(file_path):
    """Load MHC gene types from a file and return as a dictionary."""
    mhc_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            species, gene = line.strip().split('|')
            mhc_dict[gene] = species
    return mhc_dict

def identify_mhc_genotype(genotype, mhci_dict, mhcii_dict):
    """Identify the class and species of a given MHC genotype."""
    
    if genotype in mhci_dict:
        return "MHC I", genotype, mhci_dict[genotype]
    elif genotype in mhcii_dict:
        return "MHC II", genotype, mhcii_dict[genotype]
    else:
        return "unrecognized", genotype, "unrecognized"

def main(hla_genotypes, mhci_file, mhcii_file):
    """Main function to identify MHC genotypes."""
    # Load MHC gene types from files
    mhci_dict = load_mhc_file(mhci_file)
    mhcii_dict = load_mhc_file(mhcii_file)

    # Initialize MHC I and MHC II genotype dictionaries
    mhci_genotypes = {}
    mhcii_genotypes = {}
    species_genotypes = {}

    # Split and process each HLA genotype
    genotypes_array = hla_genotypes.split(',')

    for genotype in genotypes_array:
        mhc_class, gene, species = identify_mhc_genotype(genotype, mhci_dict, mhcii_dict)
        if mhc_class == "MHC I":
            if mhc_class not in mhci_genotypes:
                mhci_genotypes[mhc_class] = []
            mhci_genotypes[mhc_class].append(gene)
        elif mhc_class == "MHC II":
            if mhc_class not in mhcii_genotypes:
                mhcii_genotypes[mhc_class] = []
            mhcii_genotypes[mhc_class].append(gene)
        
        if species != "unrecognized":
            if species not in species_genotypes:
                species_genotypes[species] = []
            species_genotypes[species].append(gene)

    # Construct result dictionary
    result = {
        "mhc_genotypes": {
            "MHC I": mhci_genotypes.get("MHC I", []),
            "MHC II": mhcii_genotypes.get("MHC II", [])
        },
        "species_genotypes": species_genotypes
    }
    
    # Print results as JSON
    print(json.dumps(result, indent=4))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Identify MHC genotypes and their classes.")
    parser.add_argument('-g', '--genotypes', required=True, help="Comma-separated list of HLA genotypes.")
    parser.add_argument('--mhci-file', default="lib/MHCI_genetype_list", help="Path to the MHC I genotype list file.")
    parser.add_argument('--mhcii-file', default="lib/MHCII_genetype_list", help="Path to the MHC II genotype list file.")
    args = parser.parse_args()

    main(args.genotypes, args.mhci_file, args.mhcii_file)
