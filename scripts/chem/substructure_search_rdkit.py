#>  Copyright (c) 2020
#>  Heidelberg University and Heidelberg Institute of Theoretical Studies (HITS, www.h-its.org)
#>  Schloss-Wolfsbrunnenweg 35
#>  69118 Heidelberg, Germany
#>
#>  Please send your contact address to get information on updates and
#>  new features to "mcmsoft@h-its.org". Questions will be
#>  answered as soon as possible.
#>  References:
#>  A rapid identification of hit molecules for target proteins via physico-chemical descriptors.
#>  (2013) Phys. Chem. Chem. Phys., 15, 9107-9116.
#>  Authors: Goutam Mukherjee and B. Jayaram
#>  Version 1.0 (April 2013)

#>  RASPD+: Fast protein-ligand binding free energy prediction using simplified physicochemical features
#>  Authors: Stefan Holderbach, Lukas Adam, B. Jayaram, Rebecca C. Wade, Goutam Mukherjee
#>  Version 1.0 (June 2020)
#>  ChemRxiv preprint (https://doi.org/10.26434/chemrxiv.12636704.v1), 2020

#>  Author: Stefan Holderbach, Heidelberg University
#>  Requirements: python, rdkit 
#>  How to run:
#>  source activate /exports/scratch/anaconda3
#>  python substructure_search_rdkit.py -q query.txt -t target.txt -o output.txt
#>  conda deactivate

from __future__ import print_function
import argparse
from rdkit import Chem


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", required=True, type=str, help="File containing query string (SMILES)")
    parser.add_argument("-t", "--target", required=True, type=str, help="Path to file with smiles of ligands to be checked")
    parser.add_argument("-o", "--output", required=True, type=str, help="Path to output file")
    args = parser.parse_args()

    with open(args.target, "r") as target_file, open(args.output, 'w') as outfile, open(args.query, 'r') as query_file:
        query = query_file.read().strip()
        try:
            query_canon = Chem.MolFromSmiles(Chem.CanonSmiles(query))
        except:
            print("Invalid SMILES code for the query")
            exit(1)
        for ln, line in enumerate(target_file):
            smile = line.split("\t")[0]
            try:
                canon_smile = Chem.CanonSmiles(smile)
                mol = Chem.MolFromSmiles(canon_smile)
                out = mol.HasSubstructMatch(query_canon)
            except:
                print("invalid smile or failure during canonicalization in line {}".format(ln + 1))
                out = "error"
            finally:
                outfile.write("{}\n".format(out))

   
