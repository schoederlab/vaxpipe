import pandas as pd
from argparse import ArgumentParser

parser = ArgumentParser(prog='Design the mutations using Rosetta')
parser.add_argument('-i', '--input', help='mutation .csv file', type=str)
parser.add_argument('-o', '--output', help='mutation .txt file', type=str)
args = parser.parse_args()

aa_mapping = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}
chain = 'A'

mutation_csv = pd.read_csv(args.input, sep=',')

mutations = []

for index, row in mutation_csv.iterrows():
    position = row['Position']
    mutation_aa = row['Mutation']
    mutation_aa_3letter = aa_mapping[mutation_aa]

    mutation_string = f"{position}{chain}_{mutation_aa_3letter}"
    mutations.append(f"{mutation_string}")

# Method 2: Write to a file that bash can source
mutations_file = args.output
with open(mutations_file, 'w') as f:
    f.write(f'mutations=({" ".join([f'"{mut}"' for mut in mutations])})\n')

