import pandas as pd
from argparse import ArgumentParser
import os
import re

parser = ArgumentParser(prog='Design the mutations using Rosetta')
parser.add_argument('-i', '--input', help='mutation .csv file', type=str, required=True)
parser.add_argument('-o', '--output', help='output filename (.txt, row number inferred)', type=str, required=True)
args = parser.parse_args()

aa_mapping = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}
chain = 'A'

# Read mutation CSV
mutation_csv = pd.read_csv(args.input, sep=',')

# Ensure output directory exists
output_dir = os.path.dirname(args.output)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)

# Extract row number from output filename (e.g. "5.txt" → row_num = 5)
basename = os.path.basename(args.output)
match = re.match(r"(\d+)\.txt$", basename)
if not match:
    raise ValueError(f"Output filename must be like '5.txt', got '{basename}'")

row_num = int(match.group(1))

# Select the correct row (1-based index)
if row_num < 1 or row_num - 1 > len(mutation_csv):
    raise IndexError(f"Row {row_num} is out of range (CSV has {len(mutation_csv)} rows)")

row = mutation_csv.iloc[row_num - 2]
position = row["Position"]
mutation_aa = row["Mutation"]
mutation_aa_3letter = aa_mapping[mutation_aa]

mutation_string = f"{position}{chain}_{mutation_aa_3letter}"

# Write to specified output file
with open(args.output, 'w') as f:
    f.write(mutation_string + "\n")

