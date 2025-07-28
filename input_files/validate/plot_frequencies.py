from Bio import SeqIO
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import pandas as pd

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 2

def read_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def read_wild_type(file_path):
    for record in SeqIO.parse(file_path, "fasta"):
        return str(record.seq)

def count_mutations(sequences, wild_type):
    position_counts = {}
    for seq in sequences:
        for i, (wt_residue, mut_residue) in enumerate(zip(wild_type, seq)):
            if wt_residue != mut_residue:
                if i not in position_counts:
                    position_counts[i] = {}
                if mut_residue in position_counts[i]:
                    position_counts[i][mut_residue] += 1
                else:
                    position_counts[i][mut_residue] = 1
    return position_counts

parser = ArgumentParser(prog='Plot Frequencies for evaluation of proposed mutations')
parser.add_argument('-i', '--input', help='fastafile', type=str)
parser.add_argument('-r', '--ref', type=str, help='path to the wildtype fasta file')
parser.add_argument('-o', '--output', type=str, default='output.png', help='filename for the resulting output. default: output.*')
args = parser.parse_args()

fasta_file = args.input  # Path to your FASTA file
wild_type_file = args.ref  # Path to your wild-type FASTA file

sequences = read_fasta(fasta_file)
wild_type = read_wild_type(wild_type_file)
mutation_counts = count_mutations(sequences, wild_type)


mutation_list = []
for position, mutations in mutation_counts.items():
    for mutation, count in mutations.items():
        mutation_list.append((position + 1, wild_type[position], mutation, count))


mutation_list.sort(key=lambda x: x[3], reverse=True)

# Create DataFrame from top 20 mutations
top_mutations = mutation_list[:20]
df = pd.DataFrame(top_mutations, columns=["Position", "WildType", "Mutation", "Count"])

# Save to CSV
csv_output_file = args.output.replace('.png', '.csv')
df.to_csv(csv_output_file, index=False)

# Create x and y values
labels = [f"{orig}{pos}{mut}" for pos, orig, mut, score in mutation_list]
scores = [score for _, _, _, score in mutation_list]

# Plot the top20
plt.figure(figsize=(5, 5))
plt.bar(labels, scores, color='skyblue')
plt.xlabel("Mutation")
plt.ylabel("Count / -")
plt.xticks(rotation=90)
plt.xlim(-0.5,20.5)
plt.tight_layout()
plt.savefig(args.output, dpi=600)