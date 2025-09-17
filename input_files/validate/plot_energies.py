import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
from argparse import ArgumentParser

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 2

parser = ArgumentParser(prog='Plot Frequencies for evaluation of proposed mutations')
parser.add_argument('-i1', '--input1', help='inputpath1', type=str)
parser.add_argument('-i2', '--input2', help='inputpath2', type=str)
parser.add_argument('-o', '--output', type=str, default='output.png',
                    help='filename for the resulting output png. default: output.png')
args = parser.parse_args()

# Helper: extract numeric part of filename for sorting and labeling
def numeric_key(fname):
    match = re.search(r'(\d+)\.sc$', fname)
    return int(match.group(1)) if match else 0

# Load files with numeric sorting
files1 = sorted(glob.glob(f'{args.input1}/*.sc'), key=numeric_key)
files2 = sorted(glob.glob(f'{args.input2}/*.sc'), key=numeric_key)

# Read data
df1 = [pd.read_csv(f, sep=r'\s+', header=1) for f in files1]
df2 = [pd.read_csv(f, sep=r'\s+', header=1) for f in files2]

# Compute delta scores
delta_total_score = [df1[i]['total_score'] - df2[i]['total_score'] for i in range(len(files1))]
lowest_delta_total_score = [df1[i]['total_score'].min() - df2[i]['total_score'].min()
                            for i in range(len(files1))]

# Filter out NaNs
filtered_lists = [[x for x in sublist if not np.isnan(x)] for sublist in delta_total_score]

# Variant labels = just the numbers from filenames
ticklabels = [str(numeric_key(f)) for f in files1]

# Plot
plt.figure(figsize=(8,5))
plt.boxplot(filtered_lists)
plt.plot(range(1, 1+len(lowest_delta_total_score)),
         lowest_delta_total_score, linestyle='None', marker='.',
         color='black', markerfacecolor='orange')

plt.ylabel('Δ total score / REU')
plt.xlabel('variant')
plt.xticks(range(1, len(ticklabels)+1), labels=ticklabels, rotation=90)

plt.xlim(0.5, len(ticklabels) + 0.5)
plt.hlines(y=0, xmin=0, xmax=len(ticklabels)+1, color='black', linestyle='--')

plt.savefig(args.output, bbox_inches='tight')