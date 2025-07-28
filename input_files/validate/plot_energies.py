import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from argparse import ArgumentParser

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 2

parser = ArgumentParser(prog='Plot Frequencies for evaluation of proposed mutations')
parser.add_argument('-i1', '--input1', help='inputpath1', type=str)
parser.add_argument('-i2', '--input2', help='inputpath2', type=str)
parser.add_argument('-o', '--output', type=str, default='output.png', help='filename for the resulting output png. default: output.png')
args = parser.parse_args()


df1 = []
files = glob.glob(f'{args.input1}/score_design_*?.sc')
for f in sorted(files):
    df1.append(pd.read_csv(f, sep=r'\s+', header=1))

df2 = []
files = glob.glob(f'{args.input2}/score_design_*?.sc')
for f in sorted(files):
    df2.append(pd.read_csv(f, sep=r'\s+', header=1))

delta_total_score = []

for i in range(0, len(files)):
    delta_total_score.append(df1[i]['total_score'] - df2[i]['total_score'])

files = sorted(files)
ticklabels = [filename.split('/')[-1].replace('score_design_', '').replace('.sc', '')
              for filename in files]

lowest_delta_total_score = []

for i in range(0, len(files)):
    lowest_delta_total_score.append(df1[i]['total_score'].min() - df2[i]['total_score'].min())

filtered_lists = [[x for x in sublist if not np.isnan(x)] for sublist in delta_total_score]

plt.figure(figsize=(8,5))
plt.boxplot(filtered_lists)
plt.plot(range(1,1+len(lowest_delta_total_score)), lowest_delta_total_score, linestyle='None', marker='.', color='black', markerfacecolor='orange')
plt.ylabel('Δ total score / REU')
plt.xlim(0.5,21.5)
plt.xticks(range(1, len(ticklabels)+1), labels=ticklabels, rotation=90)
plt.xlabel('variant')
plt.hlines(y=0, xmin=0, xmax=22, color='black', linestyle='--')
plt.savefig(f'{args.output}', bbox_inches='tight')