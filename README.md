# vaxpipe

**vaxpipe** is a Snakemake-based pipeline designed to automate structural modeling of symmetric protein complexes using the [Rosetta](https://www.rosettacommons.org/) molecular modeling suite. 

It streamlines the process of generating symmetry definition files, running symmetric relaxation on input PDB structures, and prediction of possible amino acid mutations using deep learning based methods, such as ProteinMPNN and ESM and database driven approaches, such as Rosetta FastDesign, with a focus on shape complementarity of protein interfaces.

---

## Purpose

The goal of `vaxpipe` is to support **vaccine design** workflows by:

- Generating symmetry definitions using Rosetta's `make_symmdef_file.pl` to speed up the calculations
- Performing energy minimization (`relax`) on symmetric oligomers
- predicting amino acid probabilities using deep learning based methods (P-MPNN, ESM)
- designing novel sequences which are further evaluated using Rosetta
- Ensuring reproducibility and modularity using [Snakemake](https://snakemake.readthedocs.io/)

## Requirements

- Python ≥ 3.7, Biopython, tqdm
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- Rosetta compiled with pytorch and tensorflow libraries (`relax`, symmetry tools)

## Repository Structure

```
vaxpipe/
├── snakefile # the main workflow
├── input_files/ # containes the necessary input .xml files to run the rosetta jobs
├── test/ # containes the test input pdb for validation
├── validate/ # containes the necessary input for validating proposed mutations
└── output/ # Relaxed models and symmetry files

```
---

## Installation

generate snakemake conda environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
pip install biopython tqdm
```

change the `config.yaml` file with the corresponding paths

## Execution

```
snakemake -j 12 # executing vaxpipe using 12 cpu cores
```