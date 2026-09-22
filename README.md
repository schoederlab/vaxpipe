# vaxpipe

**vaxpipe** is a Snakemake-based pipeline designed to automate structural modeling of symmetric protein complexes using the [Rosetta](https://www.rosettacommons.org/) molecular modeling suite. 

It streamlines the process of the prediction of possible amino acid mutations using deep learning based methods, such as ProteinMPNN and ESM and database driven approaches, such as Rosetta FastDesign, with a focus on shape complementarity of protein interfaces.


### Operational Guidelines

- **Input PDB Structure**: The pipeline expects an input PDB file (e.g., `test/3ft7.pdb`) which represents a symmetric protein complex. **The file has to be placed inside the `workdir` configured in `config.yaml`**, and its basename (without `.pdb`) is the `samples` entry.
- **Rosetta Installation**: A functional Rosetta installation compiled with PyTorch and TensorFlow libraries is required for the deep learning-based predictions. For cluster usage, a singularity image can be downloaded from Docker Hub (see details below). Ensure that the `ROSETTA` path is correctly specified in the `config.yaml`.
- **Snakemake**: The pipeline is orchestrated using Snakemake. While users should ideally be familiar with Snakemake, this documentation aims to provide sufficient guidance for basic local and HPC execution. Detailed instructions are provided within the relevant sections.
- **Configuration**: The `config.yaml` file must be updated with appropriate paths to Rosetta binaries, the name of the pdb file, and input/output directories. This file is critical for the pipeline's correct operation. See [Configuration](#configuration) for every key.
- **Output Files**: Intermediate and final output files will be generated in a structured directory format, as defined in the `snakefile`. If the pipeline run successfully, several .png files will be present in the output directory, which represent both frequency and energy evaluation of tested mutations.
- **Resource Management**: For HPC execution, cluster submission parameters (partition, account, runtime, memory) are set in `profiles/slurm/config.yaml`; per-rule overrides live in the `resources:` blocks of `snakefile-hpc`.

### Local Execution Requirements

- `Python ≥ 3.7, biopython, tqdm, matplotlib, pandas`
- [Snakemake ≥ 9.0](https://snakemake.readthedocs.io/en/v9.3.0/)
- the BLAST+ software package: (https://blast.ncbi.nlm.nih.gov)
- Rosetta compiled with pytorch and tensorflow libraries. A detailed information on how to compile Rosetta with pytorch and tensorflow support can be found [here](https://docs.rosettacommons.org/docs/latest/build_documentation/Building-Rosetta-with-TensorFlow-and-Torch)
- ESM model (will be downloaded automatically)
- for the PROSS protocol, PSSMs are generated, which needs the `UniRef30_2020_06` database (`wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz`). Point `uniref_db` in `config.yaml` at the database prefix.


### HPC Requirements
- `Python ≥ 3.7, biopython, tqdm, matplotlib, pandas`
- Snakemake ≥ 9.0
- **singularity**
- We provide a snakemake file that relies on a Rosetta docker/singularity image. The image is available on Docker Hub and can be pulled using singularity: `singularity pull docker://rosettacommons/rosetta:ml-387` (used for cluster execution). `snakefile-hpc` expects it as `rosetta_ml.sif` inside `rosettadir`.
- ESM model [download](https://git.iwe-lab.de/moritzertelt/ML_graphs/-/tree/main/tensorflow_graphs/ESM/esm2_t33_650M_UR50D). Currently, the pipeline just accepts this ESM model. **The model needs to be downloaded and copied into the repository path.**
- for the PROSS protocol, PSSMs are generated, which needs the `UniRef30_2020_06` database (`wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz`). Point `uniref_db` in `config.yaml` at the database prefix.

### Repository Structure

```
vaxpipe/
├── snakefile # the main workflow
├── snakefile-hpc # the main workflow optimized for the HPC cluster
├── config.yaml # containes path information
├── input_files/ # containes the necessary input .xml files to run the rosetta jobs
├── profiles/slurm/ # containes the config.yaml file for cluster execution
└── test/ # containes a test input pdb for validation/benchmarking
```
---

### Configuration

All paths live in `config.yaml`:

| key | meaning |
| --- | --- |
| `rosettadir` | Rosetta root. Locally the binaries are expected under `<rosettadir>/main/source/bin/`, on the cluster the container is expected at `<rosettadir>/rosetta_ml.sif`. |
| `workdir` | Working directory. The input PDB has to live here; all results are written here. |
| `inputdir` | This repository's `input_files/` directory (RosettaScripts XMLs and helper scripts). |
| `uniref_db` | Prefix of the hhsuite `UniRef30_2020_06` database, used by `hhblits` to build the PROSS PSSM. |
| `samples` | Basename of the input PDB (`3ft7` for `3ft7.pdb`). |
| `pross_temps` | PROSS delta-score thresholds. These must match `delta_filter_thresholds` of the `FilterScan` filter in `input_files/pross/filter/filterscan.xml`, because that filter derives the resfile names from them. |

The number of design iterations and the number of mutations that are carried forward
into validation are set at the top of the snakefiles (`ITERATIONS` and `MUTATIONS`).
The local `snakefile` uses small values (5 each) so a run finishes on a workstation;
`snakefile-hpc` uses 200 iterations and 20 mutations.

---

### Installation

#### Installation using conda

generate snakemake conda environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
conda install bioconda::blast
pip install snakemake-executor-plugin-slurm
pip install biopython tqdm matplotlib pandas

```

**change the `config.yaml` file with the corresponding paths**

---

### Understanding and Executing the Pipeline with Snakemake

`vaxpipe` is built with Snakemake, a powerful workflow management system that helps you create reproducible and scalable data analysis pipelines. Here's a brief guide to get you started:

#### Core Concepts:
- **Rules**: Define how output files are generated from input files. Each rule specifies a command to be run.
- **Wildcards**: Allow rules to be applied to many files (e.g., `data/{sample}.txt` where `{sample}` is a wildcard).
- **DAG (Directed Acyclic Graph)**: Snakemake automatically builds a DAG of jobs, ensuring that all dependencies are met before a rule is executed.
- **Execution**: Snakemake determines which rules need to be run based on the desired output files and the current state of your input files.

#### Essential Commands:

1.  **Dry Run**: Always start with a dry run to see what Snakemake plans to do without actually executing any commands (locally):
    ```bash
    snakemake --dry-run 
    # or a shorter version
    snakemake -n
    ```
2.  **Execute Locally**: Run the pipeline on your local machine, utilizing `X` CPU cores and opting into Conda-based software deployment (replacement for the deprecated `--use-conda` flag):
    ```bash
    snakemake --cores X --software-deployment-method conda
    ```
3.  **Execute on HPC Cluster**: For cluster execution, use the HPC-optimized snakefile together with Snakemake’s executor interface (replacement for the deprecated `--cluster` flag). See the `hpc execution` section below for a detailed SLURM example:
    ```bash
    snakemake --profile ./profiles/slurm
    ```
For more in-depth information, please refer to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/).

---

### Pipeline Stages

1. **Preprocessing** — `score_jd2` cleans and renumbers the input PDB, `make_symmdef_file.pl`
   derives the symmetry definition, and a symmetric `relax` produces the structure that
   every design branch starts from (`relax_<sample>_new_0001_INPUT.pdb`).
2. **Sampling** — three independent branches produce `ITERATIONS` designs each:
   `esm` (ESM-2 probabilities), `pmpnn` (ProteinMPNN probabilities) and `indes`
   (Rosetta symmetric interface `FastDesign`).
3. **PROSS** — PSSM plus coordinate constraints, a `FilterScan` over every residue, and
   design/WT scoring at each threshold in `pross_temps`.
4. **Analysis** — sequences are extracted from all designs, mutation frequencies are
   counted against the wild type, and the `MUTATIONS` most frequent substitutions are
   each rebuilt twice (`design` with the mutation, `control` without it).
5. **Plotting** — a mutation-frequency bar chart and a Δ total score box plot per branch.

Rosetta names its output structures `<prefix><input>_<suffix>_00XX.pdb`, so the rule
outputs carry those `_0001` / `_00XX` infixes. When changing `-out:prefix`,
`-out:suffix` or `-nstruct`, the `output:` patterns have to be adjusted to match.

### Output Layout

```
<workdir>/
├── <sample>.pdb                          # your input
├── <sample>_clean_0001.pdb, <sample>.symm, <sample>_new.pdb
├── relax_<sample>_new_0001_INPUT.pdb     # starting point of all design branches
├── <sample>_2.symm
├── <sample>_WT.fasta                     # reference sequence
├── logs/                                 # one log file per job
├── esm/ | pmpnn/ | indes/
│   ├── <variant>_relax_<sample>_new_0001_INPUT_<iter>_0001.pdb
│   ├── <sample>_<variant>.fasta
│   ├── <sample>_<variant>_frequency.png / .csv
│   ├── <sample>_energydifference_<variant>.png
│   └── <sample>_<variant>/
│       ├── <n>.txt                       # selected mutation, e.g. "45A_LEU"
│       └── design/ | control/
│           ├── <n>.sc                    # score file that is compared
│           └── <n>_decoys/               # the 20 structures behind it
└── pross/
    ├── <sample>.hhr | .a3m | .psi | .pssm | .cst
    ├── <sample>_resfiles_pross/designable_aa_resfile.<temp>
    └── <sample>_pross_design_<temp>.sc | <sample>_pross_wt_<temp>.sc
```

Every rule writes its stdout/stderr to `<workdir>/logs/<rule>_<wildcards>.log`. If a
job fails, that file is the first place to look.
