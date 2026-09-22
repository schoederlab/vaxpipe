# vaxpipe

**vaxpipe** is a Snakemake-based pipeline designed to automate structural modeling of symmetric protein complexes using the [Rosetta](https://www.rosettacommons.org/) molecular modeling suite. 

It streamlines the process of the prediction of possible amino acid mutations using deep learning based methods, such as ProteinMPNN and ESM and database driven approaches, such as Rosetta FastDesign, with a focus on shape complementarity of protein interfaces.


### Operational Guidelines

- **Input PDB Structure**: The pipeline expects an input PDB file (e.g., `test/3ft7.pdb`) which represents a symmetric protein complex. **The file has to be placed inside the `workdir` configured in `config.yaml`**, and its basename (without `.pdb`) is the `samples` entry.
- **Rosetta**: All Rosetta steps run inside a singularity image, so no local Rosetta build is needed. Point `rosettadir` in `config.yaml` at the directory holding `rosetta_ml.sif`.
- **Snakemake**: The pipeline is orchestrated using Snakemake. There is a single `snakefile` for both local and cluster execution — they differ only in how Snakemake is invoked. While users should ideally be familiar with Snakemake, this documentation aims to provide sufficient guidance for basic local and HPC execution.
- **Configuration**: The `config.yaml` file must be updated with appropriate paths, the name of the pdb file, and input/output directories. This file is critical for the pipeline's correct operation. See [Configuration](#configuration) for every key.
- **Output Files**: Intermediate and final output files will be generated in a structured directory format, as defined in the `snakefile`. If the pipeline run successfully, several .png files will be present in the output directory, which represent both frequency and energy evaluation of tested mutations.
- **Resource Management**: For HPC execution, cluster submission parameters (partition, account, runtime, memory) are set in `profiles/slurm/config.yaml`; per-rule overrides live in the `resources:` blocks of the `snakefile`.

### Requirements

The same requirements apply to local and cluster runs, because every Rosetta call goes
through the singularity image either way.

- `Python ≥ 3.7, biopython, tqdm, matplotlib, pandas`
- [Snakemake ≥ 9.0](https://snakemake.readthedocs.io/en/v9.3.0/), plus `snakemake-executor-plugin-slurm` for cluster runs
- **singularity**
- The Rosetta image, pulled from Docker Hub with `singularity pull docker://rosettacommons/rosetta:ml-387` and placed as `rosetta_ml.sif` inside `rosettadir`.
- ESM model [download](https://git.iwe-lab.de/moritzertelt/ML_graphs/-/tree/main/tensorflow_graphs/ESM/esm2_t33_650M_UR50D). Currently, the pipeline just accepts this ESM model. **The model needs to be downloaded and copied into the repository path**, where the `run_esm` rule bind-mounts it into the container.
- Only needed for the PROSS branch (see [PROSS](#pross)): the BLAST+ package (https://blast.ncbi.nlm.nih.gov) and `hhblits`, which run on the host rather than in the container, plus the `UniRef30_2020_06` database. See [PROSS dependencies](#pross-dependencies) for the install and download commands.

### Repository Structure

```
vaxpipe/
├── snakefile # the workflow, used for both local and cluster runs
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
| `rosettadir` | Directory holding the Rosetta container, which is expected at `<rosettadir>/rosetta_ml.sif`. |
| `workdir` | Working directory. The input PDB has to live here; all results are written here. |
| `inputdir` | This repository's `input_files/` directory (RosettaScripts XMLs and helper scripts). |
| `uniref_db` | Prefix of the hhsuite `UniRef30_2020_06` database, used by `hhblits` to build the PROSS PSSM. |
| `samples` | Basename of the input PDB (`3ft7` for `3ft7.pdb`). |
| `iterations` | Designs generated per branch (esm, pmpnn, indes). |
| `mutations` | Most frequent mutations carried into the design/control validation. |
| `pross_temps` | PROSS delta-score thresholds. These must match `delta_filter_thresholds` of the `FilterScan` filter in `input_files/pross/filter/filterscan.xml`, because that filter derives the resfile names from them. |

`iterations` and `mutations` are what a local trial run and a production run on the
cluster mainly differ in. The committed values (200 and 20) are meant for the cluster;
for a quick local check override them on the command line rather than editing the file:

```bash
snakemake --cores 4 --config iterations=5 mutations=5
```

---

### Installation

#### Installation using conda

generate snakemake conda environment

```
conda create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
pip install snakemake-executor-plugin-slurm
pip install biopython tqdm matplotlib pandas

```

#### ESM model

The `run_esm` rule bind-mounts the model from the directory *above* `workdir`, so the
model has to end up at `<workdir>/../esm2_t33_650M_UR50D`. With the committed
`config.yaml` (`workdir: <repo>/test`) that is the repository root:

```
<repo>/esm2_t33_650M_UR50D/
├── saved_model.pb
├── keras_metadata.pb
└── variables/
```

A symlink works as well as a copy, and `.gitignore` already excludes the directory.

#### PROSS dependencies

Only needed if you enable the `#pross` block in `rule all` (see [PROSS](#pross)).
`hhblits` and `psiblast` run on the *host*, not in the Rosetta container, so they go
into the same conda environment:

```
conda activate snakemake
conda install -c conda-forge -c bioconda blast hhsuite
```

Then fetch the UniRef30 database (~50 GB compressed, ~86 GB unpacked) and point
`uniref_db` at the file prefix, not at the directory:

```
mkdir -p input_files/UniRef30_2020_06 && cd input_files/UniRef30_2020_06
wget -c https://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
tar -xzf UniRef30_2020_06_hhsuite.tar.gz
rm UniRef30_2020_06_hhsuite.tar.gz
```

This unpacks `UniRef30_2020_06_*.ffdata` / `*.ffindex` files, which `hhblits -d`
addresses through their shared prefix, hence the doubled name in
`uniref_db: .../UniRef30_2020_06/UniRef30_2020_06`.

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
2.  **Execute Locally**: Run the pipeline on your local machine, utilizing `X` CPU cores. Singularity still has to be available, since every Rosetta step runs in the container:
    ```bash
    snakemake --cores X
    ```
3.  **Execute on HPC Cluster**: For cluster execution, use the SLURM profile, which submits each job through Snakemake’s executor interface (replacement for the deprecated `--cluster` flag):
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

### PROSS

The PROSS rules (`generate_PSSM_and_constraints`, `filterscan`, `pross_design`,
`pross_design_wt`) are implemented but their targets are commented out in `rule all`,
so they do not run by default. They are the only part of the pipeline that needs
`hhblits`, BLAST+ and the UniRef30 database on the host. To enable them, uncomment the
`#pross` block in `rule all` and set `uniref_db` in `config.yaml`.

Note that `filterscan` runs one Rosetta job per residue and writes the resfiles for all
`pross_temps` thresholds in that single pass, which is why it is one Snakemake job
rather than one per threshold.

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
