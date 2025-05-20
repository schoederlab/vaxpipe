#!/bin/bash
#SBATCH --job-name=rosetta_design
#SBATCH --output=logs/rosetta_%A_%a.out
#SBATCH --error=logs/rosetta_%A_%a.err
#SBATCH --partition=barnard
#SBATCH --array=0-20
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00


# Define the list of mutations
mutations=("752A_ALA" "771A_TRP" "52A_VAL" "1037A_GLY" "461A_ALA" "1075A_TRP" "871A_ASN" "50A_VAL" "1065A_PHE" "335A_THR" "1136A_ALA" "405A_ILE" "1062A_ALA" "666A_LEU" "444A_ASN" "448A_ALA" "928A_TRP" "1094A_SER" "808A_GLN" "844A_VAL" "661A_ASN")

# Get the mutation for this task
mutation=${mutations[$SLURM_ARRAY_TASK_ID]}
mutpos=${mutation%_*}
mutaa=${mutation#*_}

# Design or control for protocol
protocol="design"
outdir="output-design"
nstruct=100

# Create the output directory if it doesn't exist
mkdir -p ${outdir}

# Run the command
singularity run -B /data/horse/ws/jari462g-ml-rosetta/interface-design/scoring ../../rosetta_ml.sif rosetta_scripts \
    -parser:protocol design.v02.xml \
    -parser:script_vars mutpos=${mutpos} mut_aa=${mutaa} protocol=${protocol} symfile=prefusion_relaxed.symm \
    -in:file:s prefusion_relaxed_INPUT.pdb \
    -corrections:beta_nov16 \
    -out:path:all ${outdir} \
    -overwrite \
    -out:suffix _design_${mutpos}_${mutaa} \
    -nstruct ${nstruct}






