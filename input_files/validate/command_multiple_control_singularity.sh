#!/bin/bash

# Define the list of mutations
source "$1"

# Design or control for protocol
protocol="control"
outdir="output-control"
nstruct=100

# Create the output directory if it doesn't exist
mkdir -p "$2"

# Loop through each mutation
for mut in "${mutations[@]}"; do
    # Split the mutation string into position and amino acid
    mutpos=${mut%_*}
    mutaa=${mut#*_}

    # Run the command
    singularity run -B "$6" "$7"/rosetta_ml.sif rosetta_scripts \
        -parser:protocol "$5" \
        -parser:script_vars mutpos=${mutpos} mut_aa=${mutaa} protocol=${protocol} symfile="$3" \
        -in:file:s "$4" \
        -corrections:beta_nov16 \
        -out:path:all "$2" \
        -overwrite \
        -out:suffix _design_${mutpos}_${mutaa} \
        -nstruct ${nstruct} \
        -beta
done
