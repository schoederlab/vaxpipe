#!/bin/bash

# Define the list of mutations
mutations=("1062A_ALA" "1048A_ALA" "1085A_ARG" "799A_ASN" "915A_ASN" "1004A_LEU" "973A_ASN" "1104A_ASN" "952A_ALA" "962A_GLN")

# Design or control for protocol
protocol="control"
outdir="output-control"
nstruct=20

# Create the output directory if it doesn't exist
mkdir -p ${outdir}

# Loop through each mutation
for mut in "${mutations[@]}"; do
    # Split the mutation string into position and amino acid
    mutpos=${mut%_*}
    mutaa=${mut#*_}

    # Run the command
    $ROSETTA/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
        -parser:protocol design.v02.xml \
        -parser:script_vars mutpos=${mutpos} mut_aa=${mutaa} protocol=${protocol} symfile=prefusion_relaxed.symm \
        -in:file:s prefusion_relaxed_INPUT.pdb \
        -corrections:beta_nov16 \
        -out:path:all ${outdir} \
        -overwrite \
        -out:suffix _design_${mutpos}_${mutaa} \
        -nstruct ${nstruct}
done
