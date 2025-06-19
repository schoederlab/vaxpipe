#!/bin/bash

mkdir -p output

for file in *.pdb; do
  $ROSETTA/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
  @ design.options \
  -parser:protocol run_mpnn_and_save.xml \
  -s $file \
  -parser:script_vars weights=${file}_mpnn_probs.weights \
  -out:path:all output \
  -overwrite
done
