#!/bin/bash

mkdir -p output2

$ROSETTA/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
  -parser:protocol sample_mutations.xml \
  -parser:script_vars sym=prefusion_relaxed.symm \
  -s prefusion_relaxed_INPUT.pdb \
  -out:path:all output2 \
  -nstruct 200 \
  -beta \
  -overwrite
