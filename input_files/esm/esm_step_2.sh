$ROSETTA/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
    -parser:protocol sample_mutations.xml \
    -s "prefusion_relaxed_INPUT.pdb" \
    -parser:script_vars sym="prefusion_relaxed.symm" \
    -out:prefix esm_ \
    -ex1 \
    -ex2 \
    -nstruct 200 \
    -beta \
    -overwrite
