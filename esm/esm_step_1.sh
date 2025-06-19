$ROSETTA/main/source/bin/rosetta_scripts.pytorchtensorflow.linuxgccrelease \
    -parser:protocol run_esm_and_save.xml \
    -s "prefusion.pdb" \
    -nstruct 200 \
    -ex1 \
    -ex2 \
    -beta \
    -auto_download \
    -multiple_processes_writing_to_one_directory 