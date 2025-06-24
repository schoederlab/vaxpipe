#!/bin/bash
# first flag is name
# second flag is number of relax runs

#generate-sym
$ROSETTA/main/source/src/apps/public/symmetry/make_symmdef_file.pl -p ${1}.pdb -a A -i B > prefusion.symm

#run sym relax
mkdir -p sym-relax
cd sym-relax
$ROSETTA/main/source/bin/relax.pytorchtensorflow.linuxgccrelease -s ../${1}_INPUT.pdb \
 -constrain_relax_to_start_coords \
 -nstruct ${2} -out:suffix _relax \
 -multiple_processes_writing_to_one_directory \
 -out:prefix relax_ \
 -scorefile ${1}_relax.sc \
 -symmetry_definition ../prefusion.symm
cd ..
