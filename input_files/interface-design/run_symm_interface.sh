symfile=prefusion_relaxed.symm
input=prefusion_relaxed_INPUT.pdb

singularity run ../rosetta_ml.sif rosetta_scripts @design.options \
	-parser:protocol sym_design.xml \
	-s $input \
	-parser:script_vars symfile=$symfile \
	-out:file:scorefile des_relax-design.sc \
	-out:suffix _des_relax-design \
	-overwrite > des_design.log
