#!/bin/sh
# Emit one CoordinateConstraint per CA atom, all relative to the first CA of the structure.
# Two passes over the file: the first finds the root atom, the second writes the constraints.
pdb=$1

awk '
	FNR == NR {
		# first pass: remember the first ATOM CA as the constraint root
		if (root == "" && $1 == "ATOM" && $3 == "CA") {
			root = "CA " trim(substr($0, 23, 4)) substr($0, 22, 1)
		}
		next
	}
	($1 == "ATOM" || $1 == "HETATM") && $3 == "CA" {
		res_num = trim(substr($0, 23, 4))
		chain = substr($0, 22, 1)
		x = trim(substr($0, 31, 8))
		y = trim(substr($0, 39, 8))
		z = trim(substr($0, 47, 8))
		# harmonic constraint, mean=0 std=1
		print "CoordinateConstraint", $3, res_num chain, root, x, y, z, "HARMONIC 0 1"
	}
	function trim(s) {
		gsub(/^[ \t]+|[ \t]+$/, "", s)
		return s
	}
' "$pdb" "$pdb"
