#!/usr/bin/env python3
"""Turn per-residue amino acid probabilities into a PSI-BLAST ASCII PSSM.

The ESM and ProteinMPNN rules emit "<pos> <restype> <probability>" tables, but the
PROSS protocol reads a PSSM through SeqprofConsensus and FavorSequenceProfile, which
expect PSI-BLAST log-odds. Feeding probabilities straight in would break the design
gate: SeqprofConsensus keeps every residue type scoring >= 0, and a probability is
never negative, so all 20 amino acids would pass and the filter would do nothing.

Scores are therefore 2 * log2(p / background), the half-bit scale PSI-BLAST uses,
clamped because the neural models put far more mass on their favourite residue than
an alignment does -- ESM reaches p = 2e-6, which would otherwise score about -29 and
swamp the res_type_constraint bonus that is tuned for BLAST-sized numbers.
"""

import argparse
import math

# Robinson & Robinson background amino acid frequencies, as used by BLAST.
BACKGROUND = {
    "A": 0.07805, "R": 0.05129, "N": 0.04487, "D": 0.05364, "C": 0.01925,
    "Q": 0.04264, "E": 0.06295, "G": 0.07377, "H": 0.02199, "I": 0.05142,
    "L": 0.09019, "K": 0.05744, "M": 0.02243, "F": 0.03856, "P": 0.05203,
    "S": 0.07120, "T": 0.05841, "W": 0.01330, "Y": 0.03216, "V": 0.06441,
}

# Column order of a PSI-BLAST PSSM.
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

HEADER = (
    "Last position-specific scoring matrix computed, weighted observed percentages "
    "rounded down, information per position, and relative weight of gapless real "
    "matches to pseudocounts"
)


def read_weights(path):
    """Return {position: {one_letter_aa: probability}}."""
    profile = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            pos, restype, weight = line.split()
            aa = THREE_TO_ONE.get(restype.upper())
            if aa is None:
                raise ValueError(f"unknown residue type {restype!r} in {path}")
            profile.setdefault(int(pos), {})[aa] = float(weight)
    if not profile:
        raise ValueError(f"no probabilities found in {path}")
    return profile


def read_sequence(path):
    with open(path) as handle:
        return "".join(l.strip() for l in handle if not l.startswith(">"))


def log_odds(probability, aa, clamp):
    if probability <= 0.0:
        return -clamp
    score = 2.0 * math.log2(probability / BACKGROUND[aa])
    return max(-clamp, min(clamp, int(round(score))))


def write_pssm(profile, sequence, out_path, clamp):
    positions = sorted(profile)
    if len(positions) != len(sequence):
        raise ValueError(
            f"profile has {len(positions)} positions but the sequence has "
            f"{len(sequence)} residues -- they must describe the same chain"
        )

    with open(out_path, "w") as out:
        out.write("\n")
        out.write(HEADER + "\n")
        out.write(" " * 8 + "".join(f"{aa:>4}" for aa in AA_ORDER) * 2 + "\n")

        for pos, native in zip(positions, sequence):
            probs = profile[pos]
            missing = set(AA_ORDER) - set(probs)
            if missing:
                raise ValueError(
                    f"position {pos} is missing {''.join(sorted(missing))}"
                )
            scores = "".join(f"{log_odds(probs[aa], aa, clamp):>4}" for aa in AA_ORDER)
            # PSI-BLAST reports observed percentages alongside the scores; the design
            # protocol reads only the scores, so the probabilities stand in directly.
            percents = "".join(f"{int(round(probs[aa] * 100)):>4}" for aa in AA_ORDER)
            info = sum(
                p * math.log2(p / BACKGROUND[aa])
                for aa, p in probs.items() if p > 0
            )
            out.write(f"{pos:>5} {native}{scores}{percents} {info:>5.2f} {0.00:>4.2f}\n")

        out.write("\n")
        out.write(" " * 22 + "K         Lambda\n")
        out.write("PSI Ungapped         0.1639     0.3169\n")
        out.write("PSI Gapped           0.0502     0.2670\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--weights", required=True,
                        help="probability table written by the ESM/ProteinMPNN rules")
    parser.add_argument("-s", "--sequence", required=True,
                        help="FASTA of the same chain, for the native residue column")
    parser.add_argument("-o", "--out", required=True, help="PSSM to write")
    parser.add_argument("--clamp", type=int, default=10,
                        help="largest absolute log-odds score (default: 10)")
    args = parser.parse_args()

    profile = read_weights(args.weights)
    sequence = read_sequence(args.sequence)
    write_pssm(profile, sequence, args.out, args.clamp)


if __name__ == "__main__":
    main()
