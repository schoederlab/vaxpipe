from __future__ import print_function
import os, sys
from argparse import ArgumentParser, BooleanOptionalAction
from argparse import RawTextHelpFormatter
from math import floor
import time

try:
    from tqdm import tqdm
    PROGRESS = True
except ImportError:
    PROGRESS = False
    print('[INFO] No tqdm available. Continue without progress.')

try:
    from Bio import PDB
except ImportError:
    print('[ERROR] The biopython module is needed. Exit.')
    sys.exit()

alpha_1 = list("ARNDCQEGHILKMFPSTWYV-")
states = len(alpha_1)
alpha_3 = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
           'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','X']
three_to_one = dict(zip(alpha_3, alpha_1))

def check_existense(file) -> None:
    '''Simple check file existence. Exit if not.'''
    if not os.path.exists(file):
        print(f'[Error] The provided file or path "{file}" dont exists. Exit.')
        exit(1)


def get_files_from_path(path:str, filesuffix:str):
    '''Return all files (with path) in a given path containing the filesuffix.'''
    return [os.path.join(path, file) for file in os.listdir(path) if file.endswith(filesuffix)]


def three_to_one_(resname):
    try:
        aa = three_to_one[resname]
        return aa
    except KeyError:
        return 'X'


def get_aa_seq_from_pdb(filename, chain_id:str=None, gaps:bool=False):
    parser = PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(id=filename, file=filename)
    three_code = []
    for model in struct:
        for chain in model:
            # skip chain if we only interested in a single chain
            if chain_id and chain.get_id() != chain_id:
                continue
            
            # convert to list to access elements by index
            chain = chain.get_unpacked_list()
            for i, residue in enumerate(chain):
                if gaps:
                    try:
                        last_j = chain[i-1].get_id()[1]
                        curr_j = chain[i].get_id()[1]
                        if curr_j - last_j > 1:
                            three_code.append('-' * (curr_j - last_j - 1))
                        three_code.append(three_to_one_(residue.resname))
                    except IndexError:
                        three_code.append(three_to_one_(residue.resname))
                else:
                    three_code.append(three_to_one_(residue.resname))
                    
    if chain_id and len(three_code) == 0:
        print(f'[WARNING] Chain id "{chain_id}" not found in structure.')
        return -1
    
    return ''.join(three_code)


def _write_fasta(result_seq:list, output_file:str, linebreak:bool, suffix:str):
    if suffix:
        output_file = output_file.replace('.fasta', '') + '_' + str(round(time.time())) + '_' + suffix + '.fasta'
    else:
        output_file = output_file
         
    with open(output_file, 'w') as fasta:
        for name, seq in result_seq:
            fasta.write(f'>{name}\n')
            if linebreak:
                r = floor(len(seq)/80)
                for j in range(r+1):
                    fasta.write(seq[ (j*80) : ((j+1)*80)] + '\n')
            else:
                fasta.write(seq + '\n')


def collect_pdb_files(paths: list, filesuffix: str = ".pdb"):
    """Collect PDB files from a list of files or directories."""
    files = []
    for path in paths:
        if os.path.isfile(path):
            if path.endswith(filesuffix):
                files.append(path)
            else:
                print(f'[WARNING] File "{path}" does not end with {filesuffix}, skipping.')
        elif os.path.isdir(path):
            files.extend(get_files_from_path(path=path, filesuffix=filesuffix))
        else:
            print(f'[ERROR] Path "{path}" is neither a file nor a directory. Skipping.')
    return files


def main(args):
    input_path = args.pdbpath
    output_file = args.output
    
    chain_id = args.chain if args.chain else None
    gaps = True if args.gaps else False
    linebreak = True if args.linebreak else False
    batch = args.batch
    
    # new logic
    files = collect_pdb_files(input_path, filesuffix=".pdb")
    
    if PROGRESS and len(files) > 1:
        bar = tqdm
    else:
        def f(x, **kwargs): return x
        bar = f
    
    print(f'[INFO] Reading files from "{input_path}".')
    result_seq = []
    for i in bar(range(len(files)), desc=f'[INFO] Extract fasta from "{input_path}"', dynamic_ncols=True):
        basename = files[i].split('/')[-1].replace('.pdb', '')
        basename = basename + f'_chain_{chain_id}' if chain_id else basename
        
        aa_seq = get_aa_seq_from_pdb(filename=files[i], chain_id=chain_id, gaps=gaps)
        if aa_seq == -1:
            exit(1)
        if aa_seq:
            result_seq.append((basename, aa_seq))
        
        if len(result_seq)%batch == 0:
            suffix = 'batch_' + str(len(result_seq))
            _write_fasta(result_seq, output_file, linebreak, suffix)
    
    print(f'[INFO] Writing sequences to "{output_file}".')
    _write_fasta(result_seq, output_file, linebreak, suffix='')
    print('[INFO] Done.')
    return


# =========================================================
if __name__ == '__main__':
    parser = ArgumentParser(prog='Create fasta for pdbs.', 
                            description='Create a fasta file containing the sequences of ' +
                            'all pdb files in the given directory or explicit file list.\n' +
                            'Non standard amino acids / Ligands / HETATM ' +
                            'are replaced by "X".\n\n', 
                            formatter_class=RawTextHelpFormatter)
    
    parser.add_argument('-p', '--path', type=str, dest='pdbpath', nargs='+', required=True,
                        help='List of pdb files and/or directories to create fasta file.\n' +
                        'All files with the suffix ".pdb" will be processed.\n' +
                        'Single pdbs are also allowed.')
    parser.add_argument('-c', '--chain', type=str, dest='chain', required=False,
                        help='Only get the sequence of a specific chain for all pdbs.')
    parser.add_argument('--gaps', dest='gaps', action=BooleanOptionalAction,
                        help='Insert the gaps in the structures based on pdb numbering.\n'+
                        'Default: no gaps')
    parser.add_argument('--linebreak', dest='linebreak', action=BooleanOptionalAction,
                        help='Write fasta file with linebreak for sequences at 80 column width. ' +
                        'Default: no linebreaks')
    parser.add_argument('--batch', dest='batch', type=int, default=10000,
                        help='Write fasta file in batches. Recommended for paths with\n'
                        + 'large number of pdbs. Batch files will be named OUTPUTNAME_(time)_batch_(N), where\n'
                        + 'N is the number of total sequence in file. Default: 10,000')
    parser.add_argument('-o', '--out', type=str, dest='output', required=False, default=None,
                        help='Name for output file.\nDefault: Name of the path/pdbfile with suffix ".fasta".')
    args = parser.parse_args()

    # Create output name
    if not args.output:
        if len(args.pdbpath) == 1 and os.path.isfile(args.pdbpath[0]):
            args.output = os.path.basename(args.pdbpath[0]).replace('.pdb', '.fasta')
        else:
            p = os.path.dirname(os.path.join(args.pdbpath[0], ''))
            args.output = p + '.fasta' if p != '.' else os.getcwd().split('/')[-1] + '.fasta'
    
    main(args)