#!/usr/bin/env python
''' puzzleGen is application written for Phylo to process raw sequence data 
    and to discover puzzles from their alignment.
'''
__author__ = "Rudolf Lam"
__email__ = "rudolf.lam@mail.mcgill.ca"
__status__ = "Production"

import Bio
from Bio import *
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Seq import *
from Bio.SeqRecord import *
from Bio.Alphabet import *
from Bio.Align import *

import numpy as np
import logging
import json

FORMATS = set(['json', 'fasta', 'embl', 'fastq', 'fastq-solexa', 'fastq-illumina', 'genbank', 'imgt', 'phd', 'sff', 'tab', 'qual'])
ALIGNMENTS = set(['clustalw', 'muscle', 'tcoffee'])

def writeFasta(sequences, fasta_name="out.fasta"):
    with open(fasta_name, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    return fasta_name

def makeRecords(seq_name_pairs):
    return [SeqRecord(Seq(seq, generic_dna), id=name) for seq, name in seq_name_pairs]

def readJSON(handle):
    import json
    data = json.load(handle)
    name_seq_pairs = [(name.strip(),seq.strip()) for name, seq in data.items()]
    records = makeRecords(name_seq_pairs)
    return records

def readSequences(handle, fileformat='fasta'):
    if fileformat == 'json':
        return readJSON(handle)
    return list(SeqIO.parse(handle, fileformat))

def standardizeSequences(seqs):
    maxLength = max([len(seq) for seq in seqs])
    return [seq + '-' * (maxLength - len(seq)) for seq in seqs]

def clean(msa):
    m,n=len(msa), len(msa[0])
    start,end = 0,n
    for i in range(len(msa[0])):
        if msa[:,i]!='-'*m:
            start = i
            break
    for i in range(len(msa[0])):
        if msa[:,-i]!='-'*m:
            end = -i 
            break
    return msa[:,start:end+1] if end!=0 else msa[:,start:]

def align(seqs, algorithm='clustalw'):
    msa = MultipleSeqAlignment(seqs)
    AlignIO.write(msa, "temp.fasta", "fasta")
    alignment = ''
    if algorithm=='clustalw':
        from Bio.Align.Applications import ClustalwCommandline
        cline = ClustalwCommandline("clustalw", infile="temp.fasta")
        result = cline()
        alignment = AlignIO.read("temp.aln", "clustal")
    elif algorithm=='muscle':
        from Bio.Align.Applications import MuscleCommandline
        cline = MuscleCommandline('./muscle3.8.31_i86linux64', input="temp.fasta", out="temp.aln", clw=True)
        result = cline()
        alignment = AlignIO.read("temp.aln", "clustal")
    elif algorithm=='tcoffee':
        from Bio.Align.Applications import TCoffeeCommandline
        cline = TCoffeeCommandline(infile="temp.fasta", output="clustalw", outfile="temp.aln")
        result = cline()
        alignment = AlignIO.read("temp.aln", "clustal")
    return alignment

def getStretches(seq):
    isGap = False
    gapCount = 0
    for c in seq:
        if c is '-':
            if not isGap:
                isGap = True
        else:
            if isGap:
                isGap = False
                gapCount += 1
    if isGap:
        gapCount+=1
    return gapCount
def meanSquareDiffLength(seqs):
    ungapped_lengths = [len(seq.seq.ungap('-'))for seq in seqs]
    ul = np.array(ungapped_lengths)
    m = np.mean(ul)
    return np.sqrt(np.mean((ul-m)**2))
def shannon_entropy(list_input):
    import math
    unique_base = set(list_input)                           # Get only the unique bases in a column
    #unique_base = unique_base.discard("-")
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base)                        # Number of residues of type i
        P_i = n_i/float(M)                                  # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    #print sh_entropy
    return sh_entropy
def treeLength(treefile):
    tree = Phylo.read(treefile, "newick")
    return tree.total_branch_length()
def pairwiseDiff(seqs):
    total = 0
    for i,seq in enumerate(seqs):
        for seq2 in seqs[i+1:]:
            total += (sum(c1!=c2 for c1,c2 in zip(seq.seq,seq2.seq)))
    return total
def extractFeatures(seqs):
    stretchCount = sum([getStretches(str(seq.seq)) for seq in seqs])
    msd = meanSquareDiffLength(seqs)
    tree = treeLength('temp.dnd')
    pairwise_difference = pairwiseDiff(seqs)
    entropy = shannon_entropy([str(seq.seq) for seq in seqs])
    return {'stretches': stretchCount,
            'num_seqs':len(seqs),
            'mean_square_difference':msd,
            'tree_description': tree,
            'mismatches': pairwise_difference,
            'entropy' : entropy}

def computeCol(alignment,i, mismatches={}, matches={}, gaps={}):
    if (i,i) in mismatches.keys():
        return
    col = alignment[:,i]
    nogaps = [c for c in col if c!='-']
    n_gaps = len(col) - len(nogaps)
    most_common = max(nogaps, key=lambda char: nogaps.count(char))
    n_mismatch = sum([1 for c in col if c!=most_common])
    n_match = sum([1 for c in col if c==most_common])
    mismatches[(i,i)] = n_mismatch
    matches[(i,i)] = n_match
    gaps[(i,i)] = n_gaps

def compute(alignment, i, j, mismatches={}, matches={}, gaps={}):
    if (i,j) in mismatches.keys():
        return
    if i == j:
        computeCol(alignment,i, mismatches, matches, gaps)
    elif i < j:
        compute(alignment,i,j-1, mismatches, matches, gaps)
        computeCol(alignment,j, mismatches, matches, gaps)
        mismatches[(i,j)] = mismatches[(i,j-1)] + mismatches[(j,j)]
        matches[(i,j)] = matches[(i,j-1)] + matches[(j,j)]
        gaps[(i,j)] = gaps[(i,j-1)] + gaps[(j,j)]

def interestingness(alignment, i, j, scoring, mismatches={}, matches={}, gaps={}):
    compute(alignment, i, j, mismatches, matches, gaps)
    return scoring(mismatches[(i,j)], matches[(i,j)], gaps[(i,j)])

def R(n_mismatch, n_match, n_gap):
    return n_mismatch/ (n_mismatch+ n_match+ n_gap)

def getBlock(alignment, i, n):
    return alignment[:,i:i+n]

def windowInterestingness(alignment, size, scoring, mismatches={}, matches={}, gaps={}):
    ''' Return the interestingness of all subsequences of a certain size '''
    n = len(alignment[0])
    return [interestingness(alignment, i, i+size-1, scoring, mismatches, matches, gaps) for i in range(n-size+1)]

def interestingWindows(alignment, min_size, max_size, scoring, mismatches={}, matches={}, gaps={}):
    interestingness = {}
    for i in range(min_size, max_size+1):
        interestingness[i] = windowInterestingness(alignment, i, scoring, mismatches={}, matches={}, gaps={})
    return interestingness

def selectPuzzles(alignment, min_size, max_size, scoring, threshold, base_filename, mismatches={}, matches={}, gaps={}, difficult=False):
    logging.info('Selecting puzzles')
    windows = interestingWindows(alignment, min_size, max_size, scoring, mismatches, matches, gaps)
    logging.info('Found %d puzzles' % (sum([sum([1 for v in vs if v>threshold]) for i,vs in windows.items()])))
    count = 0
    for n, iscores in windows.items():
        for i, iscore in enumerate(iscores):
            if iscore > threshold:
                block_filename = "%s_offset%d_length%d.fasta" % (base_filename, i, n)
                block = getBlock(alignment, i, n)
                AlignIO.write(block, block_filename, "fasta")
                logging.info('Wrote block %s to %s' % (str(block) , block_filename))
                features_filename = "%s_offset%d_length%d.txt" % (base_filename, i, n)
                with open(features_filename, 'w') as features_file:
                    features_file.write('%s\n' % json.dumps(extractFeatures(block), sort_keys=True))
                if difficulty:
                    with open('%s_difficult_offset%d_length%d.txt'%(base_filename, i, n), 'w') as difficulty_file:
                        difficulty_file.write('%f'%difficulty(block))
                count += 1
    print('Found %d puzzle(s)'%count)
    
def difficulty(seqs):
    x = extractFeatures(seqs)['tree_description']
    m = 2.0327509373773545
    b = -0.59127446582908016
    return m*x + b

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-f', '--format', default='fasta', help='The input file format')
    parser.add_argument('-a', '--alignment', default='clustalw', help='The algorithm to use for multiple sequence alignment')
    parser.add_argument('-t', '--threshold', type=float, default=0.3, help='Threshold for interestingness. Must be between [0,1) with 0 being least interesting')
    parser.add_argument('-v', '--verbose', action='store_true', help='Display in verbose mode')
    parser.add_argument('-d', '--difficulty', action='store_true', help='Compute difficulty of puzzles')
    parser.add_argument('--lsformats', action='store_true', help='Display all available formats')
    parser.add_argument('--lsalgorithms', action='store_true', help='Display all available alignment algorithms')
    parser.add_argument('--features', action='store_true', help='Extract features from sequences')
    parser.add_argument('--min_size', type=int, default=10, help='Smallest window size')
    parser.add_argument('--max_size', type=int, default=20, help='Largest window size')
    
    
    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.INFO)
    if args.lsformats or args.lsalgorithms:
        if args.lsformats:
            print("File formats available : %s"% str(FORMATS))
        if args.lsalgorithms:
            print("Alignments available : %s"% str(ALIGNMENTS))
        exit()
    fileformat, alignment = args.format, args.alignment
    if fileformat not in FORMATS:
        raise NotImplementedError('%s is not a compatible file format. \nUse one of:\n%s'%(fileformat, str(FORMATS)))
    if alignment not in ALIGNMENTS:
        raise NotImplementedError('%s is not an implemented algorithm. \nUse one of:\n%s'%(alignment, str(ALIGNMENTS)))
    if args.features:
        sequences = readSequences(args.infile, fileformat=fileformat)
        sequences = standardizeSequences(sequences)
        msa = align(sequences, alignment)
        extractFeatures(msa)
        args.outfile.write('%s\n' % json.dumps(extractFeatures(sequences), sort_keys=True))
        exit()
    else:
        sequences = readSequences(args.infile, fileformat=fileformat)
        sequences = standardizeSequences(sequences)
        msa = align(sequences, alignment)
        msa = clean(msa)
        mismatches,matches, gaps= {},{},{}
        basename = args.outfile.name.split('.')
        basename = basename[0] if basename[0]!='<stdout>' else 'out'
        n = len(msa[0])
        min_size, max_size = args.min_size, args.max_size
        if min_size > max_size or min_size > n or max_size < 1:
            raise ValueError('Make sure min_size/ max_size are valid')
        selectPuzzles(msa, min_size, max_size, R, args.threshold, basename, mismatches, matches, gaps, args.difficulty)
