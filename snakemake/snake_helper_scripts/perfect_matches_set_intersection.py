
from __future__ import print_function
import os,sys
import argparse
import re
# import ssw
import numpy as np
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import pandas as pd
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from matplotlib_venn import venn3, venn3_circles
import edlib
import math
import errno

def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    if ed == -1:
        ed= len(x)

    if task == "path":
        locations =  result["locations"]
        cigar =  result["cigar"]
        return (ed, locations, cigar)
    else:
        return (ed, None, None)

def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip()
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip()
        else:
            temp += line.strip().upper()
    if accession:
        yield accession, temp

# def check_inexact(a,b,c):
#     for seq1 in a:
#         best_ed = len(seq1)
#         best_cig = ""
#         best_loc = ""
#         for seq2 in c:
#             if math.fabs(len(seq1) - len(seq2)) > best_ed:
#                 continue

#             edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="NW", task="path", k=min(100, best_ed))
#             if edit_distance < best_ed:
#                 best_ed = edit_distance
#                 best_cig = cigar
#                 best_loc = locations

#         if best_ed > 0:
#             print(best_ed,best_cig, best_loc)
#             if best_ed > 500:
#                 print(seq1)

from collections import defaultdict

def main(args):

    hits = defaultdict(list) # {"ICE": [], "ISOCON" : [], "FLNC": []}

    for i, line in enumerate(open(args.tsv_file, "r")):
        if i ==0: # ignore header
            continue
        # ID      METHOD  GENE_FAMILY     ED
        transcript_id, method, family, ed, clipped = line.split("\t")
        hits[method].append(transcript_id)


    a = set(hits["FLNC"])
    b = set(hits["ICE"])
    c = set(hits["ISOCON"])

    # if args.inexact:
    #     check_inexact(a,b,c)

    print(len(a), len(b), len(c))

    a_not_b_c = a - (b | c)
    b_not_a_c = b - (a | c)
    a_b_not_c = (a & b) - c
    c_not_a_b = c - (a | b)
    a_c_not_b = (a & c) - b
    b_c_not_a = (b & c) - a
    a_b_c = a & b & c


    r = venn3([a, b, c], ("flnc", "ICE", "ISOCON"))
    plt.savefig(args.outfile)


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('tsv_file', type=str, help='Path to the consensus fasta file')
    # parser.add_argument('--inexact',  action='store_true', help='Check inexact matches')

    # parser.add_argument('--names', type=str, nargs=3, required=True, help='Set names')

    parser.add_argument('outfile', type=str, help='Output path of results')

    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)

    main(args)