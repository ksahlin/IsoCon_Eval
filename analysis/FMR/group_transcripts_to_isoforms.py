
from __future__ import print_function
import os,sys
import argparse
import re
import ssw
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
from matplotlib_venn import venn3, venn3_circles, venn2
import edlib
import math
import errno
import networkx as nx

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


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def cigar_to_seq(cigar, query, ref):
    cigar_tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        cigar_tuples.append((int(length), type_ ))

    r_index = 0
    q_index = 0
    q_aln = []
    r_aln = []
    for length_ , type_ in cigar_tuples:
        if type_ == "=" or type_ == "X":
            q_aln.append(query[q_index : q_index + length_])
            r_aln.append(ref[r_index : r_index + length_])

            r_index += length_
            q_index += length_
        
        elif  type_ == "I":
            # insertion w.r.t. reference
            r_aln.append('-' * length_)
            q_aln.append(query[q_index: q_index + length_])
            #  only query index change
            q_index += length_

        elif type_ == 'D':
            # deletion w.r.t. reference
            r_aln.append(ref[r_index: r_index + length_])
            q_aln.append('-' * length_)
            #  only ref index change
            r_index += length_
        
        else:
            print("error")
            print(cigar)
            sys.exit()

    return  "".join([s for s in q_aln]), "".join([s for s in r_aln])


from time import time
import parasail
def parasail_alignment(query, ref):
    user_matrix = parasail.matrix_create("ACGT", 2, -20)
    result = parasail.nw_trace_scan_16(query, ref, 50, 0, user_matrix)
    if result.saturated:
        print("SATURATED!")
    else:
        return cigar_to_seq(result.cigar.decode, query, ref)


def create_grouping_graph(transcripts):
    
    G = nx.Graph()
    for acc in transcripts.keys():
        G.add_node(acc, acc_id = int(acc.split("_")[1]), support = acc.split("_")[4])

    cntr = 0

    processed = set()
    # already_assigned = set()

    for acc1, seq1 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
        cntr += 1
        processed.add(acc1)
        # print("length t:", len(seq1), acc1)
        if cntr % 5 == 0:
            print(cntr, "sequences processed")

        for acc2, seq2 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
            if acc2 in processed:
                continue
            
            start = time()
            seq1_aln, seq2_aln = parasail_alignment(seq1, seq2)
            assert len(seq1_aln) == len(seq2_aln)
            total_elapsed = time() - start
            # print("time:", total_elapsed )

            del_seq1 = re.findall(r"[-]+",seq1_aln)
            del_seq2 = re.findall(r"[-]+",seq2_aln)
            mismatches = len([ 1 for n1, n2 in zip(seq1_aln,seq2_aln) if n1 != n2 and n1 != "-" and n2 != "-" ])

            # ## do not count length discrepancies in ends
            # inner_del_seq1 = re.findall(r"[AGCT][-]+[AGCT]",seq1_aln)
            # inner_del_seq2 = re.findall(r"[AGCT][-]+[AGCT]",seq2_aln)
            # total_inner = sum([len(d) - 2 for d in inner_del_seq1]) + sum([len(d) - 2 for d in inner_del_seq2])

            # by default (since all transcripts are distinct if we end up here), each transcript is its on gene member
            # if we find an alingment that contains only structural changes of > X (2) nucleotides, and no other smaller differences we classify as same family
            del_lengths1 = [len(del_) for del_ in del_seq1]
            del_lengths2 = [len(del_) for del_ in del_seq2]

            no_exon_del_in_seq1 = ((len(del_lengths1) > 0 and max(del_lengths1) <= 5) or len(del_lengths1)  == 0)
            no_exon_del_in_seq2 = ((len(del_lengths2) > 0 and min(del_lengths2) <= 5) or len(del_lengths2)  == 0)
            if no_exon_del_in_seq1 and no_exon_del_in_seq2:
                G.add_edge(acc1, acc2, 
                            alignment={ acc1 : seq1_aln,  acc2 : seq2_aln })



    list_of_maximal_cliques = list(nx.find_cliques(G))
    print("Number of possible isoforms:", len(list_of_maximal_cliques) )
    print("member clique sizes", [ len(cl) for cl in  sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True)] )

    return G



def main(args):

    transcripts = { acc : seq.upper()  for (acc, seq) in  read_fasta( open(args.transcripts, "r") ) }
    print("number of transcripts:", len(transcripts))
    G = create_grouping_graph(transcripts)
    list_of_maximal_cliques = list(nx.find_cliques(G))
    clique_id_to_accessions = {member_id : cl for member_id, cl in enumerate(sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True))}



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    # parser.add_argument('--inexact',  action='store_true', help='Check inexact matches')

    parser.add_argument('--transcripts', type=str, required=True, help='files')
    parser.add_argument('--outfolder', type=str, help='Output path of results')

    args = parser.parse_args()
    mkdir_p(args.outfolder)
    main(args)

