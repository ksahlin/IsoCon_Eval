#!/usr/bin/env python

from __future__ import print_function
import os,sys
import argparse
import pickle
import subprocess
import re
import pysam 
from collections import defaultdict



try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd


'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {"-" : "-", 'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)


def read_fasta(fasta_file):
    fasta_seqs = {}
    k = 0
    temp = ''
    accession = ''
    for line in fasta_file:
        if line[0] == '>' and k == 0:
            accession = line[1:].strip().split()[0]
            fasta_seqs[accession] = ''
            k += 1
        elif line[0] == '>':
            yield accession, temp
            temp = ''
            accession = line[1:].strip().split()[0]
        else:
            temp += line.strip()
    yield accession, temp


def cs_to_cigar_and_ed(cs_string):
    errors = []
    # p = r"[=\+\-\~\*][A-Za-z]+"
    p = r"[=\+\-\*][A-Za-z]+|~[a-z]+[0-9]+[a-z]+"
    matches = re.findall(p, cs_string)
    # matches2 = re.findall(p2, cs_string)
    #occurences_by_type = {}
    # print("NEW")
    cigar_ext= ""
    ed = 0
    for i, t in enumerate(matches):
        # print(t)
        e_type = t[0]
        length = len(t[1:])
        # print(t)
        if e_type == "=":
            cigar_ext += "{0}{1}".format(length, "=")
        elif e_type == "~":
            length = int(t[3:-2])
            cigar_ext += "{0}{1}".format( length, "N")
            ed += length
            # print(t, cs_string, matches, matches2)

        # here we store smaller errors/variations
        elif e_type == "*": # substitution
            cigar_ext += "{0}".format("*")
            ed += 1

        elif e_type == "-": # deletion
            cigar_ext += "{0}{1}".format(length, "D")
            ed += length

        elif e_type == "+": # insertion
            cigar_ext += "{0}{1}".format(length, "I")
            ed += length

        else: # reference skip or soft/hardclip "~", or match =
            print(t)
    # print(cigar_ext, ed)
    return cigar_ext, ed

def is_isoform_match(cs_string):
    errors = []
    p = r"[=\+\-\*][A-Za-z]+|~[a-z]+[0-9]+[a-z]+"
    matches = re.findall(p, cs_string)
    ed = 0
    for i, t in enumerate(matches):
        # print(t)
        e_type = t[0]
        # print(t)
        if e_type == "=":
            length = len(t[1:])
        elif e_type == "~":
            length = int(t[3:-2])
            ed += length
            if length >=5:
                return False
            # print(t, cs_string, matches, matches2)

        elif e_type == "*": # substitution
            ed += 1

        elif e_type == "-": # deletion
            length = len(t[1:])
            if length >=5:
                return False
            ed += length

        elif e_type == "+": # insertion
            length = len(t[1:])
            if length >=5:
                return False
            ed += length

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()


    return True


def print_out_tsv(nn_sequence_graph, best_exact_matches, reads, references, alignment_file, args):
    tsv_file = open(os.path.join(args.outfolder, "all_matches.tsv"), "w")
    tsv_file.write("predicted\treference\tedit_distance\tq_start_offset\tq_end_offset\tref_start_offset\tref_end_offset\tcs\n")
    
    exact_file = open(os.path.join(args.outfolder, "exact_matches.tsv"), "w")
    exact_file.write("predicted\treference\tq_len\tref_len\n")

    isoform_file = open(os.path.join(args.outfolder, "isoform_matches.tsv"), "w")
    isoform_file.write("predicted\treference\tedit_distance\tq_start_offset\tq_end_offset\tref_start_offset\tref_end_offset\tcs\n")

    all_exact_matches = set()
    all_isoform_matches = set()
    reads_with_isoform_matches = set()

    for q_acc in nn_sequence_graph:
        for ref_acc in nn_sequence_graph[q_acc]:
            read, cigar_ext, edit_distance = nn_sequence_graph[q_acc][ref_acc]
            # cs_string = read.get_tag("cs")
            # cigarstring = read.cigarstring

            ref_start_offset = read.reference_start
            ref_end_offset =  len(references[ref_acc]) - read.reference_end
            q_start_offset =  read.query_alignment_start
            q_end_offset =  len(reads[q_acc]) - read.query_alignment_end
            if ref_start_offset == ref_end_offset == q_start_offset == q_end_offset == edit_distance == 0:
                # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset))
                exact_file.write("{0}\t{1}\t{2}\t{3}\n".format(q_acc, ref_acc, len(reads[q_acc]), len(references[ref_acc])))
                all_exact_matches.add(ref_acc)

            if ref_start_offset <= 5 and ref_end_offset <= 5 and q_start_offset <= 5 and q_end_offset <= 5:
                cs_string = read.get_tag("cs")
                if is_isoform_match(cs_string):             
                    # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))
                    isoform_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))
                    all_isoform_matches.add(ref_acc)
                    reads_with_isoform_matches.add(q_acc)


            # print(cs_string)
            tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))

    tsv_file.close()
    exact_file.close()
    isoform_file.close()
    print("Unique exact matches:", len(all_exact_matches))
    print("Unique isoform matches:", len(all_isoform_matches))

    predictions_without_match = set(reads.keys()) - reads_with_isoform_matches
    predictions_without_match_file = open(os.path.join(args.outfolder, "predictions_without_match.tsv"), "w")
    for q_acc in predictions_without_match:
        predictions_without_match_file.write("{0}\t{1}\n".format(q_acc, len(reads[q_acc])))
    predictions_without_match_file.close()
    print("Predictions without isoform match:", len(predictions_without_match))
    print("total nr of predictions:", len(reads))


    all_matches = set([ref_acc for q_acc in nn_sequence_graph for ref_acc in nn_sequence_graph[q_acc]])
    no_hits = set(references.keys()) - all_matches
    print("No hits:", len(no_hits), no_hits)
    no_exact_hits = set(references.keys()) - all_exact_matches
    for q_acc in nn_sequence_graph:
        for ref_acc in no_exact_hits:
            if ref_acc in nn_sequence_graph[q_acc]:
                read = nn_sequence_graph[q_acc][ref_acc][0]
                # print(nn_sequence_graph[q_acc][ref_acc][1], read.reference_start, len(references[ref_acc]) - read.reference_end, ref_acc)
    print("No exact hits:", len(no_exact_hits), no_exact_hits)

    return tsv_file.name


def SAM_to_all_matches(SAM_file_path, acc_to_seq, ignore_ends_length, read_qualities = {}):
    all_matches = {}
    current_min_ed = {}

    SAM_file = pysam.AlignmentFile(SAM_file_path, "r", check_sq=False)
    references = SAM_file.references
    for read in SAM_file.fetch(until_eof=True):
        if read.query_name not in all_matches:
            all_matches[read.query_name] = []
        if read.is_unmapped:
            print(read.query_name, "is unmapped!")
            continue

        if read.is_reverse:
            # print("OMG")
            continue
        ref_name = references[read.reference_id]

        if read.query_name == ref_name:
            # print("Self alignment detected, skipping")
            continue
        cigar_tuples = read.cigartuples # list of tuples of (operation, length)
        # print(cigar_tuples)
        # has to score an "error-free" alignment containing an exon diff higher than the same exon structure with lots of indels/mismatches!
        # What could possibly work is to count only smaller diffs towards total ed, for example all op with length_ < 5 ? and low q values?
        # we need the quality stringe here in that case...
        # also we should keep the quality strings around when we correct the strings! reads should always be paired with the quality values
        # if fasta is sent as input -- make all quality values the same and use coverage
        ref_start_offset = read.reference_start
        ref_end_offset =  len(acc_to_seq[ref_name]) - read.reference_end
        q_start_offset =  read.query_alignment_start
        q_end_offset =  len(acc_to_seq[read.query_name]) - read.query_alignment_end

        # tot_ed = 0
        # if ref_start_offset > 0 and ref_start_offset > ignore_ends_length:
        #    tot_ed += ref_start_offset

        # if ref_end_offset > 0 and ref_end_offset > ignore_ends_length:
        #    tot_ed += ref_end_offset


        # if cigar_tuples[0][0] == 4 and cigar_tuples[0][1] <= ignore_ends_length: # ignore first region when query longer than ref 
        #     cigar_tuples.pop(0)

        # if cigar_tuples[-1][0] == 4 and cigar_tuples[-1][1] <= ignore_ends_length: # ignore first region when length difference
        #     cigar_tuples.pop(-1)
        
        # here we will compare "absolute nearest neighbors" with "probabilistic nearest neighbors"
        cs_string = read.get_tag("cs")
        ext_cigarstring, tot_ed = cs_to_cigar_and_ed(cs_string)
        

        all_matches[read.query_name].append( (ref_name, read, ext_cigarstring, tot_ed) )

    return all_matches, SAM_file


def best_matches_to_accession_graph(best_matches, acc_to_seq):
    nr_edges = 0
    tot_aln_distance = 0
    no_nn = 0
    more_than_one_nn = 0
    nearest_neighbor_graph = {}
    for s_acc in best_matches:
        nearest_neighbor_graph[s_acc] = {}
        if len(best_matches[s_acc]) == 0:
            no_nn +=1
        if len(best_matches[s_acc]) >1 :
            more_than_one_nn += 1

        for nn_acc, s_aln_obj, cigar_ext, tot_ed in best_matches[s_acc]:
            nearest_neighbor_graph[s_acc][nn_acc] = (s_aln_obj, cigar_ext, tot_ed)
            nr_edges += 1
            # tot_aln_distance += current_min_ed[s_acc]


    # print("TOTAL SW ALIGNMENT DISTANCE:", tot_aln_distance)
    print("TOTAL NR EDGES:", nr_edges)
    # if nr_edges > 0:
    #     print("Average distance per edge:", tot_aln_distance/ float(nr_edges))
    print("Number of seqs without NN:", no_nn)
    print("Number of seqs with more than one NN:", more_than_one_nn)

    # could possibly store alignments here as well by implementing cigar to alignment or -cs to alignment
    return nearest_neighbor_graph


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
    user_matrix = parasail.matrix_create("ACGT", 2, -4)
    result = parasail.nw_trace_scan_32(query, ref, 12, 1, user_matrix)
    if result.saturated:
        print("SATURATED!")
    else:
        if sys.version_info[0] < 3:
            cigar_string = str(result.cigar.decode).decode('utf-8')
        else:
            cigar_string = str(result.cigar.decode, 'utf-8')
        
        return cigar_string, cigar_to_seq(cigar_string, query, ref)


def align(targets, queries, outfolder):
    print('Aligning with parasail.')
    for q_acc, q_seq in queries.items():
        for t_acc, t_seq in targets.items():
            # q_rc, t_rc = reverse_complement(q_seq), reverse_complement(t_seq)
            cigar_string, (q_aln, r_aln) = parasail_alignment(q_seq, t_seq)
            # q_aln, r_aln = reverse_complement(q_aln), reverse_complement(r_aln)
            if q_acc == "PB.1.35|groupE" or q_acc == "PB.1.22|groupE":
                print(cigar_string)
                print(q_aln)
                print(q_acc, q_seq)

    # SAM_file = queries + ".sam"


    return #SAM_file

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def main(args):

    # map reads to references with parasail N-W
    references =  {acc : seq.upper() for acc, seq in read_fasta(open(args.references,"r"))} 
    predictions =  {acc : seq.upper() for acc, seq in read_fasta(open(args.predictions,"r"))} 
    SAM_file = align(references, predictions, args.outfolder) # unique strings are targets, not_in clusters are queries
    

    acc_to_seq  = merge_two_dicts(reads, references)
    all_matches, alignment_file = SAM_to_all_matches(SAM_file, acc_to_seq, 100, read_qualities = {})
    nn_sequence_graph = best_matches_to_accession_graph(all_matches, acc_to_seq)

    print_out_tsv(nn_sequence_graph, all_matches,  reads, references, alignment_file, args)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('predictions', type=str, help='CCS reads to map in fasta format.')
    parser.add_argument('references', type=str, help='Biological material to map to.')
    parser.add_argument('outfolder', type=str, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    # parser.add_argument('--realign', action="store_true", help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)

    