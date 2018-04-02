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
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

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
            # ed += length
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


def get_introns():
    coords = []
    return 

def detect_isoforms(ref_samfile_path, pred_samfile_path, outfile):

    ref_samfile = pysam.AlignmentFile(ref_samfile_path, "r", check_sq=False)
    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)

    introns = ref_samfile.find_introns(ref_samfile.fetch(until_eof=True))
    print(introns)
    print(len(introns))
    sys.exit()
    counter = 0
    references = SAM_file.references
    for read in ref_samfile.fetch(until_eof=True):
        if read.is_unmapped:
            print(read.query_name, "is unmapped!")
            continue

        if read.is_secondary:
            print(read.query_name, "is secondary!")
            continue

        ref_name = references[read.reference_id]

        cigar_tuples = read.cigartuples # list of tuples of (operation, length)
        # print(cigar_tuples)
        # has to score an "error-free" alignment containing an exon diff higher than the same exon structure with lots of indels/mismatches!
        # What could possibly work is to count only smaller diffs towards total ed, for example all op with length_ < 5 ? and low q values?
        # we need the quality stringe here in that case...
        # also we should keep the quality strings around when we correct the strings! reads should always be paired with the quality values
        # if fasta is sent as input -- make all quality values the same and use coverage
        # ref_start_offset = read.reference_start
        # ref_end_offset =  len(acc_to_seq[ref_name]) - read.reference_end
        q_start_offset =  read.query_alignment_start
        q_end_offset =  len(reads[read.query_name]) - read.query_alignment_end

        cs_string = read.get_tag("cs")
        ext_cigarstring, tot_ed = cs_to_cigar_and_ed(cs_string)
        read_length = float(len(reads[read.query_name]))
        alignment_id = (read_length - float(tot_ed)) / read_length # ( read.infer_read_length - float(tot_ed) ) / float(read.infer_read_length)
        alignment_coverage =  len(read.query_alignment_sequence)  / read_length # len(read.query_alignment_sequence)  / float(read.infer_read_length)

        filter_read = False
        if alignment_id < ALIGN_IDENTITY:
            print(read.query_name, "Below identity", alignment_id)  
            filter_read = True

        if alignment_coverage < ALIGN_COVERAGE:
            print(read.query_name, "Below coverage", alignment_coverage)  
            filter_read = True

        if not (ALIGNMENT_START[0] <= read.reference_start <= ALIGNMENT_START[1]):
            print(read.query_name, "Bad start", read.reference_start)  
            filter_read = True

        if not (ALIGNMENT_END[0] <= read.reference_end <= ALIGNMENT_END[1]):
            print(read.query_name, "Bad end", read.reference_start) 
            filter_read = True

        if not filter_read:
            counter += 1
            outfile.write(">{0}\n{1}\n".format(read.query_name, reads[read.query_name] ))

    print(counter, "remaining after filtering")


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def main(args):
    
    # filtered_predictions = {acc : seq for acc, seq in read_fasta(open(args.predictions,"r"))} 

    outfile = open(os.path.join(args.outfolder, args.prefix + ".fa"), "w")
    detect_isoforms(args.refsamfile, args.querysamfile, outfile)
    # nn_sequence_graph = best_matches_to_accession_graph(all_matches, acc_to_seq)
    # print_out_tsv(nn_sequence_graph, all_matches,  reads, references, alignment_file, args)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('refsamfile', type=str, help='Samfile.')
    parser.add_argument('querysamfile', type=str, help='Samfile.')
    # parser.add_argument('predictions', type=str, help='Fasta file with only filtered isoform hits to FMR region (output of "filter_hits_on_hg19" script).')
    parser.add_argument('outfolder', type=str, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')  
    parser.add_argument('prefix', type=str, help='prefix to outfile.')  

    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)

    