#!/usr/bin/env python

from __future__ import print_function
import os,sys
import argparse
import pickle
import subprocess
import re
import pysam 
from collections import defaultdict

from networkx import nx

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

# def is_isoform_match(cs_string):
    # errors = []
    # p = r"[=\+\-\*][A-Za-z]+|~[a-z]+[0-9]+[a-z]+"
    # matches = re.findall(p, cs_string)
    # ed = 0
    # for i, t in enumerate(matches):
    #     # print(t)
    #     e_type = t[0]
    #     # print(t)
    #     if e_type == "=":
    #         length = len(t[1:])
    #     elif e_type == "~":
    #         length = int(t[3:-2])
    #         ed += length
    #         if length >=5:
    #             return False
    #         # print(t, cs_string, matches, matches2)

    #     elif e_type == "*": # substitution
    #         ed += 1

    #     elif e_type == "-": # deletion
    #         length = len(t[1:])
    #         if length >=5:
    #             return False
    #         ed += length

    #     elif e_type == "+": # insertion
    #         length = len(t[1:])
    #         if length >=5:
    #             return False
    #         ed += length

    #     else: # reference skip or soft/hardclip "~", or match =
    #         print("UNEXPECTED!", t)
    #         sys.exit()


    # return True


# def print_out_tsv(nn_sequence_graph, best_exact_matches, reads, references, alignment_file, args):
#     tsv_file = open(os.path.join(args.outfolder, "all_matches.tsv"), "w")
#     tsv_file.write("predicted\treference\tedit_distance\tq_start_offset\tq_end_offset\tref_start_offset\tref_end_offset\tcs\n")
    
#     exact_file = open(os.path.join(args.outfolder, "exact_matches.tsv"), "w")
#     exact_file.write("predicted\treference\tq_len\tref_len\n")

#     isoform_file = open(os.path.join(args.outfolder, "isoform_matches.tsv"), "w")
#     isoform_file.write("predicted\treference\tedit_distance\tq_start_offset\tq_end_offset\tref_start_offset\tref_end_offset\tcs\n")

#     all_exact_matches = set()
#     all_isoform_matches = set()
#     reads_with_isoform_matches = set()

#     for q_acc in nn_sequence_graph:
#         for ref_acc in nn_sequence_graph[q_acc]:
#             read, cigar_ext, edit_distance = nn_sequence_graph[q_acc][ref_acc]
#             # cs_string = read.get_tag("cs")
#             # cigarstring = read.cigarstring

#             ref_start_offset = read.reference_start
#             ref_end_offset =  len(references[ref_acc]) - read.reference_end
#             q_start_offset =  read.query_alignment_start
#             q_end_offset =  len(reads[q_acc]) - read.query_alignment_end
#             if ref_start_offset == ref_end_offset == q_start_offset == q_end_offset == edit_distance == 0:
#                 # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset))
#                 exact_file.write("{0}\t{1}\t{2}\t{3}\n".format(q_acc, ref_acc, len(reads[q_acc]), len(references[ref_acc])))
#                 all_exact_matches.add(ref_acc)

#             if ref_start_offset <= 5 and ref_end_offset <= 5 and q_start_offset <= 5 and q_end_offset <= 5:
#                 cs_string = read.get_tag("cs")
#                 if is_isoform_match(cs_string):             
#                     # print("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))
#                     isoform_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))
#                     all_isoform_matches.add(ref_acc)
#                     reads_with_isoform_matches.add(q_acc)


#             # print(cs_string)
#             tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(q_acc, ref_acc, edit_distance, q_start_offset, q_end_offset, ref_start_offset, ref_end_offset, cigar_ext))

#     tsv_file.close()
#     exact_file.close()
#     isoform_file.close()
#     print("Unique exact matches:", len(all_exact_matches))
#     print("Unique isoform matches:", len(all_isoform_matches))

#     predictions_without_match = set(reads.keys()) - reads_with_isoform_matches
#     predictions_without_match_file = open(os.path.join(args.outfolder, "predictions_without_match.tsv"), "w")
#     for q_acc in predictions_without_match:
#         predictions_without_match_file.write("{0}\t{1}\n".format(q_acc, len(reads[q_acc])))
#     predictions_without_match_file.close()
#     print("Predictions without isoform match:", len(predictions_without_match))
#     print("total nr of predictions:", len(reads))


#     all_matches = set([ref_acc for q_acc in nn_sequence_graph for ref_acc in nn_sequence_graph[q_acc]])
#     no_hits = set(references.keys()) - all_matches
#     print("No hits:", len(no_hits), no_hits)
#     no_exact_hits = set(references.keys()) - all_exact_matches
#     for q_acc in nn_sequence_graph:
#         for ref_acc in no_exact_hits:
#             if ref_acc in nn_sequence_graph[q_acc]:
#                 read = nn_sequence_graph[q_acc][ref_acc][0]
#                 # print(nn_sequence_graph[q_acc][ref_acc][1], read.reference_start, len(references[ref_acc]) - read.reference_end, ref_acc)
#     print("No exact hits:", len(no_exact_hits), no_exact_hits)

#     return tsv_file.name


def get_exon_starts_and_stop(matches, first_exon_start):
    exon_coords = []
    ref_pos = first_exon_start
    
    current_exon_start = ref_pos
    for i, t in enumerate(matches):
        # print(t)
        e_type = t[0]
        # print(t)
        if e_type == "=":
            length = len(t[1:])
            ref_pos += length
        elif e_type == "~":
            current_exon_stop = ref_pos
            exon_coords.append( (current_exon_start, current_exon_stop) )

            length = int(t[3:-2])
            ref_pos += length
            current_exon_start = ref_pos
            # print(t)

        elif e_type == "*": # substitution
            ref_pos += 1

        elif e_type == "-": # deletion
            length = len(t[1:])
            ref_pos += length

        elif e_type == "+": # insertion
            length = len(t[1:])
            ref_pos += 0

        else: # reference skip or soft/hardclip "~", or match =
            print("UNEXPECTED!", t)
            sys.exit()

    #last exon
    current_exon_stop = ref_pos
    exon_coords.append( (current_exon_start, current_exon_stop) )

    return exon_coords

def is_same_isoform(q_isoform, ref_isoform):
    # compare cs tag at intron sites
    q_cs_string = q_isoform.get_tag("cs")
    q_start = q_isoform.reference_start
    q_end = q_isoform.reference_end

    ref_cs_string = ref_isoform.get_tag("cs")
    ref_start = ref_isoform.reference_start
    ref_end = ref_isoform.reference_end
    
    errors = []
    p = r"[=\+\-\*][A-Za-z]+|~[a-z]+[0-9]+[a-z]+"
    
    q_matches = re.findall(p, q_cs_string)
    ref_matches = re.findall(p, ref_cs_string)    
    # print(q_start, q_end, ref_start, ref_end)

    q_exons = get_exon_starts_and_stop(q_matches, q_start)
    all_q_exons = set(q_exons)
    ref_exons = get_exon_starts_and_stop(ref_matches, ref_start)
    # print(len(all_q_exons), len(ref_exons))
    if len(all_q_exons) != len(ref_exons):
        return False
    for r_start, r_stop in ref_exons:
        if (r_start, r_stop) not in all_q_exons:
            # if "transcript_846_" in q_isoform.query_name:            
            #     print(q_isoform.query_name, ref_isoform.query_name, r_start, r_stop, sorted(all_q_exons))
            #     print("start and end!", q_start, q_end, ref_start, ref_end)
            # print(False)
            return False

    # print(sorted(all_q_exons))
    # print(ref_isoform.query_name, sorted(ref_exons))
    # print()
    # print(True)
    return True




def detect_isoforms(ref_samfile_path, pred_samfile_path, outfile):

    ref_samfile = pysam.AlignmentFile(ref_samfile_path, "r", check_sq=False)
    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)

    ref_isoforms = [ref_isoform for ref_isoform in ref_samfile.fetch(until_eof=True)] 
    query_isoforms = [q_isoform for q_isoform in pred_samfile.fetch(until_eof=True)]
    counter_old = 0
    counter_new = 0
    ref_to_queries = { ref.query_name : set() for ref in ref_isoforms }
    queries_to_ref = { query.query_name : set() for query in query_isoforms }
    new_isoforms = set()
    for q_isoform in query_isoforms:
        is_new = True
        for ref_isoform in ref_isoforms:
            if is_same_isoform(q_isoform, ref_isoform) and is_same_isoform(ref_isoform, q_isoform):
                # print("YO")
                queries_to_ref[q_isoform.query_name].add(ref_isoform.query_name)
                if len(queries_to_ref[q_isoform.query_name]) > 1:
                    print("More than 1 ref")
                    print("Same", q_isoform.query_name, queries_to_ref[q_isoform.query_name] )
                is_new = False
                counter_old += 1
            
        if is_new:
            counter_new += 1
            print("New", q_isoform.query_name)
            new_isoforms.add(q_isoform.query_name)

        else:
            assert len(queries_to_ref[q_isoform.query_name]) == 1
            ref = queries_to_ref[q_isoform.query_name].pop()
            ref_to_queries[ref].add(q_isoform.query_name)

    print([ len(ref_to_queries[r]) for r in ref_to_queries] )
    total_predictions = len(query_isoforms)
    print(total_predictions, "Total predictions")
    print(counter_old, "predictions had the same isoform structure as ref")
    print(counter_new, "predictions had new isoform structure to ref")

    return queries_to_ref, new_isoforms


def group_novel_isoforms(new_isoforms, pred_samfile_path, outfile):
    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)
    query_isoforms = [q_isoform for q_isoform in pred_samfile.fetch(until_eof=True) if q_isoform.query_name in new_isoforms]
    G = nx.Graph()
    for n in new_isoforms:
        G.add_node(n)

    for i1 in query_isoforms:
        for i2 in query_isoforms:
            if i1.query_name == i2.query_name:
                continue
            else:
                if is_same_isoform(i1, i2) and is_same_isoform(i2, i1):
                    G.add_edge(i1.query_name, i2.query_name)

    maximal_cliques = [cl for cl in nx.find_cliques(G)]

    print([ len(cl) for cl in maximal_cliques] )
    print(len([ len(cl) for cl in maximal_cliques]), "unique splice sites isoforms")

    return maximal_cliques


def get_junction_novelty(q_isoform, ref_isoform):
    # compare cs tag at intron sites
    q_cs_string = q_isoform.get_tag("cs")
    q_start = q_isoform.reference_start
    q_end = q_isoform.reference_end

    ref_cs_string = ref_isoform.get_tag("cs")
    ref_start = ref_isoform.reference_start
    ref_end = ref_isoform.reference_end
    
    errors = []
    p = r"[=\+\-\*][A-Za-z]+|~[a-z]+[0-9]+[a-z]+"
    
    q_matches = re.findall(p, q_cs_string)
    ref_matches = re.findall(p, ref_cs_string)    
    # print(q_start, q_end, ref_start, ref_end)

    q_exons = get_exon_starts_and_stop(q_matches, q_start)
    ref_exons = get_exon_starts_and_stop(ref_matches, ref_start)
    all_ref_exons = set(ref_exons)

    # print(len(all_q_exons), len(ref_exons))
    if len(all_ref_exons) != len(q_exons):
        return None #("exon combination", len(all_q_exons))
    junction_novelties = []
    for q_start, q_stop in q_exons:
        if (q_start, q_stop) not in all_ref_exons:
            junction_novelties.append((q_start, q_stop))
            # if "transcript_846_" in q_isoform.query_name:            
            #     print(q_isoform.query_name, ref_isoform.query_name, r_start, r_stop, sorted(all_q_exons))
            #     print("start and end!", q_start, q_end, ref_start, ref_end)
            # print(False)
            return tuple(sorted(junction_novelties))

    else:
        print("BUG")
        sys.exit()


def get_novelty_feature(new_isoforms, pred_samfile_path, ref_samfile_path, outfile):
    pred_samfile = pysam.AlignmentFile(pred_samfile_path, "r", check_sq=False)
    ref_samfile = pysam.AlignmentFile(ref_samfile_path, "r", check_sq=False)

    novel_query_isoforms = [q_isoform for q_isoform in pred_samfile.fetch(until_eof=True) if q_isoform.query_name in new_isoforms]
    ref_isoforms = [ref_isoform for ref_isoform in ref_samfile.fetch(until_eof=True)] 
    features = {}
    for i1 in novel_query_isoforms:
        features[i1.query_name] = []
        for i2 in ref_isoforms:
            junction_novelty_to_ref = get_junction_novelty(i1, i2)
            if junction_novelty_to_ref:
                features[i1.query_name].append(junction_novelty_to_ref)

    for q_acc in features.keys():
        if not features[q_acc]:
            print()
            print("EXON NOVELTY", q_acc)
            print()
            features[q_acc] = ["exon_combination"]
        else:
            pass
            # print(q_acc, features[q_acc])

    unique_splices = {}
    for q_acc in features:
        # print(q_acc)
        assert type(q_acc) is str 
        # print(features[q_acc])
        # print( type(set(features[q_acc])) )
        # print( set(features[q_acc]) )
        # assert set(features[q_acc]) is set
        unique_splices[str(q_acc)] = set([ tup for item in features[q_acc] for tup in item ] )
    # unique_splices = { q_acc : set(features[q_acc]) for q_acc in features.keys() }
    exon_ranges = { (147014219,147014282) : 9,
                    (147018023,147018132) : 10,
                    (147018985,147019119) : 11,
                    (147019618,147019680) : 12,
                    (147022095,147022181) : 13,
                    (147024651,147024846) : 14,
                    (147026389,147026571) : 15,
                    (147027054,147027136) : 16,
                    (147030203,147030451) : 17 }

    new_tags = {}
    for q_acc in unique_splices.keys():
        tag = ""
        if unique_splices[q_acc] == set("exon_combination"):
            tag = "exon_combination"
            new_tags[q_acc] = tag
            continue

        for splice in unique_splices[q_acc]:
            # print(unique_splices[q_acc])


            for start, stop in exon_ranges:
                # print( splice)
                if splice[0] <= start and stop <= splice[1]:
                    tag += 'spans' + str(exon_ranges[(start,stop)]) + ";"
                elif start <= splice[0] and splice[1] <= stop:
                    tag += 'within' + str(exon_ranges[(start,stop)]) + ";"
                elif splice[0] <= start and splice[1] > start:
                    tag += '5"overlap' + str(exon_ranges[(start,stop)]) + ";"
                elif splice[0] >= start and splice[1] > stop:
                    tag += '3"overlap' + str(exon_ranges[(start,stop)]) + ";"

                # else:
                #     print("BUG")
                #     sys.exit()
        
        if not tag:
            print(q_acc, unique_splices[q_acc])
        assert tag != ""
        new_tags[q_acc] = tag
        # print(q_acc, unique_splices[q_acc])

    for q_acc in new_tags:
        print(new_tags[q_acc], q_acc)
    return new_tags


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def main(args):
    
    # filtered_predictions = {acc : seq for acc, seq in read_fasta(open(args.predictions,"r"))} 

    outfile = open(os.path.join(args.outfolder, args.prefix + ".fa"), "w")
    queries_to_ref, new_isoforms = detect_isoforms(args.refsamfile, args.querysamfile, outfile)
    new_clusters = group_novel_isoforms(new_isoforms, args.querysamfile, outfile)
    new_isoform_tags = get_novelty_feature(new_isoforms, args.querysamfile, args.refsamfile, outfile)

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

    