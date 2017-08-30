
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
from matplotlib_venn import venn2, venn2_circles
import edlib
import math

import dataset
from collections import defaultdict
import errno

class Prediction(object):
    """docstring for Prediction
        predicted_acc  sample1 sample2  both  acc_sample1  acc_sample2  
        
        Illumina_support  full_supp_sample1 full_supp_sample2 gene_family 
        
        perfect_match_database  category  gene_member_number

           dict[(q_name, primed_id)] = predicted_object

        """
    def __init__(self, predicted_acc, sequence, primer, batch):
        super(Prediction, self).__init__()
        self.predicted_acc = predicted_acc
        self.sequence = sequence
        self.primer = primer
        self.batch = batch

        self.sample1 = ""
        self.sample2 = ""
        self.acc_sample1 = ""
        self.acc_sample2 = ""

        self.Illumina_support = None
        self.full_supp_sample1 = None
        self.gene_family = ""

        self.perfect_match_database = ""
        self.category = "" # coding/pseudo/error
        self.gene_member_number = None #any integer/string grouping it together with other isoforms from the same number, e.g. GM1, PM1, etc..


primer_to_family_and_sample = { "barcode_1-2kb" : { 0 : ("sample1", "HSFY2"),
                                                    1 :  ("sample2", "HSFY2"),
                                                    2 : ("sample1", "RBMY"),
                                                    3 : ("sample1", "CDY1"),
                                                    4 : ("sample1", "CDY2"),
                                                    5 : ("sample1", "DAZ"),
                                                    6 : ("sample2", "RBMY"),
                                                    7 : ("sample2", "CDY1"),
                                                    8 : ("sample2", "CDY2"),
                                                    9 : ("sample2", "DAZ")
                                                },
                                "barcode_1kb" : { 0 : ("sample1", "BPY"),
                                                  1 : ("sample1", "VCY"),
                                                  2 : ("sample1", "XKRY"),
                                                  3 : ("sample1", "PRY"),
                                                  4 : ("sample1", "HSFY1"),
                                                  5 : ("sample1", "TSPY"),
                                                  6 : ("sample2", "BPY"),
                                                  7 : ("sample2", "VCY"),
                                                  8 : ("sample2", "XKRY"),
                                                  9 : ("sample2", "PRY"),
                                                  10 : ("sample2", "HSFY1"),
                                                  11 : ("sample2", "TSPY")
                                                    }    }


def edlib_ed(x, y, mode="NW", task="distance", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    if task == "path":
        cigar =  result["cigar"]
        locations =  result["locations"]
        return ed, cigar, locations
    else:
        return ed

def get_minimizers_2set_simple(querys, targets, min_query_len):
    best_edit_distances = {}
    best_cigars = {}
    i = 1
    for acc1, seq1 in querys.items():
        if i % 200 == 0:
            print("processing candidate", i)
        best_ed = 2*len(seq1)
        best_edit_distances[acc1] = {}
        best_cigars[acc1] = {}
        # if acc1 == "transcript_447_support_14_15_1.57320042782e-125_21_5":
        #     print("NOW!")
        if len(seq1) < min_query_len:
            continue
        for acc2, seq2 in targets.items():
            # edit_distance = edlib_ed(seq1, seq2, mode="HW", task = "path", k = max_ed_threshold) # seq1 = query, seq2 = target
            edit_distance, cigar, locations = edlib_ed(seq1, seq2, mode="HW", task = "path", k = best_ed)
            # if acc1 == "transcript_447_support_14_15_1.57320042782e-125_21_5" and edit_distance < 100:
            #     print(edit_distance, cigar, locations, acc2)

            if edit_distance < 0:
                continue

            if 0 <= edit_distance <= best_ed:
                # if len(locations) > 1 or locations[0][0] > 100 or (len(seq2) - locations[0][1]) > 100:
                #     # print("Skipped:", cigar, locations, len(seq1), len(seq2))
                #     continue
                if edit_distance < best_ed:
                    best_edit_distances[acc1] = {}
                    best_cigars[acc1] = {}
                best_ed = edit_distance
                best_edit_distances[acc1][acc2] = best_ed
                best_cigars[acc1][acc2] = cigar
            elif edit_distance == best_ed:
                # if len(locations) > 1 or locations[0][0] > 100 or len(seq2) - locations[0][1] > 100:
                #     # print("Skipped:", cigar, locations, len(seq1), len(seq2))
                #     continue
                best_edit_distances[acc1][acc2] = best_ed
                best_cigars[acc1][acc2] = cigar

        i += 1
        # print(best_ed)
        # if acc1 == "transcript_447_support_14_15_1.57320042782e-125_21_5":
        #     print(best_ed, best_cigars[acc1])
        #     sys.exit()
                
    # for acc in best_edit_distances:
    #     if len(best_edit_distances[acc]) >1:
    #         print(best_edit_distances[acc])
    #     for acc2 in best_edit_distances[acc]:
    #         if best_edit_distances[acc][acc2] == 0:
    #             print("OK",acc,acc2)
    # print("HERE")
    # for acc in best_edit_distances:
    #     for acc2 in best_edit_distances[acc]:
    #         print( best_edit_distances[acc][acc2])

    return best_edit_distances, best_cigars

def get_ssw_alignments(best_edit_distances, querys, targets):
    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
    best_edit_distances_ssw = {}
    best_cigars_ssw = {}
    for acc1 in best_edit_distances:
        seq1 = querys[acc1]
        best_ed = len(seq1)
        best_edit_distances_ssw[acc1] = {}
        best_cigars_ssw[acc1] = {}
        for acc2 in best_edit_distances[acc1]:
            seq2 = targets[acc2]
            result = aligner.align(seq1, seq2, revcomp=False)
            seq2_aln, match_line, seq1_aln = result.alignment
            matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
            insertion_count = seq2_aln.count("-")
            deletion_count = seq1_aln.count("-")

            sw_ed = mismatches + indels
            best_edit_distances_ssw[acc1][acc2] =  sw_ed # (deletion_count, insertion_count, mismatches )
            seq1_aln, match_line, seq2_aln = result.alignment
            best_cigars_ssw[acc1][acc2] = (result.cigar, mismatches, indels, result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1 )

            # print(acc1,acc2)
            # print(result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1, result.cigar, mismatches, indels)
            # print()


            # print(sw_ed, (deletion_count, insertion_count, mismatches ))
            # print(seq1_aln)
            # print(match_line)
            # print(seq2_aln)
            # edit_distance, locations, cigar = edlib_traceback(seq1, seq2, k =1000)
            # print(edit_distance, locations, cigar)
            # print()            
    # for acc in best_cigars_ssw:
    #     if len(best_cigars_ssw[acc]) ==0:
    #         print("!!!!", acc)
    # print(len(best_cigars_ssw))
    # sys.exit()
    return best_edit_distances_ssw, best_cigars_ssw



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


def print_fully_supported_transcripts_per_gene_member_fasta(best_cigars_ssw, fully_supported_shared, outfolder):
    targeted = set(["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    targeted_dict = defaultdict(list)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    # outfile = open(os.path.join(args.outfolder, "all_supported_and_shared.fa"), "w")
    print(len(best_cigars_ssw))
    for c in best_cigars_ssw:
        if len(best_cigars_ssw[c]) == 0:
            print(c)
        # print(len(best_cigars_ssw[c]))
        if len(best_cigars_ssw[c]) > 1:
            # print(best_cigars_ssw[c])
            smallest_SNV = 100000
            smallest_other = 100000
            for t_cand in best_cigars_ssw[c]: # pick the one with smallest number of substitutions and single nucl indels
                cigar, mismatches, indels, q_start_softcl, q_end_softcl, t_start_softcl, t_end_softcl = best_cigars_ssw[c][t_cand]
                # print(cigar,mismatches, indels, t_start_softcl, t_end_softcl, c, t_cand )
                if mismatches + indels + q_start_softcl + q_end_softcl < smallest_SNV:
                    smallest_SNV = mismatches + indels
                    smallest_other = q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl
                    t_min = t_cand
                elif mismatches + indels == smallest_SNV:
                    if q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl <  smallest_other:
                        smallest_other = q_start_softcl + q_end_softcl + t_start_softcl + t_end_softcl
                        t_min = t_cand

            t = t_min
        else:
            t = best_cigars_ssw[c].keys()[0]

        c_is_targeted = False
        for target in targeted:
            if target in t:
                c_is_targeted = True
                targeted_dict[target].append(c)

        if not c_is_targeted:
            print("!!!",t)
            print(c, best_cigars_ssw[c])
            targeted_dict["NO_DB_ALN"].append(c) 

    counter = 0
    for gene_fam in targeted_dict:
        print("Nr transcripts in {0}: {1}".format(gene_fam, len(targeted_dict[gene_fam])))
        outfile = open(os.path.join(outfolder, gene_fam + ".fa"), "w")
        for c_acc in targeted_dict[gene_fam]:
            outfile.write(">{0}\n{1}\n".format(c_acc, fully_supported_shared[c_acc]))
            counter += 1
        outfile.close()

    print("Nr transcripts in {0}: {1}".format("NO_DB_ALN", len(targeted_dict["NO_DB_ALN"])))

    print("Total written to file:", counter)

def print_shared_with_illumina_support(illumina_support_files, shared_seqs, seq_to_acc_sample1, seq_to_acc_sample2, outfolder):
    illumina_supports_sample1 = {}
    illumina_supports_sample2 = {}
    outfile = open(os.path.join(outfolder, "shared_between_samples.tsv"), "w")
    for file_ in illumina_support_files:
        # print(file_)
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        sample, gene_fam = primer_to_family_and_sample[batch_size][primer_id]
        if sample == "sample1":
            for line in open(file_, "r"):
                acc, support = line.split("\t")                
                support = float(support)
                illumina_supports_sample1[acc] = (support, sample, gene_fam)
        if sample == "sample2":
            for line in open(file_, "r"):
                acc, support = line.split("\t")                
                support = float(support)
                illumina_supports_sample2[acc] = (support, sample, gene_fam)


    outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format("acc_sample1", "acc_sample2", "Illumina_support_sample1", "Illumina_support_sample2", "gene_family"))
    cnt, cnt1, cnt2, cnt_both = 0,0,0,0
    for seq in shared_seqs:
        acc1 = seq_to_acc_sample1[seq]
        acc2 = seq_to_acc_sample2[seq]
        support1, sample, gene_fam1 = illumina_supports_sample1[acc1]
        support2, sample, gene_fam2 = illumina_supports_sample2[acc2]
        assert gene_fam1 == gene_fam2
        outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(acc1, acc2, support1, support2, gene_fam1))


def add_if_full_illumina_support(illumina_support_files, table):
    for file_ in illumina_support_files:
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        sample_id, family_id = primer_to_family_and_sample[batch_size][primer_id]

        for line in open(file_, "r"):
            acc, ratio_supported =  line.split("\t") # nr_supported, total_length =
            results_iter = table.find(predicted_acc = acc, primer= primer_id, batch=batch_size)
            results = list(results_iter)
            assert len(results) == 1
            result = results[0]
            if sample_id == "sample1":
                if result["both_samples"] == "yes":
                    data = dict(id=result["id"], Illumina_support = nr_supported/float(total_length),  
                            full_supp_sample1 = None, full_supp_sample2 =None,)

                    table.update(data, ['id'])

                else:
                    pass
            if sample_id == "sample2":
                pass


    print(len(all_acc_sample1), len(all_acc_sample2))
    print(len(full_illumina_support_sample1), len(full_illumina_support_sample2))

    return full_illumina_support_sample1, full_illumina_support_sample2

def initialize_database(files):
    """
        predicted_acc  sample1 sample2  both  acc_sample1  acc_sample2  
        
        Illumina_support  full_supp_sample1 full_supp_sample2 family, primer, batch
        
        perfect_match_database  category  gene_member_number
    """
    db = dataset.connect('sqlite:///:memory:')
    table = db['predicted_transcripts']
    for file_ in files:
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        print(primer_id, batch_size, primer_to_family_and_sample[batch_size][primer_id])
        sample_id, family_id = primer_to_family_and_sample[batch_size][primer_id]

        for (acc, seq) in read_fasta(open(file_, 'r')):

            if sample_id == "sample1":
                data = dict(predicted_acc = acc, sequence = seq, sample1 = "yes", sample2 = "", both_samples = "",
                            acc_sample1 = acc,  acc_sample2 = "", Illumina_support = None,  
                            full_supp_sample1 = None, full_supp_sample2 =None, 
                            perfect_match_database = None, category = "-",  gene_member_number = "-",
                            primer = primer_id, batch = batch_size, family = family_id)
                table.insert(data)

            elif sample_id == "sample2":
                data = dict(predicted_acc = acc, sequence = seq, sample1 = "", sample2 = "yes", both_samples = "",
                            acc_sample1 = "",  acc_sample2 = acc, Illumina_support = None,  
                            full_supp_sample1 = None, full_supp_sample2 =None, 
                            perfect_match_database = None, category = "-",  gene_member_number = "-",
                            primer = primer_id, batch = batch_size, family = family_id)
                table.insert(data)
            else:
                print("BUG in reading DATABASE")
                sys.exit()

    print(len(table), "records read")
    return table

def add_if_shared_between_samples(table):
    for record in table.all():
        q_seq = record["sequence"]
        q_id = record["id"]
        results = list(table.find(sequence = q_seq))

        if len(results) == 1: # unique to one sample
            if record["sample1"] == "yes": 
                data = dict(id=q_id, sample2 = "-", acc_sample2 = "-", both_samples = "no" )
                table.update(data, ['id'])
                # print(list(table.find(id=q_id)))
            elif record["sample2"] == "yes":
                data = dict(id=q_id, sample1 = "-", acc_sample1 = "-", both_samples = "no" )
                table.update(data, ['id'])
            else:
                print("bug!")
                sys.exit()

        elif len(results) == 2: # in both samples
            s_acc1 = results[0]["predicted_acc"] if results[0]["predicted_acc"] == results[0]["acc_sample1"]  else  results[1]["predicted_acc"]
            s_acc2 = results[0]["predicted_acc"] if results[0]["predicted_acc"] == results[0]["acc_sample2"]  else  results[1]["predicted_acc"]
            assert s_acc1 != s_acc2

            data = dict(id=q_id, sample1 = "yes", sample2 = "yes", acc_sample1 = s_acc1, acc_sample2 = s_acc2, both_samples = "yes" )
            table.update(data, ['id'])                    
            print(list(table.find(id=q_id)))

        else:
            print("More than 2 hits woot??", record )
            sys.exit()

    sys.exit()




def main(params):

    # read in objects
    table = initialize_database(args.predicted)
    
    # check if shared between samples
    add_if_shared_between_samples(table)


    # add to database the illumina support illumina support 
    add_if_full_illumina_support(args.illumina_support_files, table)



    # print("{0} predicted transcripts with full illumina support.".format(predicted_with_full_illumina_support) )

    # get shared between samples
    seq_to_acc_sample1 = {}
    seq_to_acc_sample2 = {}
    sample1_dict = {}
    sample2_dict = {}



    print(len(sample1_dict), len(sample2_dict))
    shared_seqs = set(sample1_dict.values()) &  set(sample2_dict.values())
    only_sample1_seqs = set(sample1_dict.values()) - set(sample2_dict.values())
    only_sample2_seqs = set(sample2_dict.values()) - set(sample1_dict.values())

    # shared_seqs = sample1_seqs & sample2_seqs
    # only_sample1_seqs = sample1_seqs - sample2_seqs
    # only_sample2_seqs = sample2_seqs - sample1_seqs
    print("Total shared between samples:", len(shared_seqs))
    print("Only sample1:", len(only_sample1_seqs))
    print("Only sample2:", len(only_sample2_seqs))
    print("total for both samples:", len(sample1_dict.values()) + len(sample2_dict.values()) )

    print_shared_with_illumina_support(args.illumina_support_files, shared_seqs, seq_to_acc_sample1, seq_to_acc_sample2, args.outfolder)

    groups = {}
    shared = {}
    shared1 = {}
    shared2 = {}
    cnt, cnt1, cnt2, cnt_both = 0,0,0,0
    for seq in shared_seqs:
        acc1 = seq_to_acc_sample1[seq]
        acc2 = seq_to_acc_sample2[seq]
        if acc1 in full_illumina_support_sample1 and acc2 in full_illumina_support_sample2:
            shared1[acc1] = seq
            shared2[acc2] = seq
            shared[acc1] = seq
            print("HERE")
            cnt +=1
        elif acc1 in full_illumina_support_sample1:
            print("shared and full support in sample1 but not fully spported in sample2!", acc2)
            # print(acc1, acc2, full_illumina_support_sample2[acc2])
            cnt1 += 1
        elif acc2 in full_illumina_support_sample2:
            print("shared and full support in sample2 but not fully spported in sample1!", acc1)
            cnt2 += 1
        else:
            print("shared but not fully spported in either smaple!")
            cnt_both += 1

    groups["shared"] = shared
    print("shared and full illumina support in both samples:", cnt)
    print("shared, but full illumina support only in sample1:", cnt1)
    print("shared, but full illumina support only in sample2:", cnt2)
    print("No full support in either sample:", cnt_both)
    only_sample1 = {}
    for seq in only_sample1_seqs:
        acc1 = seq_to_acc_sample1[seq]
        if acc1 in full_illumina_support_sample1:
            only_sample1[acc1] = seq
    groups["only_sample1"] = only_sample1

    only_sample2 = {} 
    for seq in only_sample2_seqs:
        acc1 = seq_to_acc_sample2[seq]
        if acc1 in full_illumina_support_sample2:
            only_sample2[acc1] = seq
    groups["only_sample2"] = only_sample2

    # outfile = open(os.path.join(args.outfolder, "all_supported_and_shared.fa"), "w")
    # already_sampled = set()
    # full_support_shared = 0
    # fully_supported_shared = {}
    # for file_ in args.illumina_support_files:
    #     for line in open(file_, "r"):
    #         is_new = False
    #         acc, support = line.split("\t")
    #         support = float(support)
    #         if support < 1.0:
    #             continue

    #         if acc in shared1:
    #             is_new = True
    #             seq = shared1[acc]
    #         if acc in shared2:
    #             is_new = True
    #             seq = shared2[acc]

    #         if is_new and (seq not in already_sampled):
    #             full_support_shared += 1
    #             fully_supported_shared[acc] = seq
    #             outfile.write(">{0}\n{1}\n".format(acc, seq))
    #             already_sampled.add(seq)
    # outfile.close()
    # print("{0} transcripts had full support and shared betweeen samples".format(full_support_shared))
    # print("{0} transcripts had full support and shared betweeen samples".format(len(fully_supported_shared)))


    # classify the fully illumina supported transcripts to their best hit to database: annotate with best hit gene member

    database = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.database, 'r'))}
    database = {acc: seq.upper() for (acc, seq) in database.items() if "UNAVAILABLE" not in seq }

    print(only_sample1["transcript_2_support_3_3_not_tested_3_-1"])

    for group in groups:
        print()
        print(group.upper())
        print()
        outfolder  = os.path.join(args.outfolder, group)
        minimizer_graph_c_to_t, best_cigars = get_minimizers_2set_simple(groups[group], database, 32)
        minimizer_graph_x_to_c, best_cigars_ssw = get_ssw_alignments(minimizer_graph_c_to_t, groups[group], database)

        # print fasta file with all the transcripts that passed the criteria separated into gene families
        print_fully_supported_transcripts_per_gene_member_fasta(best_cigars_ssw, groups[group], outfolder)   

        #print a tsv with the best db hit for each transcript and if the hit was perfect or not.
        # print_fully_supported_transcripts_per_gene_member_tsv(best_cigars_ssw, groups[group], args)   




def parse_cigar(cigar):
    # print(cigar)

    tuples = []
    result = re.split(r'[=DXSMI]+', cigar)
    # print(result)
    i = 0
    for length in result[:-1]:
        i += len(length)
        type_ = cigar[i]
        i += 1
        tuples.append((length, type_ ))
    # print(tuples)
    for length, type_ in tuples:
        if int(length) >= 10 and type_ != "M":
            # print("has exon") 
            return True
    return False

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
    parser.add_argument('--predicted', type=str, nargs="+", help='Path to consensus fasta file(s)')
    parser.add_argument('--illumina_support_files', type=str, nargs="+", help='Path to tsv files with illumina support')
    parser.add_argument('--database', type=str, help='Path to database tsv file')
    parser.add_argument('--outfile', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')
    
    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    args.outfolder = path_

    main(args)