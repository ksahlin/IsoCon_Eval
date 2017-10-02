
from __future__ import print_function
import os,sys
import argparse
import re
import itertools

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
import pickle

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
    best_clip_distances = {}
    i = 1
    for acc1, seq1 in querys.items():
        if i % 200 == 0:
            print("processing candidate", i)
        best_ed = 10
        best_clipp = 10000000
        best_edit_distances[acc1] = {}
        best_clip_distances[acc1] = {}
        if len(seq1) < min_query_len:
            continue
        for acc2, seq2 in targets.items():
            # edit_distance = edlib_ed(seq1, seq2, mode="HW", task = "path", k = max_ed_threshold) # seq1 = query, seq2 = target
            edit_distance, cigar, locations = edlib_ed(seq1, seq2, mode="HW", task = "path", k = 0)

            if edit_distance < 0 or edit_distance > best_ed:
                continue

            if len(locations) > 1:
                print("ambiguous locations exiting:", cigar, locations, acc1, acc2)
                print("q:", seq1)
                print("t:", seq2)
                sys.exit("Ambiguous best alignment!")

            if locations[0][0] > 100 or (len(seq2) - locations[0][1]) > 100:
                print("Skipped:", cigar, locations, len(seq1), len(seq2))
                continue

            # if here we have a perfect alignment that is unamgiguous (only one "location", i.e. start and stop point in best cigar) with less than 100 bp clipped in each end
            if 0 <= edit_distance < best_ed:
                best_edit_distances[acc1] = {}
                best_ed = edit_distance
                best_edit_distances[acc1][acc2] = best_ed

                best_clip_distances[acc1] = {}
                best_clipp = locations[0][0] + (len(seq2) - locations[0][1])
                best_clip_distances[acc1][acc2] = best_clipp

            elif edit_distance == best_ed:
                if locations[0][0] + (len(seq2) - locations[0][1]) < best_clipp:
                    best_edit_distances[acc1] = {}
                    best_ed = edit_distance
                    best_edit_distances[acc1][acc2] = best_ed

                    best_clip_distances[acc1] = {}
                    best_clipp = locations[0][0] + (len(seq2) - locations[0][1])
                    best_clip_distances[acc1][acc2] = best_clipp

                elif locations[0][0] + (len(seq2) - locations[0][1]) == best_clipp:
                    best_edit_distances[acc1][acc2] = best_ed
                    best_clip_distances[acc1][acc2] = best_clipp
                    # print("identical everything",acc2, best_clipp, best_ed)
                else:
                    pass
                    # print("Worse clipp distance:", acc2, locations[0][0] + (len(seq2) - locations[0][1]), best_clipp)

        i += 1
    # for a in best_edit_distances:
    #     if len(best_edit_distances[a]) > 1:
    #         print(a, best_edit_distances[a])
    return best_edit_distances, best_clip_distances



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



def add_if_full_illumina_support(illumina_support_files, table):
    for file_ in illumina_support_files:
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        sample_id, family_id = primer_to_family_and_sample[batch_size][primer_id]

        for line in open(file_, "r"):
            if len(line.split("\t")) == 3:
                acc, nr_supported, total_length =  line.split("\t")
                ratio_supported = float(nr_supported) / float(total_length)
                # print(nr_supported, nr_supported, ratio_supported)
            elif len(line.split("\t")) == 2:
                acc, ratio_supported =  line.split("\t") # nr_supported, total_length =
                ratio_supported = float(ratio_supported)

            record_iter = table.find(predicted_acc = acc, primer= primer_id, batch=batch_size)
            records = list(record_iter)
            assert len(records) == 1
            record = records[0]

            if sample_id == "sample1":
                data = dict(id=record["id"], Illumina_support = ratio_supported,  
                        full_supp_sample1 = "yes" if ratio_supported == 1.0 else "no")
                table.update(data, ['id'])

            elif sample_id == "sample2":
                data = dict(id=record["id"], Illumina_support = ratio_supported,  
                        full_supp_sample2 = "yes" if ratio_supported == 1.0 else "no")
                table.update(data, ['id'])


    print("Full Illumina support:", len(list(table.find(Illumina_support=1.0))))
    
    for record in table.all():
        q_seq = record["sequence"]
        q_id = record["id"]
        results = list(table.find(sequence = q_seq))

        if len(results) == 1: # unique to one sample
            if record["sample1"] == "yes": 
                data = dict(id=q_id, full_supp_sample2 = "-")
                table.update(data, ['id'])
                # print(list(table.find(id=q_id)))
            elif record["sample2"] == "yes":
                data = dict(id=q_id, full_supp_sample1 = "-" )
                table.update(data, ['id'])
            else:
                print("bug!")
                sys.exit()

        elif len(results) == 2: # in both samples
            is_full_supp_sample1 = results[0]["full_supp_sample1"] if results[0]["predicted_acc"] == results[0]["acc_sample1"]  else results[1]["full_supp_sample1"]
            is_full_supp_sample2 = results[0]["full_supp_sample2"] if results[0]["predicted_acc"] == results[0]["acc_sample2"]  else results[1]["full_supp_sample2"]
            assert is_full_supp_sample1 != None
            assert is_full_supp_sample2 != None
            if record["predicted_acc"] == record["acc_sample1"]: # record in sample1
                data = dict(id=q_id, full_supp_sample2 = is_full_supp_sample2 )
                table.update(data, ['id'])                    

            elif record["predicted_acc"] == record["acc_sample2"]: # record in sample2
                data = dict(id=q_id, full_supp_sample1 = is_full_supp_sample1 )
                table.update(data, ['id'])                    
            else:
                print("bug..")
                sys.exit()
        else:
            print("More than 2 hits woot??", record )
            sys.exit()

    print("Shared in samples and full Illumina support in both samples:", len(list(table.find(full_supp_sample1="yes", full_supp_sample2 = "yes")))/2)
    print("Shared in samples and full Illumina support only in sample1:", len(list(table.find(full_supp_sample1="yes", full_supp_sample2 = "no")))/2)
    print("Shared in samples and full Illumina support only in sample2:", len(list(table.find(full_supp_sample1="no", full_supp_sample2 = "yes")))/2)
    print("Shared in samples and no full Illumina support in either:", len(list(table.find(full_supp_sample1="no", full_supp_sample2 = "no")))/2)

    print("Sample1 specific and full Illumina support:", len(list(table.find(full_supp_sample1="yes", full_supp_sample2 = "-"))))
    print("Sample2 specific and full Illumina support:", len(list(table.find(full_supp_sample1="-", full_supp_sample2 = "yes"))))


def initialize_database(files):
    """
        predicted_acc  sample1 sample2  both  acc_sample1  acc_sample2  
        
        Illumina_support  full_supp_sample1 full_supp_sample2 family, primer, batch
        
        perfect_match_database  category  gene_member_number
    """
    db = dataset.connect('sqlite:///:memory:')
    table = db['predicted_transcripts']
    for file_ in files:
        print(file_)
        primer_id = int(file_.split("/")[-2])
        batch_size = file_.split("/")[-3].split("_polished")[0]
        # print(primer_id, batch_size, primer_to_family_and_sample[batch_size][primer_id])
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
                print("BUG in reading references")
                sys.exit()

    print(len(table), "records read")
    return db

def add_if_shared_between_samples(table):
    for record in table.all():
        q_seq = record["sequence"]
        q_id = record["id"]
        results = list(table.find(sequence = q_seq))

        if len(results) == 1: # unique to one sample
            if record["sample1"] == "yes": 
                data = dict(id=q_id, sample2 = "No", acc_sample2 = "-", both_samples = "no" )
                table.update(data, ['id'])
                # print(list(table.find(id=q_id)))
            elif record["sample2"] == "yes":
                data = dict(id=q_id, sample1 = "No", acc_sample1 = "-", both_samples = "no" )
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
            # print(list(table.find(id=q_id)))

        else:
            print("More than 2 hits woot??", record )
            sys.exit()

    
    print("Shared:", len(list(table.find(both_samples="yes"))))
    print("Sample1:", len(list(table.find(sample2="No"))))
    print("Sample2:", len(list(table.find(sample1="No"))))


def add_if_perfect_match_in_database(references, table):
    references = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.references, 'r'))}
    references = {acc: seq.upper() for (acc, seq) in references.items() if "UNAVAILABLE" not in seq }

    queries = { (record["predicted_acc"], record["primer"] ) : str(record["sequence"]).upper() for record in table.all()}

    minimizer_graph_c_to_t, best_clip_distances = get_minimizers_2set_simple(queries, references, 100)
    # minimizer_graph_x_to_c, best_cigars_ssw = get_ssw_alignments(minimizer_graph_c_to_t, queries, references)

    for c_acc, primer_id in  minimizer_graph_c_to_t:
        if minimizer_graph_c_to_t[(c_acc, primer_id)]:
            all_best_hits = sorted(minimizer_graph_c_to_t[(c_acc, primer_id)].keys())
            ed = list(minimizer_graph_c_to_t[(c_acc, primer_id)].values())[0]
            # best_clip_distance = list(best_clip_distances[(c_acc, primer_id)].values())[0]
            if ed == 0:
                results = list(table.find(predicted_acc=c_acc, primer = primer_id ))
                assert len(results) == 1
                record = results[0]
                data = dict(id=record["id"], perfect_match_database = "yes")
                table.update(data, ['id']) 

    for record in table.all():
        if record["perfect_match_database"] != "yes":
            data = dict(id=record["id"], perfect_match_database = "no")
            table.update(data, ['id']) 

    # for rec in table.all():
    #     print(rec["perfect_match_database"])
    # print("Total number of perfect matches to reference database:", len(list(table.find(perfect_match_database="yes"))))
    print("Shared between samples and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes")))/2)
    print("Subcategories:")
    print("Shared between samples, full illumina support in both, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "yes", full_supp_sample2 = "yes" )))/2)
    print("Shared between samples, full illumina support in sample1 only, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "yes", full_supp_sample2 = "no" )))/2)
    print("Shared between samples, full illumina support in sample2 only, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "no", full_supp_sample2 = "yes" )))/2)
    print("Shared between samples, no full illumina support in either sample but perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "no", full_supp_sample2 = "no" )))/2)

    print()
    print("Sample specific and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "no"))))    
    print("Subcategories:")
    print("Sample 1 specific, full illumina support, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "no", full_supp_sample1 = "yes"))))
    print("Sample 2 specific, full illumina support, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "no", full_supp_sample2 = "yes"))))
    print("Sample 1 specific, no illumina support, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "no", full_supp_sample1 = "no"))))
    print("Sample 2 specific, no illumina support, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "no", full_supp_sample2 = "no"))))
    print()

    for record in table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "no", full_supp_sample2 = "no" ):
        print("not supp in any:", record["predicted_acc"], record["primer"], record["batch"],  record["family"])
        print(record["sequence"])
    for record in table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "yes", full_supp_sample2 = "no" ):
        print("not supp in sample2:", record["predicted_acc"], record["primer"], record["batch"],  record["family"])
    for record in table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "no", full_supp_sample2 = "yes" ):
        print("not supp in sample1:", record["predicted_acc"], record["primer"], record["batch"],  record["family"])

    # print("Total number of predictions that does not have a perfect match to  reference database:", len(list(table.find(perfect_match_database="no"))))


def print_database_to_tsv(db):
    column_order = ['id',  'predicted_acc', 'family', 'Illumina_support', 'perfect_match_database', 'sample1', 'sample2', 'acc_sample1', 'acc_sample2', 'both_samples', 'primer', 'category', 'full_supp_sample1', 'full_supp_sample2', 'gene_member_number', 'batch', 'sequence',]

    filename = open(args.outfile, "w")
    filename.write("\t".join(column_order) + "\n")
    for record in db['predicted_transcripts'].all():
        # print(record, type(record))
        # print(record.values(), type(record))
        record_values = [str(record[field]) for field in column_order]
        fmt = "\t".join(record_values)
        fmt += "\n"
        filename.write(fmt)


def print_database_to_separate_fastas(db):
    targeted = set(["BPY", "CDY1", "CDY2", "DAZ", "HSFY1", "HSFY2", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    targeted_outfiles_dict = {}

    # create all files
    for group in ["both_samples", "sample1_specific", "sample2_specific"]:
        outfolder  = os.path.join(args.outfolder, group)
        targeted_outfiles_dict[group] = {}
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        for fam in targeted:
            if group == "both_samples":
                targeted_outfiles_dict[group][fam] = {}
                group_fam_file1 = open(os.path.join(outfolder, fam + "_sample1.fa"), "w")
                targeted_outfiles_dict[group][fam]["sample1"] = group_fam_file1
                group_fam_file2 = open(os.path.join(outfolder, fam + "_sample2.fa"), "w")
                targeted_outfiles_dict[group][fam]["sample2"] = group_fam_file2
            else:
                group_fam_file = open(os.path.join(outfolder, fam + ".fa"), "w")
                targeted_outfiles_dict[group][fam] = group_fam_file


    # print all perfect illumina predictions
    cnt_both = 0
    cnt_s1 = 0
    cnt_s2 = 0
    for record in db['predicted_transcripts'].all():
        family = record["family"]
        if record["both_samples"] == "yes":
            group = "both_samples"
            # is_fully_supported = (record["full_supp_sample1"] == record["full_supp_sample2"] == "yes")

        elif record["sample1"] == "yes":
            group = "sample1_specific"
            # is_fully_supported = (record["full_supp_sample1"] == "yes")

        elif record["sample2"] == "yes":
            group = "sample2_specific"
            # is_fully_supported = (record["full_supp_sample2"] == "yes")
        else:
            print("BUG, no category!")
            sys.exit()

        if True: # all transcripts, not just full illumina support is_fully_supported:
            seq = str(record["sequence"])
            acc = str(record["predicted_acc"])
            is_sample1 = True if record["predicted_acc"] == record["acc_sample1"] else False
            
            if record["both_samples"] == "yes":
                if is_sample1:
                    outfile = targeted_outfiles_dict[group][family]["sample1"]
                    outfile.write(">{0}\n{1}\n".format(acc, seq))
                    cnt_both += 1

                else:
                    outfile = targeted_outfiles_dict[group][family]["sample2"]
                    outfile.write(">{0}\n{1}\n".format(acc, seq))

            else:
                if is_sample1:
                    cnt_s1 += 1
                else:
                    cnt_s2 += 1

                outfile = targeted_outfiles_dict[group][family]
                outfile.write(">{0}\n{1}\n".format(acc, seq))
            

    print("Both samples written to file:", cnt_both)
    print("Sample1 written to file:", cnt_s1)
    print("Sample2 written to file:", cnt_s2)

    # remove empty files and close rest
    for group in targeted_outfiles_dict:
        for fam in targeted_outfiles_dict[group]:
            if group == "both_samples":
                for sample in ["sample1", "sample2"]:
                    outfile = targeted_outfiles_dict[group][fam][sample]
                    if outfile.tell() > 0:
                        outfile.close()
                    else:
                        outfile.close()
                        os.remove(outfile.name) 

            else:
                outfile = targeted_outfiles_dict[group][fam]
                if outfile.tell() > 0:
                    outfile.close()
                else:
                    outfile.close()
                    os.remove(outfile.name)                 


def find_longest_ORFS(table, references):
    pattern = r'ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)'
    references = {acc: seq.upper() for (acc, seq) in  read_fasta(open(args.references, 'r'))}
    references = {acc: seq.upper() for (acc, seq) in references.items() if "UNAVAILABLE" not in seq }

    for acc,seq in references.items():
        all_orfs = re.findall(pattern, seq)
        if all_orfs:
            longest_orf = max(all_orfs, key= len)
            print(len(longest_orf), len(seq), seq.index(longest_orf))
        else:
            print("no premature stop codon found, counting as protein coding")

    # for record in table.all():
    #     seq = record["sequence"]
    #     all_orfs = re.findall(pattern, seq)
    #     if all_orfs:
    #         longest_orf = max(all_orfs, key= len)
    #         print(len(longest_orf), len(seq), seq.index(longest_orf))
    #     else:
    #         print("no premature stop codon found, counting as protein coding")

def ssw_alignment(x, y, ends_discrepancy_threshold = 250 ):
    """
        Aligns two sequences with SSW
        x: query
        y: reference

    """

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=0, matrix=score_matrix)

    # for the ends that SSW leaves behind
    bio_matrix = matlist.blosum62
    g_open = -1
    g_extend = -0.5
    ######################################

    # result = aligner.align("GA", "G", revcomp=False)
    # y_alignment, match_line, x_alignment = result.alignment
    # c = Counter(match_line)
    # matches, mismatches, indels = c["|"], c["*"], c[" "]
    # alignment_length = len(match_line)
    # print("matches:{0}, mismatches:{1}, indels:{2} ".format(matches, mismatches, indels))
    # print(match_line)

    result = aligner.align(x, y, revcomp=False)
    y_alignment, match_line, x_alignment = result.alignment
    # print()
    # print(y_alignment)
    # print(match_line)
    # print(x_alignment)
    matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")

    # alignment_length = len(match_line)
    
    start_discrepancy = max(result.query_begin, result.reference_begin)  # 0-indexed # max(result.query_begin, result.reference_begin) - min(result.query_begin, result.reference_begin)
    query_end_discrepancy = len(x) - result.query_end - 1
    ref_end_discrepancy = len(y) - result.reference_end - 1
    end_discrepancy = max(query_end_discrepancy, ref_end_discrepancy)  # max(result.query_end, result.reference_end) - min(result.query_end, result.reference_end)
    # print("disc:", start_discrepancy, end_discrepancy)
    tot_discrepancy = start_discrepancy + end_discrepancy

    if 0 < start_discrepancy <= ends_discrepancy_threshold:
        print("HERE",start_discrepancy)
        matches_snippet = 0
        mismatches_snippet = 0
        if result.query_begin and result.reference_begin:
            query_start_snippet = x[:result.query_begin]
            ref_start_snippet = y[:result.reference_begin]
            alns = pairwise2.align.globalds(query_start_snippet, ref_start_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            # print(alns)
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            # print(matches_snippet, mismatches_snippet, indels_snippet)
            query_start_alignment_snippet = top_aln[0]
            ref_start_alignment_snippet = top_aln[1]
        elif result.query_begin:
            query_start_alignment_snippet = x[:result.query_begin]
            ref_start_alignment_snippet = "-"*len(query_start_alignment_snippet)
            indels_snippet = len(ref_start_alignment_snippet)
        elif result.reference_begin:
            ref_start_alignment_snippet = y[:result.reference_begin]
            query_start_alignment_snippet = "-"*len(ref_start_alignment_snippet)
            indels_snippet = len(query_start_alignment_snippet)
        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        # print(ref_start_alignment_snippet)
        # print(query_start_alignment_snippet)
        y_alignment = ref_start_alignment_snippet + y_alignment
        x_alignment = query_start_alignment_snippet + x_alignment

    if 0 < end_discrepancy <= ends_discrepancy_threshold:
        print("HERE2", end_discrepancy)
        matches_snippet = 0
        mismatches_snippet = 0
        if query_end_discrepancy and ref_end_discrepancy:
            query_end_snippet = x[result.query_end+1:]
            ref_end_snippet = y[result.reference_end+1:]
            alns = pairwise2.align.globalds(query_end_snippet, ref_end_snippet, bio_matrix, g_open, g_extend)
            top_aln = alns[0]
            mismatches_snippet = len(list(filter(lambda x: x[0] != x[1] and x[0] != '-' and x[1] != "-", zip(top_aln[0],top_aln[1]))))
            indels_snippet = top_aln[0].count("-") + top_aln[1].count("-")
            matches_snippet = len(top_aln[0]) - mismatches_snippet - indels_snippet
            query_end_alignment_snippet = top_aln[0]
            ref_end_alignment_snippet = top_aln[1]
        elif query_end_discrepancy:
            query_end_alignment_snippet = x[result.query_end+1:]
            ref_end_alignment_snippet = "-"*len(query_end_alignment_snippet)
            indels_snippet = len(ref_end_alignment_snippet)

        elif ref_end_discrepancy:
            ref_end_alignment_snippet = y[result.reference_end+1:]
            query_end_alignment_snippet = "-"*len(ref_end_alignment_snippet)
            indels_snippet = len(query_end_alignment_snippet)

        else:
            print("BUG")
            sys.exit()
        matches, mismatches, indels = matches + matches_snippet, mismatches + mismatches_snippet, indels + indels_snippet

        y_alignment = y_alignment + ref_end_alignment_snippet
        x_alignment = x_alignment + query_end_alignment_snippet

    return x_alignment, y_alignment, matches, mismatches, indels       


from networkx import nx


def create_isoform_graph(transcripts):
    
    G = nx.Graph()
    for acc, primer_id in  transcripts.keys():
        G.add_node(int(acc.split("_")[1]), accession = acc)
    # nodes = [ (int(acc.split("_")[1]), {"name": acc}) for acc, primer_id in  transcripts.keys()]
    # G.add_nodes_from(nodes)

    # family_members = defaultdict(set)
    # family_members_alignments = defaultdict(dict)

    score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
    aligner = ssw.Aligner(gap_open=2, gap_extend=0, matrix=score_matrix)
    cntr = 0

    processed = set()
    # already_assigned = set()
    for acc1, seq1 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
        cntr += 1
        processed.add(acc1)
        # print("length t:", len(seq1), acc1)
        if cntr % 5 == 0:
            print(cntr, "sequences processed")

        # if acc1 in already_assigned:
        #     # print("allready assigned to larger sequence!")
        #     continue

        for acc2, seq2 in sorted(transcripts.items(), key= lambda x: len(x[1]), reverse=True):
            if acc2 in processed:
                continue
            # if acc1 == acc2:
            #     continue
            # if seq1 == seq2:
            #     continue

            # result = aligner.align(seq1, seq2, revcomp=False)
            # seq2_aln, match_line, seq1_aln = result.alignment

            seq1_aln, seq2_aln, matches, mismatches, indels = ssw_alignment(seq1, seq2, ends_discrepancy_threshold = 2000 )
            
            del_seq1 = re.findall(r"[-]+",seq1_aln)
            del_seq2 = re.findall(r"[-]+",seq2_aln)
            mismatches = len([ 1 for n1, n2 in zip(seq1_aln,seq2_aln) if n1 != n2 and n1 != "-" and n2 != "-" ])

            # print(del_seq1)
            # print(del_seq2)
            # matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")

            # if acc1 == ('transcript_167_support_20_27_6.43088670805209e-39_264_1', 2) and acc2 == ('transcript_461_support_7_8_2.1670307017329207e-22_32_1', 2) :
            #     print(len(seq1), len(seq2), len(seq1_aln), len(seq2_aln))
            #     print(matches, mismatches, indels)
            #     print(del_seq1, del_seq2)
            #     print(seq2_aln)
            #     # print(match_line)
            #     print(seq1_aln)
            #     # print(result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1)            

            # alignment can bee different in ends
            # if (result.query_begin > 0 and result.reference_begin > 0) or ( (len(seq1) - result.query_end - 1) > 0 and (len(seq2) - result.reference_end -1)  > 0):
            #     res = edlib.align(seq1, seq2, task="path")          
            #     print(res)

            # print(match_line)
            # print(acc1, acc2)
            # print(result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1)            
            # print(seq2_aln)
            # print(match_line)
            # print(seq1_aln)
            # by default (since all transcripts are distinct if we end up here), each transcript is its on gene member
            # if we find an alingment that contains only structural changes of > X (10) nucleotides, and no other smaller differences we classify as same family
            if mismatches == 0:
                # print(del_seq1)
                # print(del_seq2)
                # if acc2 == ("transcript_20_support_3_3_9.883280159550703e-07_6_1", 2):
                #     print(len(seq1), len(seq2), len(seq1_aln), len(seq2_aln))
                #     print(matches, mismatches, indels)
                #     print(del_seq1, del_seq2)
                #     print(seq1_aln)
                #     print(seq2_aln)
                    # sys.exit()

                del_lengths1 = [len(del_) for del_ in del_seq1]
                del_lengths2 = [len(del_) for del_ in del_seq2]
                no_small_del_in_seq1 = ((len(del_lengths1) > 0 and min(del_lengths1) >= 3) or len(del_lengths1)  == 0)
                no_small_del_in_seq2 = ((len(del_lengths2) > 0 and min(del_lengths2) >= 3) or len(del_lengths2)  == 0)
                # print(no_small_del_in_seq1, no_small_del_in_seq2)
                # print((len(del_seq1) > 0 and min(del_seq1) >= 3), len(del_seq1)  == 0)
                if no_small_del_in_seq1 and no_small_del_in_seq2:
                    G.add_edge(int(acc1[0].split("_")[1]), int(acc2[0].split("_")[1]), 
                                alignment={ int(acc1[0].split("_")[1]) : seq1_aln,  int(acc2[0].split("_")[1]) : seq2_aln })

                    # print("Isoform!")
                    # seq1 is that larger reference
                    # if len(del_seq1) == 0:
                    #     family_members[acc1].add(acc2)
                    #     already_assigned.add(acc2)
                    #     family_members_alignments[acc1][acc1] = seq1_aln
                    #     family_members_alignments[acc1][acc2] = seq2_aln

                    # elif len(del_seq2) == 0:
                    #     print(len(seq1), len(seq2), len(seq1_aln), len(seq2_aln))
                    #     print(matches, mismatches, indels)
                    #     print(del_seq1, del_seq2)
                    #     print(seq2_aln)
                    #     print(seq1_aln)
                    #     print("BUG! seq2 larger!, Why did we reach here")
                    #     sys.exit()

                    # else:
                    #     pass
                    #     # print(len(seq1), len(seq2), len(seq1_aln), len(seq2_aln))
                    #     # print(matches, mismatches, indels)
                    #     # print(del_seq1, del_seq2)
                    #     # print(seq1_aln)
                    #     # print(seq2_aln)
                    #     # print("Isoforms that both contain missing exons (smaller isoform might be assigned later though)!", acc1, acc2)
                    #     # check if isoform acc2 has been added to previous sequence

                    #     # sys.exit()

                else:
                    pass
                    # print("Different only by small indel!!")

        # if acc1 not in family_members:
        #     # print(acc1, "is not assigned to larger isoform and also had no isoforms!")
        #     family_members[acc1] = set()
        #     family_members_alignments[acc1][acc1] = seq1


        # print("Nr member isoforms:", len(family_members[acc1]) + 1)
            # matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
            # insertion_count = seq2_aln.count("-")
            # deletion_count = seq1_aln.count("-")

            # print()
            # sw_ed = mismatches + indels
            # best_edit_distances_ssw[acc1][acc2] =  sw_ed # (deletion_count, insertion_count, mismatches )
            # seq1_aln, match_line, seq2_aln = result.alignment
            # best_cigars_ssw[acc1][acc2] = (result.cigar, mismatches, indels, result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1 )

    list_of_maximal_cliques = list(nx.find_cliques(G))
    print("Number of possible members:", len(list_of_maximal_cliques) )
    print("clique sizes", [ len(cl) for cl in  sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True)] )
    for cl in  sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True):
        print()
        print("NEXT CLIQUE")
        for isoform1, isoform2 in itertools.combinations(cl, 2):
            print(len(G.edge[isoform1][isoform2]["alignment"][isoform1]), len(G.edge[isoform1][isoform2]["alignment"][isoform2]) )

    # print(G.edges(data=True))
    # return family_members, family_members_alignments, G
    return G


def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d   

# def is_equivalence_class(family_members, member):
#     isoforms_to_member = family_members[member]   # this is a set
#     print("Main set:", isoforms_to_member)

#     for candidate_member in family_members[member]:
#         temp_set = isoforms_to_member - set([candidate_member])
#         print("Isoform set:", temp_set)
#         isoforms_to_candidate_member = family_members[candidate_member]
#         if not temp_set.issubset(isoforms_to_candidate_member):
#             print("here why?")
#             return False

#     return True        

def construct_consensus_progressively(G, clique_id_to_accessions, transcripts, number_to_member_acc):
    pattern = re.compile(r"[-]+")
    consensi = {}
    for member_id, list_of_acc in clique_id_to_accessions.items():
        if len(list_of_acc) == 1:
            isoform = list_of_acc[0]
            consensi[member_id] = transcripts[ number_to_member_acc[isoform] ]
        else:
            consensi[member_id] = ""

    for member_id, list_of_acc in clique_id_to_accessions.items():
        if len(list_of_acc) > 1:
            isoform = list_of_acc[0]
            c_temp = transcripts[ number_to_member_acc[isoform] ]
            for isoform_id in list_of_acc[1:]:
                (acc, primer_id) = number_to_member_acc[ isoform_id ]
                isof = transcripts[(acc, primer_id)]
                isof_aln, c_temp_aln, matches, mismatches, indels = ssw_alignment(isof, c_temp, ends_discrepancy_threshold = 2000 )

                c_list = []
                for n1, n2 in zip(c_temp_aln, isof_aln):
                    if n1 != "-":
                        c_list.append(n1)
                    else:
                        c_list.append(n2)
                c_temp = "".join([ n for n in c_list ] )
            consensi[member_id] = c_temp
        
        print("LENGHT C", len(consensi[member_id]))


    # form_consensus_from_pairwise = {} # list of tuples
    # for clique_id in clique_id_to_accessions:
    #     print("starting gapped multialignmet creation for clique ", clique_id)
    #     longest_pairwise = 0
    #     for isoform1 in clique_id_to_accessions[clique_id]:
    #         for isoform2 in clique_id_to_accessions[clique_id]:
    #             if isoform1 == isoform2:
    #                 continue 

    #             isoform1_alnmnt = G.edge[isoform1][isoform2]["alignment"][isoform1]
    #             isoform2_alnmnt = G.edge[isoform1][isoform2]["alignment"][isoform2]
    #             print(len(isoform1_alnmnt))
    #             if len(isoform1_alnmnt) > longest_pairwise:
    #                 form_consensus_from_pairwise[clique_id] = (isoform1_alnmnt, isoform2_alnmnt)
    #                 longest_pairwise = len(isoform1_alnmnt) 

    # # print(form_consensus_from_pairwise)
    # # print(list(consensi.keys()))
    # for member_id in consensi:
    #     if not consensi[member_id]: # no consensi formed yet
    #         c_list = []
    #         for n1, n2 in zip(form_consensus_from_pairwise[member_id][0], form_consensus_from_pairwise[member_id][1]):
    #             if n1 != "-":
    #                 c_list.append(n1)
    #             else:
    #                 c_list.append(n2)
    #         c = "".join([ n for n in c_list ] )
    #         consensi[member_id] = c
    #     print("C",member_id,   len(consensi[member_id]))
    return consensi

def construct_multialignment_from_consensi(clique_id_to_accessions, transcripts, number_to_member_acc, consensi):
    multialignments = {}
    for member_id, list_of_acc in clique_id_to_accessions.items():
        multialignments[member_id] = {}
        # multialignments[member_id]["consensus"] = consensi[member_id]
        cons = consensi[member_id]
        for isoform_id in list_of_acc:
            (acc, primer_id) = number_to_member_acc[ isoform_id ]
            isof = transcripts[(acc, primer_id)]
            isof_aln, cons_aln, matches, mismatches, indels = ssw_alignment(isof, cons, ends_discrepancy_threshold = 2000 )
            multialignments[member_id][acc] = isof_aln
            # print(len(isof_aln), len(cons), len(isof)) 
            # print(isof_aln)
            # print(cons)
            assert len(isof_aln) == len(cons)

    return multialignments


# def construct_multialignment_from_pairwise_alignments(G, clique_id_to_accessions, transcripts, number_to_member_acc):
    # multialignments = { member_id : { acc : "" for acc in list_of_acc}  for member_id, list_of_acc in clique_id_to_accessions.items()}
    # pattern = re.compile(r"[-]+")

    # for clique_id in multialignments:
    #     print("starting gapped multialignmet creation for clique ", clique_id)
    #     for isoform1 in multialignments[clique_id]:
    #         deletion_positions_and_length = {isoform1 : {} }
    #         for isoform2 in multialignments[clique_id]:
    #             if isoform1 == isoform2:
    #                 continue 
    #             isoform1_alnmnt = G.edge[isoform1][isoform2]["alignment"][isoform1]
    #             # deletions_in_isoform = re.findall(pattern, isoform1_alnmnt)
    #             position_offset = 0
    #             for m in pattern.finditer(isoform1_alnmnt):
    #                 del_length = len(m.group())
    #                 pos = m.start() - position_offset
    #                 position_offset += del_length
    #                 # print( m.start(), m.group(), len(m.group()) )
    #                 if pos not in deletion_positions_and_length[isoform1]:
    #                     deletion_positions_and_length[isoform1][pos] = del_length
    #                 else:
    #                     if del_length > deletion_positions_and_length[isoform1][pos]: # add the longest stretch only
    #                         deletion_positions_and_length[isoform1][pos] = del_length


    #         # insert all gaps
    #         seq = transcripts[ number_to_member_acc[isoform1] ]
    #         seq_gapped = ""
    #         prev_pos = 0
    #         for pos in sorted(deletion_positions_and_length[isoform1].keys()):
    #             gap_length = deletion_positions_and_length[isoform1][pos]
    #             seq_gapped += seq[ prev_pos: pos] + "-"*gap_length
    #             # print(pos, gap_length, len(seq[ prev_pos: pos]), len("-"*gap_length))
    #             # print(len(seq_gapped))

    #             prev_pos = pos

    #         seq_gapped += seq[ prev_pos : ]             
    #         multialignments[clique_id][isoform1] = seq_gapped

    #         # print( len(seq),  len(seq_gapped), number_to_member_acc[isoform1] ) #, number_to_member_acc[isoform1])
    #         print(len(seq_gapped), len(seq))
    # # print(multialignments)
    # return multialignments


def get_gene_member_number(table, params):

    tsv_outfile = open(os.path.join(args.outfolder, "predicted_to_genemembers.tsv"), "w")
    tsv_outfile.write("FAMILY\tMEMBER_ID\tACCESSION\n")
    dot_graphs_folder = os.path.join(args.outfolder, "dot_graphs")
    mkdir_p(dot_graphs_folder)
    # doing the gene member analysis of only shared transcripts
    all_family_consensi = {}
    for target in ["BPY", "CDY1", "CDY2", "DAZ", "HSFY1", "HSFY2", "PRY", "RBMY", "TSPY", "XKRY", "VCY"]:
        family_records = table.find(table.table.columns.predicted_acc == table.table.columns.acc_sample1, family = target, both_samples = "yes" )
        accession_to_record_id = {}
        transcripts = {}
        print(target)
        for record in family_records:
            transcripts[(record["predicted_acc"], record["primer"] )] = record["sequence"]
            accession_to_record_id[(record["predicted_acc"], record["primer"] )] = record["id"]
        print(len(transcripts))

        if params.graph_prefix:
            G = pickle.load( open( os.path.join(params.graph_prefix, target + ".p"), "rb" ) )
        else:
            G = create_isoform_graph(transcripts)
            pickle.dump( G, open( os.path.join(params.outfolder, target + ".p"), "wb" ) )
            dot_graph_out = os.path.join(dot_graphs_folder, target + ".dot")
            nx.write_dot(G, dot_graph_out)

        number_to_member_acc = { int(acc.split("_")[1]) : (acc, primer_id) for acc, primer_id in  transcripts.keys()}        
        list_of_maximal_cliques = list(nx.find_cliques(G))
        clique_id_to_accessions = {member_id : cl for member_id, cl in enumerate(sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True))}
        print(clique_id_to_accessions)
        # multialignments = construct_multialignment_from_pairwise_alignments(G, clique_id_to_accessions, transcripts, number_to_member_acc)
        consensi = construct_consensus_progressively(G, clique_id_to_accessions, transcripts, number_to_member_acc)
        all_family_consensi[ target ] = consensi
        multialignments = construct_multialignment_from_consensi(clique_id_to_accessions, transcripts, number_to_member_acc, consensi)

        ## Write output 
        family_folder = os.path.join(args.outfolder, target)
        mkdir_p(family_folder)


        for member_id, cl in clique_id_to_accessions.items(): #enumerate(sorted(list_of_maximal_cliques, key= lambda x: len(x), reverse=True)):
            members_fasta = open(os.path.join(family_folder, str(member_id) + ".fa"), "w") 
            members_aligned_fasta = open(os.path.join(family_folder, str(member_id) + "_aligned.fa"), "w")
            members_aligned_fasta.write( ">{0}\n{1}\n".format("consensus", consensi[member_id]) )
            for transcript_number in cl:
                isoform, primer_id = number_to_member_acc[transcript_number]
                # print(len(G.edge[isoform1][isoform2]["alignment"][isoform1]), len(G.edge[isoform1][isoform2]["alignment"][isoform2]) )
                tsv_outfile.write("{0}\t{1}\t{2}\n".format(target, member_id, str(isoform)))
                members_fasta.write(">{0}\n{1}\n".format(isoform, transcripts[(isoform, primer_id)]) )
                members_aligned_fasta.write( ">{0}\n{1}\n".format(isoform, multialignments[member_id][isoform]) )



        all_processed = set([ acc for member_id, list_of_acc in clique_id_to_accessions.items() for acc in list_of_acc ])
        print("Number non processed sequences(bug if more than 0): ", len( set(transcripts.keys()).difference(all_processed)) )

    return all_family_consensi


def plot_member_ed(params, all_family_consensi):

    # print all pairwise distances to tsv file
    ## fam member_id1 member_id2 len(m1) len(m2) ed
    ed_tsv_file =open( os.path.join(params.outfolder, "member_edit_distances.tsv"), "w" )
    ed_tsv_file.write("family\tmember1\tmember2\tm1_len\tm2_len\ted\n")

    for fam in all_family_consensi:
        processed = set()
        for m1_id, m1_seq in all_family_consensi[fam].items():
            ed_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(fam, m1_id, m1_id, len(m1_seq), len(m1_seq), 0 ))
            processed.add(m1_id)
            for m2_id, m2_seq in all_family_consensi[fam].items():
                if m2_id not in processed:
                    result = edlib.align(m1_seq, m2_seq, task = "distance")
                    ed = result["editDistance"]
                    ed_tsv_file.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(fam, m1_id, m2_id, len(m1_seq), len(m2_seq), int(ed) ))

    ed_tsv_file.close()

    # plot
    dtypes =  {"family": str, "ed" : int, "member1" : int, "member2" : int, "m1_len" : int, "m2_len" : int}
    data = pd.read_csv(ed_tsv_file.name, sep="\t", dtype = dtypes)
    print(data.isnull().values.any())

    def facet_heatmap(data, **kws):
        data = data.pivot("member1", 'member2', 'ed')
        pd.set_option('display.max_rows', 1000)
        print(data)
        print(data.columns)

        mask = data.isnull()
        print(mask)
        for i in range(len(mask)):
            for j in range(len(mask[i])):
                if data[i ][j] > 99:
                    # print("true")
                    mask[i][j] = True

        sns.heatmap(data, cmap='coolwarm_r', annot=True, mask = mask, **kws)  # <-- Pass kwargs to heatmap fmt = "d"


    # with sns.plotting_context(font_scale=5.5):
    #     g = sns.FacetGrid(data, col="family", col_wrap=3, size=3, aspect=1 ) #, col_order = ["TSPY13P", "HSFY2", "DAZ2"])

    # cbar_ax = g.fig.add_axes([.92, .3, .02, .4])  # <-- Create a colorbar axes

    # g = g.map_dataframe(facet_heatmap,
    #                     cbar_ax=cbar_ax,
    #                     vmin=0, vmax=10)  # <-- Specify the colorbar axes and limits
    plot_outfolder = os.path.join(params.outfolder, "edit_distance_plots")
    mkdir_p(plot_outfolder)
    for target in ["BPY", "CDY1", "CDY2", "DAZ", "HSFY1", "HSFY2", "PRY", "RBMY", "TSPY", "XKRY", "VCY"]:
        df = data.loc[data['family'] == target]
        facet_heatmap(df, vmin=0, vmax=10)

        outfile = os.path.join(plot_outfolder, target + ".png")
        plt.savefig(outfile)
        plt.close()

    # g.set_titles(col_template="{col_name}", fontweight='bold', fontsize=14)
    # g.fig.subplots_adjust(right=.9)  # <-- Add space so the colorbar doesn't overlap the plot
    # outfile = os.path.join(params.outfolder, "member_edit_distances.png")

    # plt.savefig(outfile)
    # plt.close()


def main(params):

    # read in objects
    db = initialize_database(args.predicted)
    table = db["predicted_transcripts"]
    # print(table.columns)
    # check if shared between samples
    add_if_shared_between_samples(table)


    # add to references the illumina support illumina support 
    add_if_full_illumina_support(args.illumina_support_files, table)

    # add if perfect match or not
    add_if_perfect_match_in_database(args.references, table)

    # Find longest ORFs
    # NOT WORKING YET!
    # find_longest_ORFS(table, args.references)

    # classify into gene memebers
    all_family_consensi = get_gene_member_number(table, params)

    # plot edit distances of consensi within a family
    plot_member_ed(params, all_family_consensi)

    # print to file
    print_database_to_tsv(db)

    # print fasta
    print_database_to_separate_fastas(db)


    # result = db['predicted_transcripts'].all()
    # dataset.freeze(result, mode='item', format='csv', filename=args.outfile)





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
    subparsers = parser.add_subparsers(help='sub-command help')
    create = subparsers.add_parser('create_db', help='a help')
    load = subparsers.add_parser('load_db', help='a help')

    create.add_argument('--predicted', type=str, nargs="+", help='Path to consensus fasta file(s)')
    create.add_argument('--illumina_support_files', type=str, nargs="+", help='Path to tsv files with illumina support')
    create.add_argument('--references', type=str, help='Path to references tsv file')
    create.add_argument('--outfile', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')
    create.add_argument('--graph_prefix', type=str, help='already computed alignments.')
    
    load.add_argument('--database', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')

    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)

    if path_:
        mkdir_p(path_)
    args.outfolder = path_

    main(args)

