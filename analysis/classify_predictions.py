
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


# def get_ssw_alignments(best_edit_distances, querys, targets):
#     score_matrix = ssw.DNA_ScoreMatrix(match=1, mismatch=-2)
#     aligner = ssw.Aligner(gap_open=2, gap_extend=1, matrix=score_matrix)
#     best_edit_distances_ssw = {}
#     best_cigars_ssw = {}
#     for acc1 in best_edit_distances:
#         seq1 = querys[acc1]
#         best_ed = len(seq1)
#         best_edit_distances_ssw[acc1] = {}
#         best_cigars_ssw[acc1] = {}
#         for acc2 in best_edit_distances[acc1]:
#             seq2 = targets[acc2]
#             result = aligner.align(seq1, seq2, revcomp=False)
#             seq2_aln, match_line, seq1_aln = result.alignment
#             matches, mismatches, indels = match_line.count("|"), match_line.count("*"), match_line.count(" ")
#             insertion_count = seq2_aln.count("-")
#             deletion_count = seq1_aln.count("-")

#             sw_ed = mismatches + indels
#             best_edit_distances_ssw[acc1][acc2] =  sw_ed # (deletion_count, insertion_count, mismatches )
#             seq1_aln, match_line, seq2_aln = result.alignment
#             best_cigars_ssw[acc1][acc2] = (result.cigar, mismatches, indels, result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1 )

#             # print(acc1,acc2)
#             # print(result.query_begin, len(seq1) - result.query_end - 1, result.reference_begin, len(seq2) - result.reference_end -1, result.cigar, mismatches, indels)
#             # print()


#             # print(sw_ed, (deletion_count, insertion_count, mismatches ))
#             # print(seq1_aln)
#             # print(match_line)
#             # print(seq2_aln)
#             # edit_distance, locations, cigar = edlib_traceback(seq1, seq2, k =1000)
#             # print(edit_distance, locations, cigar)
#             # print()            
#     # for acc in best_cigars_ssw:
#     #     if len(best_cigars_ssw[acc]) ==0:
#     #         print("!!!!", acc)
#     # print(len(best_cigars_ssw))
#     # sys.exit()
#     return best_edit_distances_ssw, best_cigars_ssw



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
    print("Total number of perfect matches to reference database:", len(list(table.find(perfect_match_database="yes"))))
    print("Shared between samples and perfect match to reference database::", len(list(table.find(perfect_match_database="yes", both_samples = "yes"))))
    print("Shared between samples, full illumina support in both, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "yes", full_supp_sample2 = "yes" ))))
    print("Shared between samples, full illumina support in sample1 only, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "yes", full_supp_sample2 = "no" ))))
    print("Shared between samples, full illumina support in sample2 only, and perfect match to reference database:", len(list(table.find(perfect_match_database="yes", both_samples = "yes", full_supp_sample1 = "no", full_supp_sample2 = "yes" ))))
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
            is_fully_supported = (record["full_supp_sample1"] == record["full_supp_sample2"] == "yes")

        elif record["sample1"] == "yes":
            group = "sample1_specific"
            is_fully_supported = (record["full_supp_sample1"] == "yes")

        elif record["sample2"] == "yes":
            group = "sample2_specific"
            is_fully_supported = (record["full_supp_sample2"] == "yes")
        else:
            print("BUG, no category!")
            sys.exit()

        if is_fully_supported:
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

    # print to file
    print_database_to_tsv(db)

    # print fasta
    print_database_to_separate_fastas(db)


    # result = db['predicted_transcripts'].all()
    # dataset.freeze(result, mode='item', format='csv', filename=args.outfile)

    # print_shared_with_illumina_support(args.illumina_support_files, shared_seqs, seq_to_acc_sample1, seq_to_acc_sample2, args.outfolder)




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
    
    load.add_argument('--database', type=str, help='A fasta file with transcripts that are shared between smaples and have perfect illumina support.')

    args = parser.parse_args()
    path_, file_prefix = os.path.split(args.outfile)

    if path_:
        mkdir_p(path_)
    args.outfolder = path_

    main(args)

