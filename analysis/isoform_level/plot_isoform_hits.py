import argparse
import os, sys
import re
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import numpy as np
import seaborn as sns
import pandas as pd
from collections import defaultdict
import errno

def plot_binary_membership(binary_membership_file, args):
    # sns.set_style("whitegrid")
    sns.set_color_codes("muted")
    dataset = pd.read_csv(binary_membership_file, sep="\t")
    ax = sns.countplot(x="GENE_FAMILY", hue="METHOD", data=dataset, hue_order=["ISOCON", "ICE", "FLNC"], palette={"ISOCON": "b", "ICE": "g", "FLNC" : "r"})
    plt.xlabel("Family")
    plt.ylabel("# Perfect matches to distinct transcripts")
    plt.title("Perfect matches to transcripts in ENSEMBL")
    # outfile = os.path.join(args.outfolder, "binary_memebership.pdf")
    plt.savefig(args.outprefix)

def get_best_hits_over_identity_threshold(file_names, targeted, args):

    # wrong in this function!!?
    # target can have several hits! we shoul also check for reduncdance and %identity
    # reads can have same identity to several identical targets!!

    best_hits = {}
    pattern = re.compile('BPY|CDY|DAZ|HSFY|PRY|RBMY|TSPY|XKRY|VCY')
    queries_seen = defaultdict(list)
    for file_ in file_names:
        for line in open(file_):
            # print(line)
            # print(line.strip().split("\t"))
            query, target, len_query, len_target, ed = line.strip().split("\t")
            ed = int(ed)
            # identity_ = 1.0 - ( (ed - 2*21)/float(max(len_query, len_target)) )

            if pattern.search(target):
                if ed > args.max_ed:
                    continue
                if query in queries_seen:
                    print("already seen, mapping to:", queries_seen[query], target, query )
                queries_seen[query].append(target)


                if target in best_hits:
                    if ed < best_hits[target]:
                        best_hits[target] = ed

                else:
                    best_hits[target] = ed


    return best_hits

def main(args):
    targeted = set(["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    # for plotting simple binary membership    
    flnc_hits = get_best_hits_over_identity_threshold(args.flnc, targeted, args)
    isocon_hits = get_best_hits_over_identity_threshold(args.isocon, targeted, args)
    ice_hits = get_best_hits_over_identity_threshold(args.ice, targeted, args)

    # do temporary file instead of creating folder here...
    binary_membership_outfile = open(args.outprefix +"_hit_to_db.tsv", "w")
    binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format("ID", "METHOD","GENE_FAMILY", "ED"))
    pattern = re.compile('BPY|CDY|DAZ|HSFY|PRY|RBMY|TSPY|XKRY|VCY')
    for target in flnc_hits:
        ed = flnc_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "FLNC", family, ed))
    for target in isocon_hits:
        ed = isocon_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "ISOCON", family, ed))
    for target in ice_hits:
        ed = ice_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "ICE", family, ed))

    binary_membership_outfile.close()
    print(len(flnc_hits), flnc_hits)
    print(len(isocon_hits), isocon_hits)
    print(len(ice_hits), ice_hits)
    flnc_best = 0
    isocon_best = 0
    ice_best = 0
    already_processed = set()
    for target in flnc_hits:
        already_processed.add(target)
        ed1 = flnc_hits[target]
        if target in isocon_hits:
            ed2 = isocon_hits[target]
        else:   
            ed2 = 2**32
        if target in ice_hits:
            ed3 = ice_hits[target]
        else:   
            ed3 = 2**32

        min_ed = min(ed1,ed2,ed3)
        if ed1 == min_ed:
            flnc_best += 1
        if ed2 == min_ed:
            isocon_best += 1
        if ed3 == min_ed:
            ice_best += 1
        print("FLNC:",ed1, "IsoCon:",ed2, "ICE:", ed3, target)

    print("HERE")

    for target in ice_hits:
        if target in already_processed:
            continue
        else:
            already_processed.add(target)
            ed1 = ice_hits[target]
            if target in isocon_hits:
                ed2 = isocon_hits[target]
            else:   
                ed2 = 2**32
            if target in flnc_hits:
                ed3 = flnc_hits[target]
            else:   
                ed3 = 2**32

            min_ed = min(ed1,ed2,ed3)
            if ed1 == min_ed:
                ice_best += 1
            if ed2 == min_ed:
                isocon_best += 1
            if ed3 == min_ed:
                flnc_best += 1
            print("FLNC:",ed3, "IsoCon:",ed2, "ICE:", ed1, target)

    print("HERE")

    for target in isocon_hits:
        if target in already_processed:
            continue
        else:
            already_processed.add(target)
            ed1 = isocon_hits[target]
            if target in ice_hits:
                ed2 = ice_hits[target]
            else:   
                ed2 = 2**32
            if target in flnc_hits:
                ed3 = flnc_hits[target]
            else:   
                ed3 = 2**32

            min_ed = min(ed1,ed2,ed3)
            if ed1 == min_ed:
                isocon_best += 1
            if ed2 == min_ed:
                ice_best += 1
            if ed3 == min_ed:
                flnc_best += 1
            print("FLNC:",ed3, "IsoCon:",ed1, "ICE:", ed2, target)

    print("TOTAL BEST:", "FLNC:",flnc_best, "IsoCon:",isocon_best, "ICE:", ice_best)
    plot_binary_membership(binary_membership_outfile.name, args)
    # os.remove(binary_membership_outfile.name)
    # FN = {}
    # for db_hit in read_hits:
    #     if db_hit not in predicted_hits:
    #         FN[db_hit] = read_hits[db_hit]
    #         print("not found:", db_hit, read_hits[db_hit])
    #     elif predicted_hits[db_hit] > read_hits[db_hit]:
    #         print("read had better hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
    #         FN[db_hit] = read_hits[db_hit]
    #     elif predicted_hits[db_hit] == read_hits[db_hit]:
    #         print("Same quality hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
    #     else:
    #         # print("predicted had better hit", predicted_hits[db_hit], read_hits[db_hit])
    #         continue

    # print("False negatives")
    # # for fn in FN:
    # #     print(FN[fn], fn)
    # better_hit_to_pred = {}
    # for db_hit in predicted_hits:
    #     if db_hit not in read_hits:
    #         better_hit_to_pred[db_hit] = predicted_hits[db_hit]
    #         print("not found in reads:", db_hit, better_hit_to_pred[db_hit])
    #     elif predicted_hits[db_hit] < read_hits[db_hit]:
    #         print("Predicted had better hit", predicted_hits[db_hit], read_hits[db_hit], db_hit)
    #         better_hit_to_pred[db_hit] = predicted_hits[db_hit]
    #     else:
    #         # print("read had better hit", predicted_hits[db_hit], read_hits[db_hit])
    #         continue

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
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--flnc', type=str, nargs="+", help='Path to the tsv file of best hits to database')
    parser.add_argument('--isocon', type=str, nargs="+", help='Path to the tsv file of best hits to database')
    parser.add_argument('--ice', type=str, nargs="+", help='Path to the tsv file of best hits to database')
    parser.add_argument('--max_ed', type=int, default = 0, help='Maximum local edit distance to reference in order to be considered a hit [default 0, consited only perfect matches].')
    parser.add_argument('--outprefix', type=str, help='Output path of results')
    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outprefix)
    mkdir_p(path_)
    args.outfolder = path_

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)
    
    main(args)


