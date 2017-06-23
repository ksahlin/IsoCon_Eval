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

def plot_binary_membership(binary_membership_file, args):
    sns.set_style("whitegrid")
    dataset = pd.read_csv(binary_membership_file, sep="\t")
    ax = sns.countplot(x="GENE_FAMILY", hue="METHOD", data=dataset)
    plt.xlabel("Family")
    plt.ylabel("# Best Hits")
    plt.title("Unique transcript ID's found in consensus of each method")
    outfile = os.path.join(args.outfolder, "binary_memebership.pdf")
    plt.savefig(outfile)

def get_best_hits_over_identity_threshold(file_name, targeted, args):

    # wrong in this function!!?
    # target can have several hits! we shoul also check for reduncdance and %identity
    # reads can have same identity to several identical targets!!

    best_hits = {}
    pattern = re.compile('BPY|CDY|DAZ|HSFY|PRY|RBMY|TSPY|XKRY|VCY')
    queries_seen = defaultdict(list)
    for line in open(file_name):
        # print(line)
        # print(line.strip().split("\t"))
        query, target, len_query, len_target, ed = line.strip().split("\t")
        ed = int(ed)
        # identity_ = 1.0 - ( (ed - 2*21)/float(max(len_query, len_target)) )

        if pattern.search(target):
            # if args.min_percentage > identity_:
            #     continue
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

    binary_membership_outfile = open(os.path.join(args.outfolder, "hit_to_db.tsv"), "w")
    binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format("ID", "METHOD","GENE_FAMILY", "ED"))
    pattern = re.compile('BPY|CDY|DAZ|HSFY|PRY|RBMY|TSPY|XKRY|VCY')
    for target in flnc_hits:
        ed = flnc_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "FLNC", family, ed))
    for target in isocon_hits:
        ed = flnc_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "ISOCON", family, ed))
    for target in ice_hits:
        ed = flnc_hits[target]
        family = pattern.search(target).group(0)
        binary_membership_outfile.write("{0}\t{1}\t{2}\t{3}\n".format(target, "ICE", family, ed))

    binary_membership_outfile.close()
    print(len(flnc_hits), flnc_hits)
    print(len(isocon_hits), isocon_hits)
    print(len(ice_hits), ice_hits)
    flnc_best = 0
    isocon_best = 0
    ice_best = 0
    for target in flnc_hits:
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
    print("TOTAL BEST:", "FLNC:",flnc_best, "IsoCon:",isocon_best, "ICE:", ice_best)
    plot_binary_membership(binary_membership_outfile.name, args)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--flnc', type=str, help='Path to the predicted transcript fasta file')
    parser.add_argument('--isocon', type=str, help='Path to the predicted transcript fasta file')
    parser.add_argument('--ice', type=str, help='Path to the predicted transcript fasta file')
    parser.add_argument('--min_percentage', type=float, default = 0.95, help='Minimum identity threshold to be considered')
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    

    main(args)


