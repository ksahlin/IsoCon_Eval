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

def plot_binary_membership(binary_membership_file, args):
    sns.set_style("whitegrid")
    dataset = pd.read_csv(binary_membership_file, sep="\t")
    ax = sns.countplot(x="GENE_FAMILY", hue="METHOD", data=dataset)
    plt.xlabel("Family")
    plt.ylabel("# Best Hits")
    plt.title("Unique transcript ID's found in consensus of each method")
    outfile = os.path.join(args.outfolder, "binary_memebership.pdf")
    plt.savefig(outfile)

def get_best_hits(file_name, targeted):
    best_hits = {}
    pattern = re.compile('BPY|CDY|DAZ|HSFY|PRY|RBMY|TSPY|XKRY|VCY')
    for line in open(file_name):
        # print(line)
        # print(line.strip().split("\t"))
        query, target, len_query, len_target, ed = line.strip().split("\t")
        ed = int(ed)

        if pattern.search(target):
            if target in best_hits:
                if ed < best_hits[target]:
                    best_hits[target] = ed

            else:
                best_hits[target] = ed

    return best_hits

def main(args):
    targeted = set(["BPY", "CDY", "DAZ", "HSFY", "PRY", "RBMY", "TSPY", "XKRY", "VCY"])
    # for plotting simple binary membership    
    flnc_hits = get_best_hits(args.flnc, targeted)
    isocon_hits = get_best_hits(args.isocon, targeted)
    ice_hits = get_best_hits(args.ice, targeted)

    binary_membership_outfile = open(os.path.join(args.outfolder, "hit_to_db.tsv"), "w")
    binary_membership_outfile.write("{0}\t{1}\t{2}\n".format("ID", "METHOD","GENE_FAMILY"))
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
    print(len(flnc_hits))
    print(len(isocon_hits))
    print(len(ice_hits))

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
    parser.add_argument('--outfolder', type=str, help='Output path of results')
    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)
    

    main(args)


