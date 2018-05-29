import argparse
import os
import random
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import seaborn as sns
import pandas as pd

from collections import defaultdict
def transpose(dct):
    d = defaultdict(dict)
    for key1, inner in dct.items():
        for key2, value in inner.items():
            d[key2][key1] = value
    return d  


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

def swarmplot(args):
    plt.clf()
    with sns.plotting_context("paper"): #, font_scale=1.0):
        indata = pd.read_csv(args.tsv_input, sep="\t") #, dtype = { "alignment_id": float, "read_support" : int})
        # indata["alignment_id"] = pd.to_numeric(indata["alignment_id"],  errors='ignore' ) #, downcast='float')
        # indata["read_support"] = pd.to_numeric(indata["read_support"],  errors='ignore') #, downcast='signed')
        print(indata)
        ax = sns.swarmplot(x="read_support", y="alignment_id", hue="sample", order = ["1-5", "6-10", "11-20", "21-50", ">50"], hue_order = ["shared", "sample1_specific", "sample2_specific"], data=indata)
        # ax = sns.heatmap(indata, annot=True, mask = mask, fmt="d", cmap="hot_r", vmin=1, vmax=500 )
        plt.legend(loc='lower right')
        plt.savefig(args.outfile)
        plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('tsv_input', type=str, help='tsv file with hits.')
    parser.add_argument('outfile', type=str, help='Output filename')
    args = parser.parse_args()

    swarmplot(args)