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

def heatmap(args):
    ref_hits =  {acc : seq for acc, seq in read_fasta(open(args.references,"r"))} 

    data_dict = {r : {"1009": 0, "1015": 0, "1018": 0, "4549": 0, "5123": 0, "5248": 0 } for r in ref_hits}
    for i, line in enumerate(open(args.tsv_input, "r")):
        # if i == 0:
        #     continue
        # else:
        q, r, sample = line.split()
        read_support = int(q.split("_")[4])
        data_dict[r][sample] += read_support

    mask = {r : {"1009": False, "1015": False, "1018": False, "4549": False, "5123": False, "5248": False } for r in ref_hits}
    for r in ref_hits:
        for s in ["1009", "1015", "1018", "4549", "5123", "5248"]:
            if data_dict[r][s] == 0:
                mask[r][s] = True
            # else:
            #     mask[r][s] = False

    ################
    tmp_outfile = open("/Users/kxs624/tmp/isocon_fmr.tsv", "w")
    for r in sorted(ref_hits):
        row = str(r) + "\t"
        for s in ["1018", "1015", "1009", "4549", "5123", "5248"]:
            row += str(data_dict[r][s]) + "\t"
        row += "\n"
        tmp_outfile.write(row)
    ################

    data_dict = transpose(data_dict)
    mask = transpose(mask)

    plt.clf()
    with sns.plotting_context("paper"): #, font_scale=1.0):
        indata = pd.DataFrame(data_dict)
        mask = pd.DataFrame(mask)
        # indata = pd.read_csv(args.tsv_input, sep="\t")
        # indata = indata.pivot("Tseng2017", "sample", "IsoCon")
        # grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
        # f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
        ax = sns.heatmap(indata, annot=True, mask = mask, fmt="d", cmap="hot_r", vmin=1, vmax=500 )
        plt.savefig(args.outfile)
        plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('tsv_input', type=str, help='tsv file with hits.')
    parser.add_argument('references', type=str, help='The fasta file references.')
    parser.add_argument('outfile', type=str, help='Output filename')
    args = parser.parse_args()

    heatmap(args)