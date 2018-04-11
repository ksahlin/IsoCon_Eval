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
        if i == 0:
            continue
        else:
            q, r, sample = line.split()
            data_dict[r][sample] = 1

    data_dict = transpose(data_dict)
    plt.clf()
    with sns.plotting_context("paper", font_scale=1.0):
        indata = pd.DataFrame(data_dict)
        # indata = pd.read_csv(args.tsv_input, sep="\t")
        # indata = indata.pivot("Tseng2017", "sample", "IsoCon")
        ax = sns.heatmap(indata, annot=True, fmt="d", vmin=1)
        plt.savefig(args.outfile)
        plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('tsv_input', type=str, help='tsv file with hits.')
    parser.add_argument('references', type=str, help='The fasta file references.')
    parser.add_argument('outfile', type=str, help='Output filename')
    args = parser.parse_args()

    heatmap(args)