
import sys,os,argparse

try:
    # import matplotlib
    # matplotlib.use('agg')
    # import matplotlib.pyplot as plt
    import pylab as plt
    import seaborn as sns
    sns.set_palette("husl", desat=.6)
    sns.set(font_scale=1.6)
    plt.rcParams.update({'font.size': 22})
except:
    pass

def histogram(data, args, name='histogram.png', x='x-axis', y='y-axis', x_cutoff=None, title=None):
    plt.xscale('log')
    fig, ax = plt.subplots()
    import numpy as np
    MIN, MAX = 1.0e-20 , 1.0 #max(data)
    print(sorted(data))
    data= [max(1.1e-20, p) for p in data]
    print(sorted(data))

    plt.hist(data, bins=np.logspace(np.log10(MIN),np.log10(MAX), 50))
    plt.gca().set_xscale("log")


    # if x_cutoff: 
    #     plt.hist(data, range=[0, x_cutoff], bins = 100)
    # else:
    #     plt.hist(data, bins = 100)
    # plt.xlabel(x)
    # plt.ylabel(y)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(title)
    fig.tight_layout()
    # if title:
    #     plt.title(title)

    plt.savefig(os.path.join(args.outfolder, args.prefix))
    plt.clf()

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


def main(args):
    if args.fasta:
        candidate_dict = {acc : seq for acc, seq in read_fasta(open(args.fasta,"r"))} 
    elif args.tsv:
        candidate_dict = {line.split()[0] : line.split()[1] for line in open(args.tsv,"r")} 

    p_values = []
    for acc in candidate_dict:
        p_val = acc.split("_")[5]
        if p_val == "not":
            pass
            # p_values.append(0)
        else:
            p_values.append( float( acc.split("_")[5]) )

    # p_values = [float( acc.split("_")[5]) for acc in candidate_dict]
    histogram(p_values, args, title=args.prefix, x="p-value", y="Count" )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('outfolder', type=str, help='outfolder.')  
    parser.add_argument('prefix', type=str, help='prefix to outfile.') 
    parser.add_argument('--fasta', type=str, help='fasta.')
    parser.add_argument('--tsv', type=str, help='tsv.') 


    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)
