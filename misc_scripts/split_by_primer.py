
import os
import argparse
import re


def reverse_complement(string):
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

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
            temp += line.strip()
    yield accession, temp


def main(args):
    reads = {acc: seq for (acc, seq) in  read_fasta(open(args.reads, 'r'))}
    primers = {acc: seq for (acc, seq) in  read_fasta(open(args.primer_file, 'r'))}
    print(primers)
    outfiles_dict = {}
    path_, file_prefix = os.path.split(args.reads)

    for acc in primers:
        m = re.search("[0-9]+", acc)
        if m:
        # print("lol", m.group(0))
            primer_outfolder = os.path.join(args.outfolder, m.group(0))
            if not os.path.exists(primer_outfolder):
                os.makedirs(primer_outfolder)
            # print(primer_outfolder)
            p = m.group(0)
            # print(p)
            outfiles_dict[p] = open(os.path.join(primer_outfolder, file_prefix ), "w")

    print(outfiles_dict.keys())
    # log primer sequence plus if exact primer was fount in accession of sequences
    for acc, seq in reads.items():
        # print(acc)
        m = re.search("primer=[0-9]*", acc)
        if m.group(0).split("=")[1]:
            primer = m.group(0).split("=")[1]
            outfiles_dict[primer].write(">{0}\n{3}\n".format(acc, seq))

            # print(primer, acc)
            # if primers["F" + primer] in seq and primers["R" + primer] in seq:
            #     tag = "both_exact"
            #     outfiles_dict[primer].write(">{0}_{1}_{2}\n{3}\n".format(acc, tag, primer, seq))
            #     print("B")
            # elif primers["F" + primer] in seq:
            #     tag = "F_exact"
            #     outfiles_dict[primer].write(">{0}_{1}_{2}\n{3}\n".format(acc, tag, primer, seq))
            #     print("F")

            # elif primers["R" + primer] in seq:
            #     tag = "R_exact"
            #     outfiles_dict[primer].write(">{0}_{1}_{2}\n{3}\n".format(acc, tag, primer, seq))
            #     print("R")

            # else:
            #     tag = "None_exact"
            #     outfiles_dict[primer].write(">{0}_{1}_{2}\n{3}\n".format(acc, tag, primer, seq))


if __name__ == '__main__':


    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('--reads', type=str, help='Fasta file ')
    parser.add_argument('--primer_file', type=str, help='Fasta file ')
    parser.add_argument('--outfolder', type=str, help='outfile folder to put output in. ')

    args = parser.parse_args()

    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)
