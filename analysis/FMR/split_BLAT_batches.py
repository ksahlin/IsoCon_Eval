
import argparse
import sys
import os

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

def main(args):

    predictions =  {acc.replace("|", "_") : seq.upper() for acc, seq in read_fasta(open(args.predictions,"r"))} 

    batch_size = 25
    batch_number = 0
    # outfile = open(os.path.join(args.outfolder, args.prefix + _"{0}".format(batch_number)))
    for i, (acc, seq) in enumerate(predictions.items()):
        if i % batch_size == 0:
            batch_number += 1
            outfile = open(os.path.join(args.outfolder, args.prefix + "_{0}.fa".format(batch_number)), "w")
            outfile.write(">{0}\n{1}\n".format(acc, seq))            
        else:
            outfile.write(">{0}\n{1}\n".format(acc, seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('predictions', type=str, help='predictions in fasta format.')
    parser.add_argument('outfolder', type=str, help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')    
    parser.add_argument('prefix', type=str, help='Outfile prefix.')
    # parser.add_argument('--realign', action="store_true", help='A fasta file with transcripts that are shared between samples and have perfect illumina support.')
    
    args = parser.parse_args()


    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    if args.outfolder and not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)


    main(args)
