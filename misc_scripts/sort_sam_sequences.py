import argparse
import os
import pysam

def main(args):

    #read SAM
    samfile = pysam.AlignmentFile(args.sam, "r")
    sam_records = {}
    for read in samfile.fetch():
        seq_id = read.query_name.split("_")[1]
        sam_records[seq_id] = read

    #print in order according to seq_order file
    sorted_file = pysam.AlignmentFile(args.outfile, "w", template=samfile)
    for line in open(args.seq_order, "r"):
        seq_id = line.strip()
        read = sam_records[seq_id]
        sorted_file.write(read)

    # for seq_id, read in sam_records.items():


if __name__ == '__main__':

# Take care of input

    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('--seq_order', type=str, help='Tsv file ')
    parser.add_argument('--sam', type=str, help='SAM file ')
    parser.add_argument('--outfile', type=str, help='outfile path. ')

    args = parser.parse_args()


    main(args)