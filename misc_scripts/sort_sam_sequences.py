import argparse
import os
import errno
import pysam

def main(args):

    #read SAM
    samfile = pysam.AlignmentFile(args.sam, "r")
    sam_records = {}
    for read in samfile.fetch():
        if read.query_name == "super_consensus":
            seq_id = read.query_name
            # print("here")
        else:
            seq_id = read.query_name.split("_")[1]
        sam_records[seq_id] = read


    cluster_id = 1
    clusters = { 1: set()}
    transcript_to_cluster = {}
    for line in open(args.seq_order, "r"):
        seq_id = line.strip()
        if seq_id == "-":
            cluster_id += 1
            clusters[cluster_id] = set()
        else:
            clusters[cluster_id].add(seq_id)
            transcript_to_cluster[seq_id] = cluster_id
    # print(sam_records)
    for cluster_id, cluster in clusters.items():
        #print in order according to seq_order file
        cluster_sam_file = pysam.AlignmentFile(os.path.join(args.outfolder, str(cluster_id) + ".sam" ), "w", template=samfile)
        read = sam_records["super_consensus"]
        cluster_sam_file.write(read)
        for seq_id in cluster:
            read = sam_records[seq_id]
            cluster_sam_file.write(read)

    # for seq_id, read in sam_records.items():

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

# Take care of input

    parser = argparse.ArgumentParser(description = "This script does: \n\n {0}".format(main.__doc__)) 
    parser.add_argument('--seq_order', type=str, help='Tsv file ')
    parser.add_argument('--sam', type=str, help='SAM file ')
    parser.add_argument('--outfolder', type=str, help='outfolder path. ')

    args = parser.parse_args()

    mkdir_p(args.outfolder)
    main(args)