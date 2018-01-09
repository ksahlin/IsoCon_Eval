import argparse
import os
import errno
from math import *
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")
# import matplotlib.pyplot as plt
# import matplotlib

import seaborn as sns
import pysam

def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def get_ccs(ccs_file):  
    passes = []
    for read in ccs_file.fetch(until_eof=True):
        # print(read.query_name)
        # ccs_id = int(read.query_name.strip().split("/")[-2])
        # ccs_read = CCS(ccs_id, read.query_alignment_sequence, read.query_qualities, read.get_tag("np"))
        passes.append(int(read.get_tag("np")))        
    return passes


def main(args):
    ccs_file = pysam.AlignmentFile(args.bamfile, "rb", check_sq=False)
    num_passes = get_ccs(ccs_file)

    print("mean:", sum(num_passes)/float(len(num_passes)), "median:", sorted(num_passes)[int(len(num_passes)/2)], "min:", min(num_passes), "max:", max(num_passes) )

    plt.hist(num_passes, 50)
    plt.ylabel('Count')
    plt.xlabel('Number of passes')
    plt.savefig(args.outfile)
    plt.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="plot nr passes of reads.")
    parser.add_argument('--bamfile', type=str, help='Path to ccs bam or sam file.')
    parser.add_argument('--outfile', type=str, help='Output path of plot')
    args = parser.parse_args()

    path_, file_prefix = os.path.split(args.outfile)
    mkdir_p(path_)
    args.outfolder = path_

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)
    
    main(args)