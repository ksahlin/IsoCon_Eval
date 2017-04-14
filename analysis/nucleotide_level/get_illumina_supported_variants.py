


import pysam
samfile = pysam.AlignmentFile("ex1.bam", "rb" )
for pileupcolumn in samfile.pileup("chr1", 100, 120):
    print ("\ncoverage at base %s = %s" %
           (pileupcolumn.pos, pileupcolumn.n))
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            # query position is None if is_del or is_refskip is set.
            print ('\tbase in read %s = %s' %
                  (pileupread.alignment.query_name,
                   pileupread.alignment.query_sequence[pileupread.query_position]))

samfile.close()


import argparse
import os
from collections import Counter
import re

def main(args):
    output_file = open(os.path.join(args.outfolder, "illumina_consensus.tsv" ), "w")

    FP_positions = []
    for line in open(args.support_file, "r"):
        consensus_name, ref_pos, ref_site, total_support, read_bases, base_qualities =  line.strip().split()

        # # here we discover unsopported positions (FP)
        # if int(total_support) < 2:
        #     FP_positions.append(int(ref_pos))

        # get all alternate loci and their support over positions

        ins_pattern = "\+[0-9]+[ACGTNacgtn]+"
        ins_res = re.findall(ins_pattern, read_bases.upper())
        c_ins = Counter(ins_res)


        del_pattern =  "-[0-9]+[ACGTNacgtn]+"
        del_res = re.findall(del_pattern, read_bases.upper())
        c_del = Counter(del_res)


        subs_pattern = "(?<![0-9'^'])[ACGTN]"
        subs_res = re.findall(subs_pattern, read_bases.upper())
        subs_res_string = "".join([subs for subs in subs_res])
        c_subs = Counter(subs_res)
        

        c = Counter(read_bases.upper())
        inferred_site_count = c["."] + c[","] + len(ins_res)  + len(subs_res) 
        print(total_support, inferred_site_count, c["."], c[","] , len(ins_res) , len(del_res) , len(subs_res)  )
        assert int(total_support) == inferred_site_count

        output_file.write("{0}\t{1}\t{2}:{3}\t".format(consensus_name, ref_pos, ref_site, c["."] + c[","], ))

        for site, count in  c_subs.items():
            output_file.write("{0}:{1}\t".format(site,count))

        for site, count in  c_ins.items():
            output_file.write("{0}:{1}\t".format(site,count))

        for site, count in  c_del.items():
            output_file.write("{0}:{1}\t".format(site,count))
        

if __name__ == '__main__':

# Take care of input


    parser = argparse.ArgumentParser(description = "Get illumina read support on consensus transcripts")
    parser.add_argument('support_file', type=str, help='A tsv file which is the output from samtools mpileup -a -f. ')
    parser.add_argument('outfolder', type=str, help='outfile folder to put output in. ')

    args = parser.parse_args()
    if not os.path.exists(args.outfolder):
        os.makedirs(args.outfolder)

    main(args)