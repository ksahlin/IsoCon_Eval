from __future__ import print_function
import sys
import argparse
import re
import itertools
import os
import argparse
import numpy as np
import misc_functions
import random
from collections import defaultdict
# import math

def mutate(sequence, positions):
    new_sequence = list(sequence)
    choices = [("m", 0.3333), ("i", 0.3333), ("d", 0.33333)]
    for p in positions:
        r = random.uniform(0,1)
        if r < 0.3333:
            new_sequence[p] = random.choice( list(set(["A","G","C","T"]) - set(new_sequence[p])) )
            assert new_sequence[p] != sequence[p]
        elif 0.3333 <= r < 0.6666:
            new_sequence[p] = new_sequence[p] + random.choice(["A","G","C","T"])
            assert new_sequence[p] != sequence[p]
        else:
            new_sequence[p] = ""

    return "".join([n for n in new_sequence])


def generate_reads(ref_seq, mutated_ref_seq, params):
    sequence_transcripts = {"ref_seq": ref_seq, "mutated_ref_seq" : mutated_ref_seq}
    # read lengths ~ according to P6-C4 chemistry histogram from here 
    # http://www.slideshare.net/GenomeInABottle/jan2016-pac-bio-giab   slide 13
    # this looks like it can be well approximated by triangiular distributions with parameters
    # 0 (base start), 10000 (peak), ~45000 (base end)
    # http://docs.scipy.org/doc/numpy/reference/generated/numpy.random.triangular.html
    # pacbios own distribution is here:
    # http://www.pacb.com/blog/new-chemistry-boosts-average-read/
    
    # read_lengths = np.random.triangular(0, 10000, 45000, config["read_count"])

    # Get average quality based on subread length and the length of the transcript
    # while read count is less than the red count we want:
    #   1. Draw read length from distribution for each read length in simulated read lengths
    #   2. Randomly select a transcript from pool
    #   3. Get average quality based on the number of passes = floor(read_length/transcript_length)
    #       Avg quality is derived from this plot: https://speakerdeck.com/pacbio/specifics-of-smrt-sequencing-data
    #          slide 21, the P4 C2 chemistry, in case of one pass we chose 13 percent error rate from here : http://www.sciencedirect.com/science/article/pii/S1672022915001345.
    #          Out of the errors we follow the data here: http://bib.oxfordjournals.org/content/17/1/154.full.pdf
    #           and here http://www.homolog.us/Tutorials/index.php?p=2.8&s=2
    #           that indicates that for a pacbio genomic read, we have roughly 13-15 percent error rate (older chemistry)
    #           we choose 13. Out of the total errors, we let 11/16 = 68.75 be insertions
    #           4/16= 25% be deletions and 1/16 = 6.25% be substitutions  (given here http://www.homolog.us/Tutorials/index.php?p=2.8&s=2 and http://bib.oxfordjournals.org/content/17/1/154.full.pdf)
    #   4. generate the read

    quality_function = {1: 0.87, 2:0.95, 3: 0.957, 4:0.969, 5:0.981, 6:0.985, 7:0.99, 
                8:0.992, 9:0.994, 10: 0.995,11: 0.995,12: 0.995, 13: 0.996, 14: 0.996, 15: 0.996,
                16: 0.999, 17: 0.999, 18: 0.999} # passes : quality
    read_count = 1
    # just generate all numbers at once and draw from this 5x should be enough
    it = 0
    lengths = np.random.triangular(0, 10000, 45000, 5*params.read_count)
    ref_choice = np.random.choice(["ref_seq", "mutated_ref_seq" ], size = 5*params.read_count, p = [1.0-params.abundance_ratio, params.abundance_ratio])
    pacbio_reads = {}
    # reads_generated_log = defaultdict(int)
    # errors = []
    while read_count <= params.read_count:
        if it >= len(lengths):
            lengths = np.random.triangular(0, 10000, 45000, 5*params.read_count)
            it = 0

        read_len = lengths[it]
        acc = ref_choice[it]
        transcript = sequence_transcripts[acc]
        passes =  int(read_len/ len(transcript))
        # print(passes, read_len, len(transcript))
        if passes > 0:
            if passes < 18:
                quality = quality_function[passes]
            else:
                quality = 0.999

            subs_rate = (1.0 - quality)*0.0625
            ins_rate = (1.0 - quality)*0.6875
            del_rate = (1.0 - quality)*0.25
            read, error_log, total_error_length, total_indel_length, total_del_length = misc_functions.mutate_sequence(transcript, subs_rate, ins_rate, del_rate)
            read_acc = "{0}_read_{1}_error_rate_{2}_total_errors_{3}".format(acc, str(read_count), total_error_length/float(len(read) + total_del_length), total_error_length)
            # reads_generated_log[acc.split(":copy")[0]] += 1
            # errors.append(total_error_length)

            pacbio_reads[read_acc] = read
            read_count += 1

        it += 1


    # for acc, abundance in misc_functions.iteritems(reads_generated_log):
    #     params.logfile.write("{0}\t{1}\n".format(acc, abundance))

    return pacbio_reads



def main(params):
    # config = misc_functions.read_config(params.config)

    reference = {acc : seq for acc, seq in misc_functions.read_fasta(open(params.transcript,"r"))} 
    ref_seq = list(reference.values())[0]
    while True:
        random_positions = set([random.randint(0, len(ref_seq) - 1) for i in range(params.ed)] )
        if len(random_positions) == params.ed:
            break

    mutated_ref_seq = mutate(ref_seq, random_positions)
    out_file_ref = open(params.ref_outfile, "w")
    out_file_ref.write(">{0}\n{1}".format("mutated_ref_seq", mutated_ref_seq))
    out_file_ref.close()

    # simulate reads
    reads =  generate_reads(ref_seq, mutated_ref_seq, params)
    out_file_reads = open(params.reads_outfile, "w")
    for acc, seq in misc_functions.iteritems(reads):
        out_file_reads.write(">{0}\n{1}\n".format(acc,seq))
    out_file_reads.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate a gene family and isoforms from a set of original exons.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--transcript', type=str, help='The fasta file with the full transcript.')
    parser.add_argument('--ed', type=int, help='Mutation rate')
    parser.add_argument('--abundance_ratio', type=float, help='Mutation rate')
    parser.add_argument('--read_count', type=int, help='Number of reads to simulate.')
    parser.add_argument('reads_outfile', type=str, help='Generated reads output file')
    parser.add_argument('ref_outfile', type=str, help='generated ref output file ')

    params = parser.parse_args()
    path_, file_prefix = os.path.split(params.reads_outfile)
    misc_functions.mkdir_p(path_)
    main(params)

