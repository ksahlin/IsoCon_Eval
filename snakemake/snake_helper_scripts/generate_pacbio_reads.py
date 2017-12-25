import os
import argparse
import numpy as np
import misc_functions
import random
from collections import defaultdict

# log how many reads are simulated for each variant.


# def 

def main(params):
    # config = misc_functions.read_config(params.config)
    sequence_transcripts = {}
    # for acc, seq in misc_functions.read_fasta(open(params.sequence_material, "r")):
    #     sequence_transcripts[acc] = seq  
    sequence_transcripts = dict(misc_functions.read_fasta(open(params.sequence_material,"r")))
    # print(sequence_transcripts)
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
    pacbio_reads = {}
    reads_generated_log = defaultdict(int)
    errors = []
    while read_count <= params.read_count:
        if it >= len(lengths):
            lengths = np.random.triangular(0, 10000, 45000, 5*params.read_count)
            it = 0

        read_len = lengths[it]
        acc = random.choice( list(sequence_transcripts.keys()))
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
            # params.logfile.write("{0}, error places: {1}\n".format(read_acc, error_log))
            reads_generated_log[acc.split(":copy")[0]] += 1
            errors.append(total_error_length)

            pacbio_reads[read_acc] = read
            read_count += 1

        it += 1


    for acc, abundance in misc_functions.iteritems(reads_generated_log):
        params.logfile.write("{0}\t{1}\n".format(acc, abundance))

    n = float(len(errors))
    mu =  sum(errors) / n
    sigma = (sum(list(map((lambda x: x ** 2 - 2 * x * mu + mu ** 2), errors))) / (n - 1)) ** 0.5
    min_error = min(errors)
    max_error = max(errors)
    errors.sort()
    if len(errors) %2 == 0:
        median_error = (errors[int(len(errors)/2)-1] + errors[int(len(errors)/2)]) / 2.0
    else:
        median_error = errors[int(len(errors)/2)]

    params.logfile.write("mean error: {0}, sd error:{1}, min_error:{2}, max_error:{3}, median_error:{4}\n".format(mu, sigma, min_error, max_error, median_error))

    out_file = open(params.outfile, "w")
    for acc, seq in misc_functions.iteritems(pacbio_reads):
        out_file.write(">{0}\n{1}\n".format(acc,seq))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate pacbio reads from a set of transcripts.")
    parser.add_argument('sequence_material', type=str, help='The fasta file with sequences to be sequenced.')
    parser.add_argument('outfile', type=str, help='Output path to fasta file')
    parser.add_argument('read_count', type=int, help='Number of reads to simulate.')
    # parser.add_argument('config', type=str, help='config file')


    params = parser.parse_args()
    path_, file_prefix = os.path.split(params.outfile)
    misc_functions.mkdir_p(path_)
    params.logfile = open(os.path.join(path_, file_prefix + ".log"), "w")
    main(params)
