#!/usr/bin/env bash
python --version
root="/Users/kxs624/tmp/FMR1/analysis/Pauls"
transcript_folder="/Users/kxs624/tmp/FMR1"
psl_folder="/Users/kxs624/tmp/FMR1/analysis/Pauls/blat_psl_files"
sam_folder="/Users/kxs624/tmp/FMR1/analysis/Pauls/sam_files"
fasta_batches="/Users/kxs624/tmp/FMR1/analysis/Pauls/fasta_batches"
output="/Users/kxs624/tmp/FMR1/analysis/Pauls/analysis_output"
fasta="/Users/kxs624/tmp/FMR1/"

samples="1009 1015 1018 4549 5123 5248"

# samples="5123 5248"  #"1009 1015 1018 4549"  #"53_isoforms 1009"

for sample in $samples
do
    # python split_BLAT_batches.py $transcript_folder/$sample/final_candidates.fa $fasta_batches/$sample
    rm $psl_folder/$sample/all.psl

    awk 'FNR==1{print ""}1' $psl_folder/$sample/*.psl > $psl_folder/$sample/all.psl

    uncle_psl.py -N 50 -f $transcript_folder/$sample/final_candidates.fa $psl_folder/$sample/all.psl $sam_folder/$sample\_blat.sam

    cat $root/sam_header_tabs.sam $sam_folder/$sample\_blat.sam > $sam_folder/$sample\_blat_w_header.sam

    echo python get_distinct_isoforms.py $sam_folder/53_isoforms_blat_w_header.sam $sam_folder/$sample\_blat_w_header.sam $output $sample --fasta $fasta/$sample/final_candidates.fa
    python get_distinct_isoforms.py $sam_folder/53_isoforms_blat_w_header.sam $sam_folder/$sample\_blat_w_header.sam $output/$sample $sample --fasta $fasta/$sample/final_candidates.fa
done

rm $output/all_samples.tsv
cat $output/*/*.tsv > $output/all_samples.tsv
python plot_hits.py $output/all_samples.tsv $transcript_folder/53_isoforms/final_candidates.fa ~/tmp/heatmap_test_53xxxx.pdf
