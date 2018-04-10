#!/usr/bin/env bash

root="/Users/kxs624/tmp/FMR1/analysis/Pauls"
transcript_folder="/Users/kxs624/tmp/FMR1"
psl_folder="/Users/kxs624/tmp/FMR1/analysis/Pauls/blat_psl_files"
sam_folder="/Users/kxs624/tmp/FMR1/analysis/Pauls/sam_files"
output="/Users/kxs624/tmp/FMR1/analysis/Pauls/analysis_output"

# samples= "53_isoforms 1009 1015 1018 4549 5123 5248"
samples="53_isoforms 1009"

for sample in $samples
do
    rm $psl_folder/$sample\_all.psl

    awk 'FNR==1{print ""}1' $psl_folder/$sample\_*.psl > $psl_folder/$sample\_all.psl

    uncle_psl.py -N 50 -f $transcript_folder/$sample/final_candidates.fa $psl_folder/$sample\_all.psl $sam_folder/$sample\_blat.sam

    cat $root/sam_header_tabs.sam $sam_folder/$sample\_blat.sam > $sam_folder/$sample\_blat_w_header.sam

    python get_distinct_isoforms.py $sam_folder/53_isoforms_blat_w_header.sam $sam_folder/$sample\_blat_w_header.sam $output $sample
done

#cat $output/*.tsv > $output/all_samples.tsv
#python plot_hits.py
