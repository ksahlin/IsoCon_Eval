#!/usr/bin/env bash
python --version
root="/Users/kxs624/Dropbox/IsoCon/REVISION_DATA/ampliconic_analysis/PIPELINE_INTERMEDIATE_DATA"
reffolder="/Users/kxs624/Dropbox/IsoCon/REVISION_DATA/ampliconic_analysis"
# transcript_folder="/Users/kxs624/Dropbox/IsoCon/REVISION_DATA/ampliconic_analysis/PIPELINE_INTERMEDIATE_DATA"
psl_folder=$root/"blat_psl_files"
sam_folder=$root/"sam_files"
fasta_batches=$root/"fasta_batches"
output="/Users/kxs624/tmp/ampliconic_analysis/analysis_output_test"
# fasta="/Users/kxs624/tmp/FMR1/"
        
samples="shared" #"sample1_specific sample2_specific" #   "both_DB shared" #"sample1_specific sample2_specific shared"

for sample in $samples
do
    echo $sample
    # python split_BLAT_batches.py $root/$sample.fa $fasta_batches/$sample
    
    # rm $psl_folder/$sample/all.psl
    # awk 'FNR==1{print ""}1' $psl_folder/$sample/*.psl > $psl_folder/$sample/all.psl
    # echo uncle_psl.py -N 10 -f $root/$sample.fa $psl_folder/$sample/all.psl $sam_folder/$sample\_blat.sam
    # uncle_psl.py -N 10 -f $root/$sample.fa $psl_folder/$sample/all.psl $sam_folder/$sample\_blat.sam
    # cat $root/sam_header_tabs.sam $sam_folder/$sample\_blat.sam > $sam_folder/$sample\_blat_w_header.sam

    echo python get_distinct_isoforms_ampl.py $sam_folder/both_DB_blat_w_header.sam $sam_folder/$sample\_blat_w_header.sam $output $sample --query_fasta $root/$sample.fa --ref_fasta $reffolder/both_DB.fa
    python get_distinct_isoforms_ampl.py $sam_folder/both_DB_blat_w_header.sam $sam_folder/$sample\_blat_w_header.sam $output/$sample $sample --query_fasta $root/$sample.fa --ref_fasta $reffolder/both_DB.fa
done

rm $output/all_samples_best_hits.tsv
echo -e "q_acc\tread_support\talignment_id\talignment_coverage\tsample" > $output/all_samples_best_hits.tsv
cat $output/*/*_best_hits.tsv >> $output/all_samples_best_hits.tsv
# python ../ampliconic/plot_quality.py $output/all_samples_best_hits.tsv $output/dotplot_test.pdf
