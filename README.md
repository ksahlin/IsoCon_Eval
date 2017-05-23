# IsoCon_Eval
Evaluation scripts for IsoCon + targeted sequencing approach.

Folder structure should mimic the format in [evaluation google doc](https://docs.google.com/document/d/1Qx_jhRixbQi0Y1_QCYXpn9T5Z8NXLsAPb6-_bd7LDdw/edit) under "Proposed Analyses
"

# "Nucleotide-level accuracy":

Below is short documentation of the scripts for this analysis.

## Running analysis "captured varinats" (see gdoc)

The script does the following

```
1. Get illumina-supported varinats from alignments of illumina reads to consensus.
2. For every variant of exp1:
    Identify illumina reads supporting the variant.
    Check their alignments in exp 2 (google doc).
3. Summarise output:
    A variant is captured if at least one illumina read aligns perfect (match) the position that signaled a varinat in the consensus, otherwise it's non-captured.
    Non-captured variants are FN  (caveat: FN that are not present in the reference database are not counted.) 
    We store the coverage for captured/non-captured variants and plot it as well as output information to stdout.
```


Example run:
```
python get_illumina_supported_variants.py \
-illumina_to_ref /PATH/TO/<GENE_NAME>_sample1_aln_after_rmdup.bam \
-illumina_to_pred /PATH/TO/<METHOD>/<GENE_NAME>_sample1_aln_after_rmdup.bam \
-reference /PATH/TO/<GENE_NAME>_consensus.fa \
-predicted /PATH/TO/<METHOD>/<GENE_NAME>_sample1.fa \
-outfolder  /PATH/TO/OUTPUT/<METHOD>_<GENE_NAME>_sample1 \
--variant_cutoff 2
```
    
Here, it is expected that the folder sturcture will have appropriate names. For example, `<METHOD>` should be something like, Tofu, IsoCon-v0.1.0. `<GENE_NAME>` is the gene consensus name.  The output will be in the folder `/ROOT/PATH/TO/OUTPUT/TOFU_TSPY_sample1/`. In addition, important summary stats are written to output. An example output is shown below.


```
FINAL STATS FROM RUN


Coverage over positions: mean 3877.15021008, stddev: 2939.5857172, median: 3636.0, min: 2, max:7999
Based on coverage distribution, suggested variant support cutoff is: 16.0 


Total Illumina reads that supported a variant on consensus: 8742
Total variant positions on consensus with coverage of at least 2 illumina reads: 1478
Number of Illumina reads that upported a varinat on consensus but were not mapped to any predicted transcript: 5


Total deletion sites on consensus with coverage of at least 2 illumina reads: 27:
Total insertion sites on consensus with coverage of at least 2 illumina reads: 9:
Total substitution sites on consensus with coverage of at least 2 illumina reads: 1442:


Deletion sites captured in predicted with coverage of at least 2 illumina reads: 2:
Insertion sites captured in predicted with coverage of at least 2 illumina reads: 2:
Substitution sites captured in predicted with coverage of at least 2 illumina reads: 10:


Deletion sites not captured in predicted with coverage of at least 2 illumina reads: 25:
Insertion sites not captured in predicted with coverage of at least 2 illumina reads: 7:
Substitution sites not captured in predicted with coverage of at least 2 illumina reads: 1432:


deletions captured depths on reference:
[8, 3]
deletions captured depths on predicted:
[6, 1]
deletions not captured depths:
[4, 5, 2, 2, 13, 4, 6, 2, 3, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2, 2, 7, 7, 4]


insertions captured depths on reference:
[94, 5]
insertions captured depths on predicted:
[90, 5]
insertions not captured depths:
[2, 2, 2, 5, 2, 2, 2]


substitutions captured depths on reference:
[58, 233, 73, 4, 1498, 11, 1519, 6, 32, 15]
substitutions captured depths on predicted:
[34, 178, 22, 2, 1285, 1, 1311, 1, 17, 11]
substitutions not captured depths:
[2, 20, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 4, 3, 2, 2, 8, 8, 2, 3, 2, 3, 2, 2, 2, 2, 2, 3, 3, 2, 4, 2, 2, 4, 2, 2, 2, 3, 2, 4, 2, 2, 3, 3, 2, 4, 2, 4, 10, 2, 2, 6, 2, 3, 3, 2, 3, 5, 2, 4, 3, 6, 3, 7, 2, 2, 5, 2, 4, 2, 3, 2, 3, 4, 35, 4, 28, 2, 17, 6, 2, 3, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 4, 2, 2, 4, 3, 2, 2, 2, 3, 2, 3, 6, 2, 72, 46, 94, 2, 136, 5, 122, 5, 8, 9, 3, 13, 2, 2, 108, 22, 9, 4, 2, 9, 28, 7, 3, 12, 5, 5, 5, 2, 34, 2, 53, 20, 2, 2, 6, 29, 5, 3, 28, 63, 6, 2, 4, 2, 2, 2, 23, 3, 2, 4, 155, 7, 4, 2, 2, 2, 52, 2, 4, 3, 2, 3, 2, 86, 6, 19, 2, 2, 4, 9, 4, 8, 8, 4, 2, 2, 4, 4, 106, 4, 6, 10, 4, 2, 10, 2, 2, 6, 4, 93, 2, 5, 14, 6, 7, 2, 5, 4, 3, 27, 2, 6, 2, 3, 498, 3, 6, 2, 2, 2, 72, 2, 5, 3, 6, 5, 3, 83, 4, 3, 5, 24, 3, 10, 17, 6, 5, 3, 2, 16, 5, 6, 2, 2, 5, 10, 2, 4, 2, 3, 3, 2, 6, 75, 10, 4, 13, 4, 28, 4, 12, 5, 2, 2, 3, 12, 98, 3, 2, 3, 6, 93, 3, 3, 2, 18, 26, 22, 2, 6, 2, 81, 8, 3, 3, 4, 16, 11, 2, 2, 7, 7, 2, 10, 87, 7, 10, 5, 17, 7, 13, 4, 11, 3, 2, 7, 3, 7, 5, 3, 5, 11, 40, 7, 3, 2, 2, 32, 6, 3, 3, 2, 3, 4, 34, 3, 9, 2, 4, 3, 12, 49, 7, 2, 19, 4, 2, 2, 6, 4, 4, 43, 3, 3, 49, 5, 16, 7, 6, 7, 5, 6, 2, 42, 2, 4, 17, 7, 4, 10, 2, 3, 14, 3, 6, 2, 4, 4, 28, 3, 9, 4, 2, 12, 11, 2, 2, 4, 2, 4, 3, 16, 83, 2, 3, 12, 4, 2, 13, 3, 2, 16, 9, 3, 4, 2, 28, 20, 10, 9, 7, 13, 7, 50, 2, 2, 3, 21, 7, 4, 2, 2, 6, 2, 53, 2, 4, 9, 3, 2, 7, 2, 5, 35, 8, 3, 2, 2, 4, 14, 25, 2, 12, 7, 7, 13, 2, 3, 2, 5, 3, 7, 6, 2, 4, 2, 2, 47, 4, 13, 6, 4, 4, 54, 2, 14, 2, 8, 3, 3, 5, 13, 13, 5, 3, 6, 5, 3, 15, 6, 2, 2, 7, 26, 3, 7, 2, 2, 3, 11, 2, 23, 7, 4, 7, 4, 12, 10, 17, 8, 5, 5, 7, 6, 5, 4, 9, 8, 2, 2, 20, 44, 20, 4, 21, 7, 6, 13, 6, 7, 6, 6, 4, 16, 13, 9, 35, 8, 2, 19, 3, 2, 3, 4, 6, 10, 4, 2, 4, 16, 33, 6, 10, 3, 12, 7, 7, 3, 46, 6, 3, 4, 7, 16, 4, 3, 2, 2, 4, 3, 56, 10, 10, 5, 7, 4, 25, 3, 2, 4, 10, 4, 5, 6, 7, 2, 3, 5, 6, 9, 12, 2, 2, 50, 4, 12, 4, 2, 3, 5, 4, 3, 4, 11, 18, 7, 3, 10, 5, 9, 8, 16, 8, 6, 10, 10, 4, 7, 2, 24, 10, 5, 2, 13, 10, 3, 16, 11, 8, 10, 8, 17, 4, 22, 3, 7, 9, 10, 5, 4, 6, 7, 3, 19, 3, 19, 4, 10, 15, 17, 3, 14, 4, 3, 4, 6, 34, 22, 6, 9, 24, 4, 19, 9, 5, 5, 11, 3, 6, 2, 17, 4, 29, 6, 2, 23, 11, 5, 20, 9, 32, 4, 17, 9, 14, 2, 14, 2, 3, 13, 4, 12, 11, 3, 21, 8, 14, 3, 3, 9, 13, 3, 8, 4, 3, 14, 2, 9, 6, 8, 18, 5, 7, 3, 3, 2, 13, 15, 8, 18, 4, 33, 10, 10, 5, 2, 29, 5, 3, 4, 2, 16, 10, 8, 2, 8, 2, 13, 4, 4, 18, 3, 5, 8, 18, 7, 22, 5, 3, 3, 5, 5, 25, 2, 2, 8, 6, 8, 4, 16, 3, 20, 3, 7, 14, 4, 4, 17, 7, 13, 2, 14, 2, 6, 4, 3, 5, 12, 3, 7, 6, 17, 2, 5, 8, 4, 3, 9, 8, 27, 3, 6, 7, 2, 9, 21, 45, 7, 6, 2, 2, 5, 5, 437, 4, 5, 18, 6, 2, 6, 14, 3, 3, 12, 2, 4, 6, 2, 21, 5, 4, 3, 18, 2, 22, 4, 2, 2, 4, 11, 5, 2, 4, 17, 6, 5, 3, 4, 9, 10, 3, 5, 5, 3, 26, 6, 3, 17, 26, 12, 4, 2, 10, 5, 15, 8, 5, 6, 14, 13, 5, 12, 2, 4, 9, 3, 9, 8, 2, 8, 8, 3, 19, 6, 2, 13, 6, 6, 15, 4, 21, 3, 14, 5, 5, 20, 3, 12, 3, 3, 3, 3, 9, 4, 4, 4, 2, 5, 5, 3, 20, 3, 12, 5, 2, 20, 9, 10, 3, 9, 5, 6, 4, 9, 3, 7, 5, 9, 11, 2, 2, 22, 7, 3, 2, 10, 4, 2, 9, 9, 15, 3, 6, 11, 6, 7, 2, 4, 3, 3, 8, 2, 6, 5, 10, 8, 15, 7, 9, 5, 13, 9, 8, 5, 2, 11, 7, 18, 10, 3, 3, 3, 3, 3, 8, 2, 5, 4, 12, 4, 8, 8, 6, 2, 34, 14, 6, 3, 3, 7, 2, 3, 13, 4, 7, 2, 8, 9, 2, 5, 12, 2, 10, 11, 7, 5, 2, 18, 4, 3, 11, 6, 10, 4, 6, 2, 4, 4, 15, 5, 12, 7, 5, 21, 3, 7, 2, 11, 8, 4, 18, 4, 8, 2, 30, 2, 6, 3, 23, 2, 5, 3, 10, 4, 21, 4, 7, 7, 5, 3, 19, 3, 2, 3, 10, 4, 8, 2, 3, 26, 4, 4, 2, 9, 8, 2, 9, 2, 5, 15, 4, 9, 3, 4, 2, 4, 36, 6, 7, 22, 6, 6, 3, 2, 13, 9, 10, 3, 5, 3, 9, 11, 10, 4, 3, 15, 6, 10, 4, 7, 6, 8, 2, 14, 3, 3, 9, 14, 4, 9, 4, 4, 5, 7, 4, 9, 10, 3, 3, 3, 2, 15, 2, 10, 6, 21, 4, 6, 12, 5, 2, 5, 6, 3, 2, 2, 10, 4, 13, 5, 12, 2, 14, 14, 7, 4, 4, 14, 2, 2, 2, 4, 2, 4, 2, 9, 9, 3, 2, 5, 2, 3, 11, 2, 27, 3, 7, 3, 5, 3, 5, 3, 3, 2, 7, 8, 4, 28, 10, 5, 6, 3, 5, 10, 13, 5, 135, 8, 8, 3, 5, 10, 6, 15, 6, 11, 10, 4, 6, 4, 6, 5, 11, 2, 2, 2, 2, 20, 5, 17, 9, 2, 2, 2, 4, 2, 4, 5, 5, 3, 4, 8, 4, 8, 2, 15, 3, 2, 8, 6, 2, 2, 2, 9, 4, 18, 8, 3, 3, 7, 2, 3, 8, 6, 2, 9, 8, 2, 9, 3, 9, 14, 5, 5, 6, 5, 3, 2, 4, 3, 9, 9, 3, 2, 8, 6, 7, 2, 8, 8, 6, 4, 5, 2, 10, 2, 4, 10, 14, 4, 3, 6, 3, 3, 5, 2, 6, 7, 18, 2, 17, 2, 7, 15, 7, 8, 2, 14, 2, 8, 2, 12, 3, 3, 5, 7, 8, 2, 5, 3, 12, 7, 2, 3, 3, 17, 2, 2, 6, 5, 3, 2, 4, 7, 4, 2, 2, 5, 3, 4, 2, 10, 2, 10, 6, 3, 6, 2, 7, 4, 4, 6, 6, 2, 3, 4, 2, 5, 2, 2, 7, 5, 5, 4, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 6, 7, 2, 5, 7, 3, 2, 2, 2, 4, 2, 2, 3, 5, 7, 3, 3, 4, 4, 6, 12, 3, 3, 3, 3, 6, 2, 2, 4, 10, 4, 11, 5, 4, 7, 2, 5, 4, 5, 2, 5, 3, 2, 2, 10, 3, 4, 11, 2, 13, 11, 2, 5, 3, 3, 4, 5, 2, 3, 2, 7, 4, 5, 3, 5, 2, 3, 8, 5, 4, 3, 3, 4, 2, 8, 2, 8, 5, 2, 4, 2, 3, 3, 2, 3, 4, 2, 2]
```

## Running analysis "Unsupported positions" (see gdoc)

The script does the following

```
1. For each position in the predicted transcripts, check illumina support.
2. If there are less that 2 reads that matches the predicted transcript at a given position: clasify the position as unsupported.
```

Example run:
```
python get_illumina_supported_variants.py \
-illumina_to_pred /PATH/TO/<METHOD>/<GENE_NAME>_sample1_aln_after_rmdup.bam \
-predicted /PATH/TO/<METHOD>/<GENE_NAME>_sample1.fa \
-outfolder  /PATH/TO/OUTPUT/<METHOD>_<GENE_NAME>_sample1 \
--unsupported_cutoff 2
```

Example uptput

```
('number of references:', 22)
('References not seen in pileup:', 0)
Total number of positions in predicted transcripts:21245
('TOTAL UNSUPPORTED POSITIONS (Illumina support < 2, masking first and last 21 positions in predicted transcripts due to barcodes):', 660)
Percentage of supported bases: 97.0

```

# Exon level accuracy

# Isoform level accuracy
