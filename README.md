# IsoCon_Eval
Repository for running experiments with snakemake pipeline. Evaluation scripts for IsoCon.


### Installation
Install snakemake (python3 only) [link](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html?highlight=installation). 


### Simulated data

This experiment pipeline assumes you have IsoCon and ICE installed. You can however remove one of the tools as described in next section

#### Change experiment settings

You can remove tools or experiments by changing the fields in the config file `simulated_experiments.json`. To remove an experiment, change `"EXPERIMENTS" : ["MEMBER_EXPERIMENT", "ISOFORM_EXPERIMENT", "COMBINED_EXPERIEMENT" ],`, similarly, you can alter these fields to reduce the number or size of the experiments:

```
    "GENE_MEMBERS" : ["HSFY2", "DAZ2", "TSPY13P"],
    "READ_COUNTS" : [ "20", "100", "500", "2500", "12500"],
    "RUNS" : ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
    "ABUNDANCE" : ["constant", "exponential"],
    "TOOLS" : ["ISOCON", "ICE"],
``` 

Modify all paths in `simulated_experiments.json`, then try a dry run as

```
snakemake  --configfile simulated_experiments.json --latency-wait 100 --verbose  -n
```

To start the actual pipeline, run the above command without parameter `-n`.  This pipeline should be run on a cluster, below is a command to distribute jobs. 

```
snakemake  --keep-going -j 999999 --cluster "sbatch --exclude={cluster.exclude} -c {cluster.ntasks} -N {cluster.Nodes}  -t {cluster.runtime} -J {cluster.jobname} --mail-type={cluster.mail_type} --mail-user={cluster.mail}" --cluster-config cluster.json --configfile simulated_experiments.json --latency-wait 100 --verbose  -n
```

### Experimental data

The experiment pipeline assumes you have the software ToFU, ccs, and bax2bam installed to preprocess the Iso-Seq data (these three are all included in pacbio SMRTLink suite). IsoCon, ICE, and proovread will be needed to run downstream processing of CCS reads but any of these softwar can be removed from the pipeline (therefore does not have to be installed) by mdifying the `TOOLS` field in the config file `experimental_experiments.json` . 

```
    "TOOLS" :["ISOCON", "ICE_QUAL", "PROOVREAD"]
```