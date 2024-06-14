# WES Alignment and mutation calling:

This folder contains the code and processed data required to reproduce the WES analysis done in this study. 

The code in **1-fq_to_bam_alignment.sh** begins by downloading fastq files from Clarity (the laboratory management system used by the Genomics Core facility in the CRUK Cambridge Institute). It then runs the alignment on the Combined human-mouse Reference Genome using novocraft. Next, in **2-merge_bams.sh** the BAM files belonging to the same sample are merged, mouse contamination is removed, bamstats are generated. Lastly, HaplotypeCaller is run on all bams. 

## Requirements
The following software is required to reproduce the analysis:

* cutadapt
* TrimGalore-0.4.5
* samtools-1.8
* novocraft
* gatk-4.0.4.0

We use **Picard** to generate and collect alignment stats. If you wish to reproduce that step we have provided the required interval_list BED files in data: 

* Nextera37_hg19_mmu10.target.intervals - contains target regions for the combined human and mouse genome. 
* Nextera_v1.2_target.interval_list - contains target regions for the human genome only

## Supplied data

The *data* folder contains the following:

* FastqToBam.txt - table containing information regarding flowcell, pool, barcode. This is used as input for alignment and variant calling. 
* bamstats.txt - table containing alignment stats of the merged bams.
* Target intervals files regarding the target regions for teh combined human and mouse genome, and human genome only, respectively, as described above. 

The *intermediate_files* folder contains the following:

* Intermediate bam stats generated from single bam files aligned with Novoalign prior to merging. 
* /Alignment_stats contains stats generated using Picard's CollectAlignmentSummaryMetrics (alignMetrics.txt) and Novocraft's novoalign (novostats.txt). 
* /BAM_stats contains stats generated using Picard's CollectHsMetrics.

