#!/bin/bash

### ANNOTATION FOR PDX TUMOUR ONLY SAMPLES
## This script runs annovar on all the vcfs of the PDXs

## INPUT DIRECTORIES:
WD=/scratchb/cclab/masina01/PDX_annotate_Jul2022
VCF_DIR=/scratchb/cclab/masina01/PDX_Jul2022_WES/Germline_calling_merged_Novoalign_bam
SPREADSHEET=$WD/sample_list.txt


## OUTPUT DIRECTORIES:
OD=$WD
VCF_ANNOVAR_DIR=$OD/Vcf_annovar
VCF_FINALCALLS_DIR=$OD/Vcf_finalcalls
TEMP=/tmp/PDX_anno_tmp_$RANDOM
mkdir -p $VCF_ANNOVAR_DIR      #create dir if it does not exist
mkdir -p $VCF_FINALCALLS_DIR   #create dir if it does not exist
mkdir -p $TEMP                 #create dir if it does not exist


## REQUIRED FILES:
PON=/scratchb/cclab/masina01/common_files/PON_BCs_and_agnese_no_replicates.vcf.gz
SUFFIX='germline.vcf.gz'

## SOFTWARE:
ANNOVAR_PATH=/scratchb/cclab/masina01/Software/annovar/
PATH=$ANNOVAR_PATH:$PATH
R_SCRIPTS=$WD/Scripts

FILENAMES=$(cat $SPREADSHEET)
FILENAME=$(echo $FILENAMES | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)
VCF=$VCF_DIR/$FILENAME.$SUFFIX

### ANNOTATE WITH ANNOVAR:

# CHANGE VCF NAME TO REMOVE () - OTHERWISE ANNOVAR RUNS ERROR:
FILENAME_NO_BRACKETS=$(echo $FILENAME | tr '(' '-' | tr ')' '-' )
cp $VCF $TEMP/$FILENAME_NO_BRACKETS.$SUFFIX
cp $VCF.tbi $TEMP/$FILENAME_NO_BRACKETS.$SUFFIX.tbi
VCF_NO_BRACKETS=$TEMP/$FILENAME_NO_BRACKETS.$SUFFIX

# CREATE LOCAL TEMP FOR ANNOVAR TEMPORARY FILES:
LOCAL_TEMP="/tmp/temp_$RANDOM"
mkdir "$LOCAL_TEMP"

# RUN ANNOVAR:
table_annovar.pl $VCF_NO_BRACKETS \
		"$ANNOVAR_PATH"humandb/ \
		-buildver hg19 \
		-out $VCF_ANNOVAR_DIR/$FILENAME_NO_BRACKETS.annovar \
		-remove \
		-protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,snp142,dbnsfp33a,clinvar_20160302,icgc21,cosmic81_coding_reformatted,cosmic81_noncoding_reformatted \
		-operation g,r,r,f,f,f,f,f,f,f,f \
		-nastring . \
		--onetranscript \
		--tempdir $LOCAL_TEMP \
		-vcfinput

VCF_ANNOVAR=$VCF_ANNOVAR_DIR/$FILENAME_NO_BRACKETS.annovar.hg19_multianno.vcf

# INFORMATION ABOUT THE DATABASES:
#dbnsfp33a > non synonymous variant annotation (updated version of ljb26_all)
#refGene > gene
#cytoBand > cytoband
#genomicSuperDups > annotation to genomic variants in segmental duplications
#esp6500siv2_all > exome sequencing project annotation
#snp142 > dbSNP annotation
#Additional disease-specific variants:
# clinvar_20160302 (f) > clinically relevant information about the variant
# icgc21 (f) International Cancer Genome Consortium version 21
# cosmic81_coding_reformatted (f)
# cosmic81_noncoding_reformatted (f)> COSMIC is no longer supported by annovar so I had to add it manually following the instructions of the ANNOVAR website. I have split it into 2 databases, and added the annotation. "COSMIC includes somatic mutations reported in literature in various types of cancers"
