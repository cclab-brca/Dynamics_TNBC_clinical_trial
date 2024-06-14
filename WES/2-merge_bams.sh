#!/bin/bash

# MERGED BAM TO FINAL BAM FOR SAMPLE

#RUN JOB
#SBATCH --array=1-20
#SBATCH --job-name=MergedBam
#SBATCH --mem=24000
#SBATCH --ntasks-per-node=4
#SBATCH --output=/scratchb/cclab/giovan01/log/MergedBam/MergedBam.%A.%a.out.txt
#SBATCH --error=/scratchb/cclab/giovan01/log/MergedBam/MergedBam.%A.%a.err.txt

#################################### INPUT INFO ####################################

# DIRECTORIES
WD=/scratcha/cclab/giovan01/PDX/Exome
FD=/home/giovan01/PDX/Exome/FastqToSample
BAM_NOVO=$WD/BAM_files/Novoalign
TEMP_BAM=/tmp/BAM_aligned_merged
REPORTS_NOVO_ALIGNMENT_MERGED=$WD/REPORTS_NOVOALIGN/Alignment_stats_merged_bam
REPORTS_NOVO_BAM_STATS_MERGED=$WD/REPORTS_NOVOALIGN/BAM_stats_merged_bam
REPORTS_NOVO_INSERT_SIZE=$WD/REPORTS_NOVOALIGN/Insert_size_merged_bam
REPORTS_NOVO_READCOUNT=$WD/REPORTS_READCOUNT_MERGED_FILES/NOVOALIGN
REPORTS_MIRNA_NOVO_STATS=$WD/REPORTS_MIRNA_NOVO_STATS
GERM_DIR=$WD/Germline_calling_merged_Novoalign_bam
EXOME_DIRS=$BAM_NOVO:$REPORTS_NOVO_ALIGNMENT_MERGED:$REPORTS_NOVO_BAM_STATS_MERGED:$REPORTS_NOVO_INSERT_SIZE:$REPORTS_NOVO_READCOUNT:$REPORTS_MIRNA_NOVO_STATS:$GERM_DIR

# Ensure all paths exist. Create directory structure if not present
ARR=$(echo $EXOME_DIRS | tr ":" "\n")
for x in $ARR
do
        if [ ! -d "$x" ]; then
        	mkdir -p $x
        fi

done

# REQUIRED FILES
TARGET_BED=/scratchb/cclab/giovan01/common_files/Nextera_v1.2_target.interval_list    				#Pay attention to the extension required by picard. It should be target.interval_list. File created through BedToIntervals (provided by Picard)
TARGET_MIRNA_BED=/scratchb/cclab/giovan01/common_files/nextera_v1.2_plus_miRBase_V20_target.interval_list	#Pay attention to the extension required by picard.	It should be target.interval_list. File created through BedToIntervals (provided by Picard)
MIRNA_BED=/scratchb/cclab/giovan01/common_files/miRBase_v20_coords_target.interval_list				#Pay attention to the extension required by GATK. It should be target.interval_list. File created through BedToIntervals (provided by Picard)
GENOME_COMPLETE=/scratchb/cclab/giovan01/Genomes/human_g1k_v37_decoy_UCSC.mmu10.fa
GENOME=/scratchb/cclab/giovan01/Genomes/human_g1k_v37_decoy.fasta
SPREADSHEET=$FD/FastqToSample_july_2021_no_header.txt
KNOWN_SITES_GOLD_STANDARD=/scratchb/cclab/giovan01/Software/gatk/File_gatk/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
KNOWN_SITES_DBSNP=/scratchb/cclab/giovan01/Software/gatk/File_gatk/dbsnp_138.b37.vcf.gz
KNOWN_SITES_1K_GENOMES=/scratchb/cclab/giovan01/Software/gatk/File_gatk/1000G_phase1.indels.b37.sites.vcf.gz

# SOFTWARE PATH: check that they are in $PATH
JAVA_PATH=/scratchb/cclab/giovan01/Software/java/jre1.8.0_171/bin/
SAMTOOLS_PATH=/scratchb/cclab/giovan01/Software/samtools/samtools-1.8/
NOVOALIGN_PATH=/scratchb/cclab/giovan01/Software/novocraft/novocraft/
PICARD_PATH=/scratchb/cclab/giovan01/Software/picard/
GATK_PATH=/scratchb/cclab/giovan01/Software/gatk/gatk-4.0.4.0/

# OTHER INFO
#PLATFORM="Illumina"
MEMORY="16g"

####################################################################################

#Select the array (file without header)
SAMPLES_ID=$(cut -f6 $SPREADSHEET | sort | uniq)				
SAMPLE_ID=$(echo $SAMPLES_ID | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)
echo $SAMPLE_ID
FILENAMES=$(grep "$SAMPLE_ID" $SPREADSHEET | cut -f1)                          #Defines filename list from first column of spreadsheet
FILENAME=$(echo $FILENAMES | tr " " "\n")
FILES_NUMBER=$(echo $FILENAMES | tr " " "\n" | wc -l)				#Counts the number of filenames for each sample

#Generate folder on /tmp/                                                       #Creates a new folder in /tmp/ where data will be stored
TEMP_FULL_NAME_BAM="${TEMP_BAM}_${SAMPLE_ID}"
mkdir $TEMP_FULL_NAME_BAM

########### Merging Novoalign files and Haplotype Caller ###########

echo "[exome:single_bam_to_merged_bam $(date +"%Y-%m-%d %T")] Starting Novoalign merging pipeline"

#Copy BAM - Novoalign files
FILENAMES_ARRAY=($FILENAMES)
BAM_FILES=${FILENAMES_ARRAY[@]/%/".sorted.ba*"}
BAM_FILES_ARRAY=($BAM_FILES)
BAM_FILES_ALL=${BAM_FILES_ARRAY[@]/#/$BAM_NOVO/}
cp $BAM_FILES_ALL $TEMP_FULL_NAME_BAM

#Merge BAM files Novoalign
echo "[exome:single_bam_to_merged_bam $(date +"%Y-%m-%d %T")] Merging sample $SAMPLE_ID (Novoalign)"
BAM_FILENAMES=${FILENAMES_ARRAY[@]/%/".sorted.bam"}
BAM_FILENAMES_ARRAY=($BAM_FILENAMES)
BAMLIST=${BAM_FILENAMES_ARRAY[@]/#/$TEMP_FULL_NAME_BAM/}

"$NOVOALIGN_PATH"novosort -m $MEMORY \
                          -t $TEMP_FULL_NAME_BAM \
                          -c $SLURM_TASKS_PER_NODE \
                          --keeptags \
                          $BAMLIST \
                          -i \
                          -o $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam

FILENAMES_TO_BE_REMOVED=${FILENAMES_ARRAY[@]/%/".sorted.ba*"}
FILENAMES_TO_BE_REMOVED_ARRAY=($FILENAMES_TO_BE_REMOVED)
rmBAMLIST=${FILENAMES_TO_BE_REMOVED_ARRAY[@]/#/$TEMP_FULL_NAME_BAM/}

FILES_MERGED=$(samtools view -H $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam | grep "^@RG" | wc -l)

rm $rmBAMLIST

#Check for number files merged

[[ "$FILES_NUMBER" == "$FILES_MERGED" ]] && echo "[exome: merged all Novoalign bam files selected by quality]" || echo "[Error: not all Novoalign files selected have been merged to generate the final bam]"

#Read count
echo "[exome:read count on Novoalign merged files $(date +"%Y-%m-%d %T")] Compute Read Count"
READCOUNT=$(echo $("$SAMTOOLS_PATH"samtools view -F 4 -q 1 $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam | grep -v $'\t''m.chr' | wc -l) \
$("$SAMTOOLS_PATH"samtools view -F 4 -q 1 $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam | grep $'\t''m.chr' | wc -l))
echo "$READCOUNT" > $REPORTS_NOVO_READCOUNT/$SAMPLE_ID.readcount.txt

#Selecting subset human reads
echo "[exome:subset of human reads on merged files $(date +"%Y-%m-%d %T")] Starting subsetting of human reads"
"$SAMTOOLS_PATH"samtools view -h $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam \
                              | grep -v $'\t''m.chr' \
                              | grep -v 'SN:m.chr' \
                              | "$SAMTOOLS_PATH"samtools view -b - > $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam
                              "$SAMTOOLS_PATH"samtools index $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam

rm $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.bam

#Alignment stats	                                                                   #Calculates alignment stats on reference genome
echo "[exome:bamstats $(date +"%Y-%m-%d %T")] Starting CollectAlignmentMetrics on Novoalign files"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                       CollectAlignmentSummaryMetrics \
                       R=$GENOME \
                       I=$TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam \
                       O=$REPORTS_NOVO_ALIGNMENT_MERGED/$SAMPLE_ID.merged.hg.alignMetrics.txt

#Bam stats                                                                                 #Calculates BAM stats on reference genome and target regions
echo "[exome:bamstats on merged bam $(date +"%Y-%m-%d %T")] Starting CollectHsMetrics on Novoalign merged files - target regions"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                   CollectHsMetrics \
                   I=$TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam \
                   O=$REPORTS_NOVO_BAM_STATS_MERGED/$SAMPLE_ID.merged.hg.bamMetrics.txt \
                   REFERENCE_SEQUENCE=$GENOME \
                   BAIT_INTERVALS=$TARGET_BED \
                   TARGET_INTERVALS=$TARGET_BED

#Bam stats on miRNA sequences                                                              #Calculates BAM stats on reference genome and target regions
echo "[exome:bamstats on merged bam $(date +"%Y-%m-%d %T")] Starting CollectHsMetrics on Novoalign merged files - miRNA regions"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                   CollectHsMetrics \
                   I=$TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam \
                   O=$REPORTS_MIRNA_NOVO_STATS/$SAMPLE_ID.merged.hg.bamMetrics.txt \
                   REFERENCE_SEQUENCE=$GENOME \
                   BAIT_INTERVALS=$MIRNA_BED \
                   TARGET_INTERVALS=$MIRNA_BED

#Bam stats Insert Size
echo "[exome:insertSize on merged bam $(date +"%Y-%m-%d %T")] Starting CollectInsertSizeMetrics on Novoalign merged files"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                   CollectInsertSizeMetrics \
                   I=$TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam \
                   H=$REPORTS_NOVO_INSERT_SIZE/$SAMPLE_ID.merged.hg.InsertSizeMetrics_histogram.pdf \
                   O=$REPORTS_NOVO_INSERT_SIZE/$SAMPLE_ID.merged.hg.InsertSizeMetrics.txt

# Checking reads on chromosome X 
numberInX=$("$SAMTOOLS_PATH"samtools idxstats $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam | grep "^X" | cut -f3)
threshold=10000
if [ "$numberInX" -lt "$threshold" ] 
then 
    echo "[chrX has less reads than expected in Novoalign file]"
fi

# Samtool reheader of Novoalign files

"$SAMTOOLS_PATH"samtools view -H $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam > $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg_header.sam  			# extract header only
sed "s/	SM:.*/	SM:$SAMPLE_ID/" $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg_header.sam > $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg_header_corrected.sam
"$SAMTOOLS_PATH"samtools reheader $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg_header_corrected.sam $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam > $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.corrected.hg.bam
samtools index $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.corrected.hg.bam

mv $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.corrected.hg.bam $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam
mv $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.corrected.hg.bam.bai $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam.bai

# HaplotypeCaller
echo "[exome:germline variant calling on merged bam $(date +"%Y-%m-%d %T")] Starting HaplotypeCaller on Novoalign merged files"
"$GATK_PATH"gatk --java-options "-Xmx$MEMORY" HaplotypeCaller \
                 -R $GENOME \
                 -I $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.bam \
                 -L $TARGET_MIRNA_BED \
                 --min-base-quality-score 20 \
                 -stand-call-conf 30.0 \
                 --annotate-with-num-discovered-alleles  \
                 -O $GERM_DIR/$SAMPLE_ID.germline.vcf.gz

mv $TEMP_FULL_NAME_BAM/$SAMPLE_ID.dedup.merged.hg.ba* $MERGED_BAM_NOVO

rm -r "$TEMP_FULL_NAME_BAM"

echo "[exome:Merged BAM $(date +"%Y-%m-%d %T")] Merge BAM Novoalign complete"

