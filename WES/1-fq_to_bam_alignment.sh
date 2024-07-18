#!/bin/bash 

# FASTQ TO FINAL BAM FOR SINGLE LANES

#RUN JOB

#SBATCH --array=1-60
#SBATCH --job-name=FqtoBam
#SBATCH --mem=24000
#SBATCH --ntasks-per-node=8
#SBATCH --output=/scratchb/cclab/giovan01/log/FastqtoBam/FastqtoBam.%A.%a.out.txt
#SBATCH --error=/scratchb/cclab/giovan01/log/FastqtoBam/FastqtoBam.%A.%a.err.txt

#################################### INPUT INFO ####################################

# DIRECTORIES
WD=/scratcha/cclab/giovan01/PDX/Exome
FD=/home/giovan01/PDX/Exome/FastqToSample
BAM_NOVO=$WD/BAM_files/Novoalign
TEMP_BAM=/tmp/BAM_aligned
REPORTS_NOVOALIGN_ALIGNMENT=$WD/REPORTS_NOVOALIGN/Alignment_stats
REPORTS_NOVOALIGN_BAM_STATS=$WD/REPORTS_NOVOALIGN/BAM_stats
EXOME_DIRS=$BAM_NOVO:$REPORTS_NOVOALIGN_ALIGNMENT:$REPORTS_NOVOALIGN_BAM_STATS

# Ensure all paths exist. Create directory structure if not present
ARR=$(echo $EXOME_DIRS | tr ":" "\n")
for x in $ARR
do
        if [ ! -d "$x" ]; then
        	mkdir -p $x
        fi

done


# REQUIRED FILES	
TARGET=/scratchb/cclab/giovan01/common_files/Nextera37_b37_mmu10.target.intervals
GENOME=/scratchb/cclab/giovan01/Genomes/human_g1k_v37_decoy_UCSC.mmu10.fa
INDEXED=indexed.*.human_g1k_v37_decoy_UCSC.mmu10
GENOME_INDEXED=/scratchb/cclab/giovan01/Genomes/Indexed/$INDEXED.*
SPREADSHEET=$FD/FastqToSample.txt

# SOFTWARE PATH: check that they are in $PATH
JAVA_PATH=/scratchb/cclab/giovan01/Software/java/jre1.8.0_171/bin/
CLARITY_PATH=/scratchb/cclab/giovan01/Software/
CUTADAPT_PATH=/scratchb/cclab/giovan01/Software/cutadapt/
TRIM_PATH=/scratchb/cclab/giovan01/Software/trim_galore/TrimGalore-0.4.5/
SAMTOOLS_PATH=/scratchb/cclab/giovan01/Software/samtools/samtools-1.8/
NOVOALIGN_PATH=/scratchb/cclab/giovan01/Software/novocraft/novocraft/
PICARD_PATH=/scratchb/cclab/giovan01/Software/picard/

# OTHER INFO
FASTQ1_SUFFIX=".r_1.fq.gz"
FASTQ2_SUFFIX=".r_2.fq.gz" 
PLATFORM="Illumina"
MEMORY="24g"

####################################################################################	

#Select the array (file without header)
FILENAMES=$(cut -f1 $SPREADSHEET)                                                   #Defines filename list from first column of spreadsheet
FILENAME=$(echo $FILENAMES | cut -d ' ' -f $SLURM_ARRAY_TASK_ID)                    #Defines filename that will be analysed considering the array task id (assigned to the job)
SLX_id="$(cut -d '.' -f1 <<< "$FILENAME")"                                          #Parses filename and defines SLX-id
FLOWCELL="$(cut -d '.' -f3 <<< "$FILENAME")"                                        #Parses filename and defines FLOWCELL
BARCODE="$(cut -d '.' -f2 <<< "$FILENAME")"                                         #Parses filename and defines BARCODE
LANE="$(cut -d '.' -f4 <<< "$FILENAME")"                                            #Parses filename and defines LANE
FL_BAR_LN="${FLOWCELL}_${BARCODE}.${LANE}"                                          #Creates new id as defined in SAMv1.pdf (it will be used in BAM header)

#Generate folder on /tmp/                                                          #Creates a new folder in /tmp/ where data will be stored
TEMP_FULL_NAME_BAM="${TEMP_BAM}_${FILENAME}"                                              
mkdir $TEMP_FULL_NAME_BAM 

#Download fastq file from spreadsheet                                               #Downloads fastq files for $FILENAME (the filename corresponding to array task id) in $TEMP_FULL_NAME_BAM
echo "[Downoloading "$FILENAME" files for $SLX_id]"
"$JAVA_PATH"java -jar "$CLARITY_PATH"clarity-tools.jar -l $SLX_id -f "$FILENAME""$FASTQ1_SUFFIX" -f "$FILENAME""$FASTQ2_SUFFIX" -d "$TEMP_FULL_NAME_BAM"

##STARTING ALIGNMENT PIPELINE

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Starting alignment pipeline"     #Written in output (defined from sbatch command)

#Copy genome                                                                        #Copies genome files to $TEMP_FULL_NAME_BAM
cp $GENOME_INDEXED $TEMP_FULL_NAME_BAM              

#Trimming (if you want to generate a fastqc report insert --fastqc before --paired) #Performs trimming on paired end data
"$TRIM_PATH"trim_galore --output_dir $TEMP_FULL_NAME_BAM --nextera --paired "$TEMP_FULL_NAME_BAM"/"$SLX_id"/"$FILENAME""$FASTQ1_SUFFIX" "$TEMP_FULL_NAME_BAM"/"$SLX_id"/"$FILENAME""$FASTQ2_SUFFIX"
                                 

##STARTING NOVOALIGN PIPELINE

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Starting Novoalign alignment pipeline"         #Wrote in output (defined from sbatch command)

#Alignment Novoalign
"$NOVOALIGN_PATH"novoalign -d $TEMP_FULL_NAME_BAM/indexed.novoalign.human_g1k_v37_decoy_UCSC.mmu10.fa \
		           -f "$TEMP_FULL_NAME_BAM"/"$SLX_id"/"$FILENAME""$FASTQ1_SUFFIX" "$TEMP_FULL_NAME_BAM"/"$SLX_id"/"$FILENAME""$FASTQ2_SUFFIX" \
		           -o SAM "@RG\tID:"$FILENAME"\tPL:"$PLATFORM"\tPU:"$FL_BAR_LN"\tLB:"$SLX_id"\tSM:"$BARCODE"" \
		           -i PE 30-500 \
		           -t 250 \
		           -c $SLURM_TASKS_PER_NODE \
 		           -a CTGTCTCTTATA \
		           -H 20 \
                           -k \
		           -o FullNW \
		           2>$TEMP_FULL_NAME_BAM/$FILENAME.novostats.txt | \
                           "$SAMTOOLS_PATH"samtools view -1 -bS - > $TEMP_FULL_NAME_BAM/$FILENAME.bam
	# reference genome
	# input files
	# output format and @RG annotation
	# paired end option
	# absolute alignment score threshold
	# multithreading
	# adapter trimming
	# hard trimming
	# soft clipping
	# novoalign report

mv $TEMP_FULL_NAME_BAM/$FILENAME.novostats.txt $REPORTS_NOVOALIGN_ALIGNMENT         

rm -r "$TEMP_FULL_NAME_BAM"/"$SLX_id"

#Sort and index                                                                     #Sorts and indexes of bam file created by aligner
echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Running Novosort on Novoalign files"
"$NOVOALIGN_PATH"novosort  -m $MEMORY \
                           -t $TEMP_FULL_NAME_BAM \
                           -c $SLURM_TASKS_PER_NODE \
                           --markDuplicates \
                           --kt $TEMP_FULL_NAME_BAM/$FILENAME.bam \
                           -i \
                           -o $TEMP_FULL_NAME_BAM/$FILENAME.sorted.bam
rm $TEMP_FULL_NAME_BAM/$FILENAME.ba*

#Alignment stats                                                                   #Calculates alignment stats on reference genome
echo "[exome:bamstats $(date +"%Y-%m-%d %T")] Starting CollectAlignmentMetrics on Novoalign files"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                   CollectAlignmentSummaryMetrics \
                   R=$GENOME \
                   I=$TEMP_FULL_NAME_BAM/$FILENAME.sorted.bam \
                   O=$REPORTS_NOVOALIGN_ALIGNMENT/$FILENAME.alignMetrics.txt

#Bam stats                                                                         #Calculates BAM stats on reference genome and target regions
echo "[exome:bamstats $(date +"%Y-%m-%d %T")] Starting CollectHsMetrics on Novoalign files"
"$JAVA_PATH"java -Xmx$MEMORY -jar "$PICARD_PATH"picard.jar \
                   CollectHsMetrics \
                   BAIT_INTERVALS=$TARGET \
                   TARGET_INTERVALS=$TARGET \
                   I=$TEMP_FULL_NAME_BAM/$FILENAME.sorted.bam \
                   O=$REPORTS_NOVOALIGN_BAM_STATS/$FILENAME.bamMetrics.txt \
                   REFERENCE_SEQUENCE=$GENOME

mv $TEMP_FULL_NAME_BAM/$FILENAME.sorted.ba* $BAM_NOVO 

rm -r $TEMP_FULL_NAME_BAM

echo "[exome:fq_to_bam $(date +"%Y-%m-%d %T")] Alignment completed"
