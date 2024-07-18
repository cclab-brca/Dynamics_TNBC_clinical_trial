
rm(list=ls()) # clearing del workspace
graphics.off()

library(VariantAnnotation)
args <- commandArgs()

annovar_vcf_dir <- args[2]
sample_name <- args[3]
output_dir <- args[4]
PON_path <- args[5]

## CONVERT VCF INTO HUMAN READABLE:
# READ THE VCF as table:
samples <- sample_name  
vcf_input <- paste0(annovar_vcf_dir, "/", sample_name, ".annovar.hg19_multianno.vcf")   
param <- ScanVcfParam()
vcf_combined_as_table <-tryCatch(read.table(vcf_input, comment.char = "#", header=F, stringsAsFactors=F), error=function(e) NULL)
colnames(vcf_combined_as_table) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_name)

# DO THIS ONLY IF THERE ARE MULTIALLELIC SITES:
if (length(grep(',', vcf_combined_as_table$ALT))>0) { #if there is at least 1 multiallelic site
  
  ## MANIPULATE VCF IN ORDER TO 'UNPACK' THE MULTIALLELIC SITES: 
  vcf_multiallelic <- vcf_combined_as_table[grep(',', vcf_combined_as_table$ALT),]
  vcf_single_allele <- vcf_combined_as_table[grep(',', vcf_combined_as_table$ALT, invert = T),]
  nrow(vcf_multiallelic) + nrow(vcf_single_allele) == nrow(vcf_combined_as_table)
  
  # EXPAND MULTIALLELIC VCF INTO A SINGLE MUTATION PER ROW:
  vcf_multiallelic_expanded <- data.frame()
  for (i in 1:nrow(vcf_multiallelic)) {       #for each row of multiallelic
    vcf_row <- vcf_multiallelic[i,]
    nalleles <- length(gregexpr(',', vcf_row$ALT )[[1]]) + 1          #determine the number of alleles of this multiallelic site (they should always be 2, but just in case it's more) #grep does not work, I had to use gregexpr. (the [[1]] is because it returns an list of 1 instead of a simple value)
    ALT_per_allele <- unlist(apply(as.matrix(vcf_row$ALT),1,function(x){unlist(strsplit(x,","))}))           #ALT_per_allele= ALT of each allele
    info_per_allele <- unlist(apply(as.matrix(vcf_row$INFO),1,function(x){unlist(strsplit(x,"ANNOVAR_DATE="))}))[-1]  #V8 IS INFO COLUMN OF THE VCF: INFO for each allele
    FORMAT_per_allele <- character()
    for (k in 1:length(samples)){                           #For each sample:
      geno<-strsplit(as.character(vcf_row[9+k]),":")[[1]]   #geno= FORMAT of that sample
      FORMAT_per_allele_of_sample <- character()            #FORMAT_per_allele_of_sample= FORMAT of that sample, expanded per allele (one row per allele)
      for (j in 1:nalleles) {
        ALT_allele <- paste(strsplit(geno[2], ",")[[1]][c(1,j+1)], collapse = ",")    #ALT_allele= REF,ALT_j. (geno[2] is AD, which is in the format "REF,ALT1,ALT2,ALT3,..." | [[1]] is because strsplit produces a list of length 1, [1=REF, j=ALT_j] | collapse (join them with , as separator))
        geno_per_allele <- geno                                                      
        geno_per_allele[2] <- ALT_allele                    #geno_per_alllele= FORMAT information of that specific allele
        FORMAT_per_allele_of_sample <- rbind(FORMAT_per_allele_of_sample, as.character(paste(geno_per_allele, collapse = ":"))) 
      }
      colnames(FORMAT_per_allele_of_sample) <- samples[k]   #add the name of the sample
      FORMAT_per_allele <- cbind(FORMAT_per_allele, FORMAT_per_allele_of_sample)
    }
    expanded_vcf_row <- cbind(data.frame(CHROM=rep(vcf_row$CHROM, nalleles), POS=rep(vcf_row$POS, nalleles), ID=rep(vcf_row$ID, nalleles), REF=rep(vcf_row$REF, nalleles)),
                              ALT=ALT_per_allele,
                              data.frame(QUAL=rep(vcf_row$QUAL, nalleles), FILTER=rep(vcf_row$FILTER, nalleles)), 
                              INFO=info_per_allele,
                              data.frame(FORMAT= rep(vcf_row$FORMAT, nalleles)), 
                              FORMAT_per_allele)
    vcf_multiallelic_expanded <- rbind(vcf_multiallelic_expanded, expanded_vcf_row) 
  }
  
  # JOIN EXPANDED_VCF WITH VCF_SINGLE_ALLELE AND SORT:
  vcf_combined_expanded <- rbind(vcf_single_allele, vcf_multiallelic_expanded)
} else {
  vcf_combined_expanded <- vcf_combined_as_table
}

#vcf_combined_expanded_sorted <- vcf_combined_expanded [with(vcf_combined_expanded, order(CHROM, POS)),]
ChrOrder <- unique(vcf_combined_as_table$CHROM)
vcf_combined_expanded_sorted <- data.frame()
for (i in 1:length(ChrOrder)) {
  chr_subset <- vcf_combined_expanded[vcf_combined_expanded$CHROM==ChrOrder[i],]
  chr_subset_sorted <- chr_subset[order(chr_subset$POS), ]
  vcf_combined_expanded_sorted <- rbind(vcf_combined_expanded_sorted, chr_subset_sorted)
}

# EXTRACT CHR POS REF ALT:
mut<-vcf_combined_expanded_sorted[,c(1,2,4:5)]
colnames(mut)<-c("Chr","Position","REF","ALT")
head(mut)
mutID<-paste(mut$Chr,mut$Position,mut$REF,mut$ALT,sep="_")

# EXTRACT MUTATION INFO
annot<-unlist(apply(as.matrix(vcf_combined_expanded_sorted$INFO),1,function(x){unlist(strsplit(x,";"))}))  
Func.refGene<-gsub("Func.refGene=","", grep("^Func.refGene=",annot,value=T))
Gene.refGene<-gsub("Gene.refGene=","", grep("^Gene.refGene=",annot,value=T))
GeneDetail.refGene<-gsub("GeneDetail.refGene=","", grep("^GeneDetail.refGene=",annot,value=T))
ExonicFunc.refGene<-gsub("ExonicFunc.refGene=","", grep("^ExonicFunc.refGene=",annot,value=T))
AAChange.refGene<-gsub("AAChange.refGene=","", grep("^AAChange.refGene=",annot,value=T))
cytoBand<-gsub("cytoBand=","", grep("^cytoBand=",annot,value=T))
all.1000g2015aug<-gsub("1000g2015aug_all=","", grep("^1000g2015aug_all=",annot,value=T))
snp142<-gsub("snp142=","", grep("^snp142=",annot,value=T))
SIFT_score<-gsub("SIFT_score=","", grep("^SIFT_score=",annot,value=T))
SIFT_pred<-gsub("SIFT_pred=","", grep("^SIFT_pred=",annot,value=T))
Polyphen2_HDIV_score<-gsub("Polyphen2_HDIV_score=","", grep("^Polyphen2_HDIV_score=",annot,value=T))
Polyphen2_HDIV_pred<-gsub("Polyphen2_HDIV_pred=","", grep("^Polyphen2_HDIV_pred=",annot,value=T))
Polyphen2_HVAR_score<-gsub("Polyphen2_HVAR_score=","", grep("^Polyphen2_HVAR_score=",annot,value=T))
Polyphen2_HVAR_pred<-gsub("Polyphen2_HVAR_pred=","", grep("^Polyphen2_HVAR_pred=",annot,value=T))
CLINSIG<-gsub("CLINSIG=","", grep("^CLINSIG=",annot,value=T))
CLNDBN<-gsub("CLNDBN=","", grep("^CLNDBN=",annot,value=T))
CLNACC <- gsub("CLNACC=","", grep("^CLNACC=",annot,value=T))
CLNDSDB <- gsub("CLNDSDB=","", grep("^CLNDSDB=",annot,value=T))
CLNDSDBID <- gsub("CLNDSDBID=","", grep("^CLNDSDBID=",annot,value=T))
ICGC_Id <- gsub("ICGC_Id=","", grep("^ICGC_Id=",annot,value=T))
ICGC_Occurrence <- gsub("ICGC_Occurrence=","", grep("^ICGC_Occurrence=",annot,value=T))
COSMICID_coding <- gsub("COSMICID_coding=","", grep("^COSMICID_coding=",annot,value=T))
COSMICOCCURENCE_coding <- gsub("COSMICOCCURENCE_coding=","", grep("^COSMICOCCURENCE_coding=",annot,value=T))
COSMICID_noncoding <- gsub("COSMICID_noncoding=","", grep("^COSMICID_noncoding=",annot,value=T))
COSMICOCCURENCE_noncoding <- gsub("COSMICOCCURENCE_noncoding=","", grep("^COSMICOCCURENCE_noncoding=",annot,value=T))

#ADD REF.COUNT, ALT.COUNT, VAF TO EACH SAMPLE:
sample_count_VAF_all<-character()
grange_combined_calls <- makeGRangesFromDataFrame(mut, seqnames.field = "Chr", start.field = "Position", end.field = "Position")
for(k in 1:length(samples)) {                           #for every sample:
  #COMPUTE REF.COUNT, ALT.COUNT, VAF:
  vcf.k<-as.matrix(vcf_combined_expanded_sorted[,9+k])      #samples start from column 10
  #vcf.k[vcf.k=="./."]<-"0/0:0,0:0:0:0,0,0"                 #replace empty cells with 0s, this was mau's comment but it does not work with the current grep
  vcf.k[grep('./.', vcf.k, fixed = T)]<-"0/0:0,0:0:0:0,0,0" #this one does, you have to set fixed=T to tell grep to not consider ./. as special characters
  geno<-apply(vcf.k,1,function(x){unlist(strsplit(x,":"))})
  AD<-geno[2,]
  cov<-t(apply(as.matrix(AD),1,function(x){strsplit(x,",")[[1]][1:2]}))
  class(cov)<-"numeric"
  VAF<-cov[,2]/(cov[,1]+cov[,2])	
  
  #Put all info together
  sample_count_VAF <- cbind(cov,VAF)
  colnames(sample_count_VAF)<-paste(samples[k],c("ref.counts","alt.counts","VAF"),sep=".")
  sample_count_VAF_all <- cbind(sample_count_VAF_all, sample_count_VAF)
}

## ADD PON COLUMN:
param <- ScanVcfParam()
PON_vcf <- readVcf(PON_path, "hg19", param)
pon_fo <- findOverlaps(grange_combined_calls, rowRanges(PON_vcf))
PON <- rep("NOT_IN_PON", length(grange_combined_calls))
PON[queryHits(pon_fo)] <- "PON"

data<-cbind(mutID,mut,PON, sample_count_VAF_all,
            Func.refGene,Gene.refGene,GeneDetail.refGene,ExonicFunc.refGene,AAChange.refGene,cytoBand,all.1000g2015aug,snp142,
            SIFT_score,SIFT_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,
            CLINSIG,CLNDBN,CLNACC,CLNDSDB,CLNDSDBID,ICGC_Id,ICGC_Occurrence,COSMICID_coding,COSMICOCCURENCE_coding,COSMICOCCURENCE_coding,COSMICID_noncoding)

write.table(data, paste0(output_dir, "/" ,sample_name, ".finalcalls.txt"), sep = "\t", quote = F, row.names = F)




