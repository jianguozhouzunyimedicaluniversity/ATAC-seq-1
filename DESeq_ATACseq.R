source("https://bioconductor.org/biocLite.R")

#biocLite("ChIPseeker")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(ChIPseeker)
library(org.Hs.eg.db)
library(ChIPQC)
library(soGGi)
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(limma)
library(Rsubread)
library(DESeq2)
#library(clusterProfiler)
library(rGREAT)
library(Rsamtools)
library(ggplot2)
library(rtracklayer)
library(Biostrings)
library(GenomicAlignments)
library(tracktables)

setwd("/athena/khuranalab/scratch/kac2053/projects/ctDNA/plots/")

# Get narrowpeak files.
sample_dir_fanying <- "/athena/khuranalab/scratch/kac2053/projects/ctDNA/Fanying_NEPC_CRPC_ATAC-seq/macs2"
sample_dir_park <- "/athena/khuranalab/scratch/kac2053/projects/ctDNA/GSE118207/macs2"
sample_dir_GM12878 <- "/athena/khuranalab/scratch/kac2053/projects/ctDNA/GSE47753/macs2"

peaks_fanying <- dir(sample_dir_fanying, pattern = "*header.narrowPeak", full.names = TRUE)
peaks_park <- dir(sample_dir_park, pattern = "*header.narrowPeak", full.names = TRUE)
peaks_GM12878 <- dir(sample_dir_GM12878, pattern = "*header.narrowPeak", full.names = TRUE)

peaks <- c(peaks_fanying, peaks_park, peaks_GM12878)

# Obtain sample names.
peak.names <- strsplit(peaks, "/")
peak.names <- lapply(peak.names, tail, n = 1)
peak.names <- unlist(peak.names)
peak.names <- strsplit(peak.names,"_")
peak.names <- lapply(peak.names, head, n = 1)
peak.names <- unlist(peak.names)

# Get GRanges in narrowpeak files and merge all GRanges as set of non-redundant (no overlapping regions) open regions using soGGi.
myPeaks <- lapply(peaks,ChIPQC:::GetGRanges,simple=TRUE)
names(myPeaks) <- peak.names
Group <- factor(c("CRPC", "CRPC", "NEPC", "NEPC", "CRPC", "CRPC", "NEPC", "NEPC", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid"))
prostatevslymph <- factor(c("prostate", "prostate", "prostate", "prostate", "prostate", "prostate", "prostate", "prostate", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid", "lymphoblastoid"))
consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks),"none")

###
library(tidyr)

savePlot <- function(myPlot, plotname) {
  pdf(paste0(plotname, ".pdf"))
  print(myPlot)
  dev.off()
}

pca <- as.data.frame(elementMetadata(consensusToCount)) %>% 
  dplyr::select(-consensusIDs) %>% 
  as.matrix %>% t %>% prcomp %>% .$x %>% data.frame %>% 
  mutate(Samples=rownames(.)) %>% 
  mutate(Group=gsub("_\\d","",Samples)) %>% 
  ggplot(aes(x=PC1,y=PC2,colour=Group))+geom_point(size=5)

savePlot(pca, "pca_occurence")

###
library(Rsubread)

# Count the number of occurences in each genomic range and store as a vector.
occurrences <- elementMetadata(consensusToCount) %>% as.data.frame %>% dplyr::select(-consensusIDs) %>% rowSums
table(occurrences) %>% rev %>% cumsum

# Only save the genomic ranges that occur in 2+ samples.
consensusToCount <- consensusToCount[occurrences >= 2,]


###
bam_fanying <- dir("/athena/khuranalab/scratch/kac2053/projects/ctDNA/Fanying_NEPC_CRPC_ATAC-seq/original_data", full.names = TRUE,pattern = "*.\\.bam$")
bam_park <- dir("/athena/khuranalab/scratch/kac2053/projects/ctDNA/GSE118207/bam", full.names = TRUE,pattern = "*_sorted.bam$")
bam_GM12878 <- dir("/athena/khuranalab/scratch/kac2053/projects/ctDNA/GSE47753/bam", full.names = TRUE,pattern = "*_sorted.bam$")

bamsToCount <- c(bam_fanying, bam_park, bam_GM12878)
#indexBam(bamsToCount)
regionsToCount <- data.frame(GeneID=paste("ID",seqnames(consensusToCount),start(consensusToCount),end(consensusToCount),sep="_"),Chr=seqnames(consensusToCount),Start=start(consensusToCount),End=end(consensusToCount),Strand=strand(consensusToCount))

fcResults <- featureCounts(bamsToCount,annot.ext=regionsToCount,isPairedEnd = TRUE,countMultiMappingReads = FALSE,maxFragLength=100)
myCounts <- fcResults$counts
#save(myCounts,file=paste0("countsFromATAC.RData"))
colnames(myCounts) <- c("C4.2_1","C4.2_2","H660_1","H660_2","VCaP_1", "VCaP_2", "MSKCC-EF1", "H660_Park", "GM12878_50k_1", "GM12878_50k_2", "GM12878_50k_3", "GM12878_50k_4", "GM12878_500_1", "GM12878_500_2", "GM12878_500_3")
save(myCounts,file=paste0("countsFromATAC.RData"))


###
library(DESeq2)
library(ggrepel)
load(paste0("countsFromATAC.RData"))
metaData <- data.frame(Group,row.names=colnames(myCounts))
atacDDS <- DESeqDataSetFromMatrix(myCounts,metaData,~Group,rowRanges=consensusToCount)
atacDDS <- DESeq(atacDDS) #view atacDDS with > counts(atacDDS, normalized=TRUE)
atac_Rlog <- rlog(atacDDS) 
save(atacDDS,file=paste0("atacDDS_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))
save(atac_Rlog,file=paste0("atac_Rlog_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))
pca_signal <- plotPCA(atac_Rlog,intgroup="Group",ntop=nrow(atac_Rlog)) + geom_label_repel(label = rownames(colData(atac_Rlog)), nudge_y=50, nudge_x=50)
savePlot(pca_signal, "pca_signal")
#save(atacDDS,file=paste0("atacDDS_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))
#save(atac_Rlog,file=paste0("atac_Rlog_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))


###
library(ComplexHeatmap)

load(paste0("atacDDS_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))
load(paste0("atac_Rlog_Park_Fanying_GM12878_CRPCvsNEPCvslymphoblastoid.RData"))

atac_heatmap <- counts(atacDDS, normalized=TRUE)
atac_var_heatmap <- apply(atac_heatmap,1,var)
var_regions <- sort(atac_var_heatmap, decreasing = TRUE)
var_regions <- var_regions[1:1000] #Top 1000 variance regions.

# Subset atac-seq data to contain only top 1000 variance regions.
atac_1000var_heatmap <- subset(atac_heatmap, rownames(atac_heatmap) %in% names(var_regions) )

atac_log_heatmap <- log2(atac_1000var_heatmap+1)

pdf(paste0("differential_accessibility_CRPCvsNEPCvslymphoblastoid.pdf"))
Heatmap(atac_log_heatmap,column_title = "top 1000 variant regions", name = "normalized peak signal", show_row_names = FALSE)
dev.off()
#savePlot(atac_heatmap, "diffaccess_heatmap_prostatevslymph")

###
library(tracktables)
CRPCMinusNEPC <- results(atacDDS,c("Group","CRPC","NEPC"),format="GRanges")
head(CRPCMinusNEPC)
CRPCMinusNEPC <- CRPCMinusNEPC[order(CRPCMinusNEPC$pvalue)]
head(CRPCMinusNEPC)




