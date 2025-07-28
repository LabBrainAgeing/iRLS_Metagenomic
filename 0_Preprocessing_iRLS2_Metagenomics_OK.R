#DADA2 Pipeline Preprocess####
library(dada2) # v1.26
setwd("")
exp <- "./240513_MMISEQ_Metagenomica_iRLS_iSNSIBO"

inpath <- ""
#command <- paste0("cd ",inpath)
#cat(paste0("Type in your terminal: ", command))
#command <- paste0("gunzip *.fastq.gz")
#cat(paste0("Type in your terminal: ", command))
list.files(inpath)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(inpath, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(inpath, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##Inspect Quality ####
dir.create("./FastQC")
pdf(paste0("./FastQC/",exp,"_QualityCheck_fastq.pdf"))
for (i in 1:length(fnFs)) {
print(plotQualityProfile(fnFs[i]))
print(plotQualityProfile(fnRs[i]))
}
dev.off()

## Filter and Trim ####
dir.create(paste0(inpath,"/filtered"))
filtFs <- file.path(inpath, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(inpath, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(10,10), truncLen=c(260,215),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

## Learn Error Rate ####
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf(paste0(exp,"_errorRates.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ = TRUE)
dev.off()

##Dereplication ####
#Incluso nello step di Inference dada2 v1.14
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

## Learn Error Rate ####
errF <- learnErrors(derepFs, multithread=TRUE)
errR <- learnErrors(derepRs, multithread=TRUE)

pdf(paste0(exp,"_errorRates_derep.pdf"))
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

##Sample Inference #### 

#dadaFs <- dada(derepFs, err=errF, multithread = TRUE) v1.4
#dadaRs <- dada(derepRs, err=errR, multithread = TRUE) v1.4
dadaFs <- dada(derepFs, err=errF, multithread = TRUE, pool=TRUE, verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread = TRUE, pool=TRUE, verbose=TRUE)

## PairEnd Merging ####
#mergers <- mergePairs(dadaFs,derepFs,dadaRs,derepRs, verbose=TRUE) #v1.4
mergers <- mergePairs(dadaFs,filtFs,dadaRs,filtRs, verbose=TRUE) #v1.14.1

## Construct Sequence Table ####
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab))) # distribution amplicon length - which is the expected amplicon length?
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% #####]
saveRDS(seqtab, file=paste0(exp,"_seqtab.rds"))

### Remove Chimeras ####
seqtab_noC <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
saveRDS(seqtab_noC, file=paste0(exp,"_seqtab_nochim.rds"))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(derepFs,getN),sapply(derepRs,getN),sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab_noC))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered","derepF","derepR","denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
write.table(track, file=paste0(exp,"_readCount.txt"), sep="\t")

# Assign Taxonomy ####
#taxa <- assignTaxonomy(seqtab_noC, "/standard/users/m.bacalini/METAGENOMICS/IRLS_Analysis/rdp_train_set_18.fa.gz", multithread=TRUE)
taxa <- assignTaxonomy(seqtab_noC, "./silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, file=paste0(exp,"_tax.rds"))

taxa_df<-data.frame(taxa, stringsAsFactors = T)
rownames(taxa_df) <- NULL
summary(taxa_df)
#write.table(summary(taxa_df),file=paste0(exp,"_summary_TaxAssign_RDP.txt"), sep="\t",quote=F)
write.table(summary(taxa_df),file=paste0(exp,"_summary_TaxAssign_silva.txt"), sep="\t",quote=F)
