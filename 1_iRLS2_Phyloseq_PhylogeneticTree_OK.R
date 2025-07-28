# Phyloseq+Phylogenetic Tree ####
rm(list=ls())
suppressMessages(library("knitr"))
suppressMessages(library("BiocStyle"))
suppressMessages(library("ggplot2"))
suppressMessages(library(gridExtra))
suppressMessages(library(dada2))
suppressMessages(library(phyloseq))
suppressMessages(library(DECIPHER))
suppressMessages(library(phangorn))
suppressMessages(library(RcppParallel))
suppressMessages(library(grid))
suppressMessages(library(lattice))
theme_set(theme_bw())

#Loading Station ####
seqtab <- readRDS("./240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_Assign_Species_seqtab.rds")
tax.plus <- readRDS("./240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_Assign_Species_tax.rds")

taxTab <- tax.plus
write.table(taxTab,file=paste0(output_folder, output_file,"_taxTab.txt"),row.names=T,sep="\t")
seqtabNoC <- seqtab
dim(seqtabNoC)
rm(tax)
rm(tax.plus)
rm(seqtab)

table(is.na(taxTab[,1]))
table(taxTab[,1])

samdf <- read.table("./SampleSheet_IRLS_Full_BigData.txt",sep="\t",header=TRUE,quote="")
str(samdf)
dim(samdf)
colnames(samdf)[1] <- "SampleID"
rownames(samdf) <- samdf$SampleID
table(rownames(seqtabNoC)==samdf$SampleID)
samdf<-samdf[match(rownames(seqtabNoC),samdf$SampleID),]
table(rownames(seqtabNoC)==samdf$SampleID)
samdf$Group3<-factor(samdf$Group3)
samdf$Group4<-factor(samdf$Group4)

# Costruct Phylogenetic Tree ####
seqs<-getSequences(seqtabNoC)
names(seqs)<-seqs
alignment<-AlignSeqs(DNAStringSet(seqs, anchor=NA))

phang.align<-phyDat(as(alignment,"matrix"),type="DNA")
dm<-dist.ml(phang.align)
treeNJ<-NJ(dm)
fit<-pml(treeNJ, data=phang.align)

##negative edge length to 0

fitGTR<-update(fit, k=4, inv=0.2)
fitGTR<-optim.pml(fitGTR,model="GTR",optInv=T,optGamma=T,rearrangement="NNI", control=pml.control(trace=0))
saveRDS(fitGTR, "240618_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_PhyloTree.rds")

fitGTR<-readRDS("240618_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_PhyloTree.rds")

# Create phyloseq object ####
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxTab),phy_tree(fitGTR$tree))
save(ps, file="240617_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq_PhyloTree.RData")
str(ps)
