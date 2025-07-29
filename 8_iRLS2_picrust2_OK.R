#R version 4.2.2 (2022-10-31)
#iRLS2 prev10 PICRUST Analysis
#Loading Data ####
rm(list=ls())
setwd("")
library("nlme")
library("reshape2")
library("dplyr")
library("phyloseq")
library("ggplot2")
library("seqinr")
library(biomformat)
library(compositions)

load("./240617_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq_PhyloTree.RData")

ps <- prune_samples(rowSums(otu_table(ps)) > 20000, ps)
str(ps)

ps <- subset_taxa(ps, !is.na(Phylum) & Kingdom %in% c("Archaea","Bacteria"))
table(tax_table(ps)[, "Kingdom"], exclude = NULL) 

# Compute prevalence of each feature (read)
# Prevalnce is the number of samples in which a taxa appears at least once

prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 24)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))
nrow(otu_table(ps))
ncol(otu_table(ps))
nrow(prevdf)

prevalenceThreshold =  10
#prevalenceThreshold =  0.10 * nsamples(ps)
#keep all 
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)


# # Abundance filtering
# ps_abundance_filter  = transform_sample_counts(ps2, function(x) x / sum(x) )
# ps_abundance_filter = filter_taxa(ps_abundance_filter, function(x) mean(x) > 1e-5, TRUE)
# keepTaxa = colnames(otu_table(ps_abundance_filter))
# ps2 <- prune_taxa(keepTaxa, ps2)
# str(ps2)

#filter for samples of interest
sampledata <- data.frame(sample_data(ps2))
sampledata_filt<-sampledata[sampledata$Group3%in%c("CTRL","iRLS","insonne"),]
sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, ps2)

#Analysis on ASV table
counts <- data.frame(t(otu_table(physeq_sub_filt)))

#Make ASV fasta ref
ID1 <- paste0("ID",c(1:nrow(counts)))
OTU_list <- rownames(counts) 
names(OTU_list) <- paste0(">",ID1)
OTU_list <- as.list(OTU_list)
write.fasta(sequences=OTU_list, file.out="iRLS2_ASV_sequence_char_picrust.fna", names= ID1)

#Make counts biom
table_f1 <- as.data.frame(counts)
rownames(table_f1) <- ID1
table_b1 <- make_biom(table_f1)
write_biom(table_b1, biom_file="iRLS2_ASV_table.biom")

#Write ASV table and samplesheet
write.table(counts,file="iRLS2_Metagenomics_OTU_Counts_prevfilt10_picrust.txt", sep="\t")
write.table(sampledata_filt,file="iRLS2_Metagenomics_Samplesheet_picrust.txt",sep="\t")

#Now perform analysis using picrust2 function from picrust2 environment
# picrust2_pipeline.py -s /standard/users/m.bacalini/METAGENOMICS/IRLS2_picrust2/iRLS2_ASV_sequence_char_picrust.fna -i /standard/users/m.bacalini/METAGENOMICS/IRLS2_picrust2/iRLS2_ASV_table.biom -o ./picrust2_output -p 10

#add_descriptions.py -i ./pred_metagenome_unstrat.tsv.gz -m EC -o /standard/users/m.bacalini/METAGENOMICS/IRLS2_picrust2/picrust2_output/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
#add_descriptions.py -i ./path_abun_unstrat.tsv.gz -m METACYC -o /standard/users/m.bacalini/METAGENOMICS/IRLS2_picrust2/picrust2_output/pathways_out/path_abun_unstrat_descrip.tsv.gz

#pathway_Differential Analysis####
rm(list=ls())
library(data.table)
setwd("./IRLS2_picrust2/")
p_AbsAbundance <- data.frame(fread("./path_abun_unstrat_descrip.tsv", header=T,sep="\t"))
rownames(p_AbsAbundance) <- p_AbsAbundance$pathway
p_AbsAbundance <- p_AbsAbundance[,-1]
p_AbsAbundance <- droplevels(p_AbsAbundance)
head(p_AbsAbundance)

p_RelAbundance <-p_AbsAbundance[,-1]
p_RelAbundance<-clr(p_RelAbundance+1)
#for (i in 2:ncol(p_AbsAbundance)){
#  p_RelAbundance[,i] <- p_AbsAbundance[,i]/sum(p_AbsAbundance[,i],na.rm=T)
#}

##Samplesheet and Sample Filter####
ss <- read.table("./iRLS2_Metagenomics_Samplesheet_picrust.txt", sep="\t", header=T)
ss_filt <- ss[ss$Group3%in%c("CTRL","iRLS","insonne"),]
ss_filt <- droplevels(ss_filt)

p_RelAbundance_filt <- p_RelAbundance[,colnames(p_RelAbundance)%in%ss_filt$SampleID]
dim(p_RelAbundance_filt)

ss_filt <- ss_filt[ss_filt$SampleID%in%colnames(p_RelAbundance_filt),]
dim(ss_filt)

table(ss_filt$SampleID==colnames(p_RelAbundance_filt)) #sono in ordine

##Analysis####
library(limma)

Group <- factor(ss_filt$Group3,levels=c("CTRL","iRLS","insonne"))
Sesso <- factor(ss_filt$Sesso)
Age <- ss_filt$Age
BMI <- ss_filt$BMI1
Psichiatrico <- factor(ss_filt$Psichiatrico)
Batch <-factor(ss_filt$Bach)

design <- model.matrix(~0+Group+Sesso+Age+BMI+Psichiatrico+Batch, data=p_RelAbundance_filt)
head(design)
fit <- lmFit(p_RelAbundance_filt,design=design)
contrast.matrix <- makeContrasts(GroupiRLS-GroupCTRL,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl1 <- topTable(fit2, coef="GroupiRLS - GroupCTRL", num=Inf, sort.by="p")
head(topl1)
topl1 <- data.frame(rownames(topl1),topl1)
topf<- merge(p_AbsAbundance[,"description",drop=F],topl1,by="row.names")
topf<-topf[order(topf$adj.P.Val),]
topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
#write.table(topf,file="20250723_iRLSvsCTR_Limma_pathway.txt", sep="\t")

#iRLS vs insonne
contrast.matrix <- makeContrasts(GroupiRLS-Groupinsonne,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl2 <- topTable(fit2, coef="GroupiRLS - Groupinsonne", num=Inf, sort.by="p")
head(topl2)
topl2 <- data.frame(rownames(topl2),topl2)
topf2<- merge(p_AbsAbundance[,"description",drop=F],topl2,by="row.names")
topf<-merge(topf,topf2,by="Row.names")
#topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
#write.table(topf2,file="20250723_iRLSvsinsonne_Limma_pathway.txt", sep="\t")

#insonne vs iRLS
contrast.matrix <- makeContrasts(Groupinsonne-GroupCTRL,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl3 <- topTable(fit2, coef="Groupinsonne - GroupCTRL", num=Inf, sort.by="p")
head(topl3)
topl3 <- data.frame(rownames(topl3),topl3)
topf3<- merge(p_AbsAbundance[,"description",drop=F],topl3,by="row.names")
topf<-merge(topf,topf3,by="Row.names")
#topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
write.table(topf,file="20250723_ALL_Limma_pathway_.txt", sep="\t")

#EC_Differential Analysis####
rm(list=ls())
setwd("./IRLS2_picrust2/")
EC_AbsAbundance <- data.frame(fread("./pred_metagenome_unstrat_descrip.tsv", header=T,sep="\t"))
rownames(EC_AbsAbundance) <- EC_AbsAbundance$function.
EC_AbsAbundance <- EC_AbsAbundance[,-1]
EC_AbsAbundance <- droplevels(EC_AbsAbundance)
head(EC_AbsAbundance)

EC_RelAbundance <-EC_AbsAbundance[,-1]
EC_RelAbundance<-clr(EC_RelAbundance+1)
#for (i in 2:ncol(EC_AbsAbundance)){
#  EC_RelAbundance[,i] <- EC_AbsAbundance[,i]/sum(EC_AbsAbundance[,i],na.rm=T)
#}

#colSums(EC_AbsAbundance)

##Samplesheet and Sample Filter####
ss <- read.table("./iRLS2_Metagenomics_Samplesheet_picrust.txt", sep="\t", header=T)
ss_filt <- ss[ss$Group3%in%c("CTRL","iRLS","insonne"),]
ss_filt <- droplevels(ss_filt)

EC_RelAbundance_filt <- EC_RelAbundance[,colnames(EC_RelAbundance)%in%ss_filt$SampleID]
dim(EC_RelAbundance_filt)

ss_filt <- ss_filt[ss_filt$SampleID%in%colnames(EC_RelAbundance_filt),]
dim(ss_filt)

table(ss_filt$SampleID==colnames(EC_RelAbundance_filt)) #sono in ordine

##Analysis####
library(limma)

Group <- factor(ss_filt$Group3,levels=c("CTRL","iRLS","insonne"))
Sesso <- factor(ss_filt$Sesso)
Age <- ss_filt$Age
BMI <- ss_filt$BMI1
Psichiatrico <- factor(ss_filt$Psichiatrico)
Batch <-factor(ss_filt$Bach)

design <- model.matrix(~0+Group+Sesso+Age+BMI+Psichiatrico+Batch, data=EC_RelAbundance_filt)
head(design)
fit <- lmFit(EC_RelAbundance_filt,design=design)
contrast.matrix <- makeContrasts(GroupiRLS-GroupCTRL,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl1 <- topTable(fit2, coef="GroupiRLS - GroupCTRL", num=Inf, sort.by="p")
head(topl1)
topl1 <- data.frame(rownames(topl1),topl1)
topf<-merge(EC_AbsAbundance[,"description",drop=F],topl1,by="row.names")
topf<-topf[order(topf$adj.P.Val),]
#topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
#write.table(topf,file="20250723_iRLSvsCTR_Limma_EC.txt", sep="\t")

#iRLS vs insonne
contrast.matrix <- makeContrasts(GroupiRLS-Groupinsonne,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl2 <- topTable(fit2, coef="GroupiRLS - Groupinsonne", num=Inf, sort.by="p")
head(topl2)
topl2 <- data.frame(rownames(topl2),topl2)
topf2<-merge(EC_AbsAbundance[,"description",drop=F],topl2,by="row.names")
topf<-merge(topf,topf2,by="Row.names")
#topf<-topf[order(topf$adj.P.Val),]
#topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
#write.table(topf,file="20250723_iRLSvsinsonne_Limma_EC.txt", sep="\t")

#insonne vs iRLS
contrast.matrix <- makeContrasts(Groupinsonne-GroupCTRL,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topl3 <- topTable(fit2, coef="Groupinsonne - GroupCTRL", num=Inf, sort.by="p")
head(topl3)
topl3 <- data.frame(rownames(topl3),topl3)
topf3<-merge(EC_AbsAbundance[,"description",drop=F],topl3,by="row.names")
topf<-merge(topf,topf3,by="Row.names")
#topf<-topf[order(topf$adj.P.Val),]
#topf_sig <- topf[topf$P.Value<0.05&!is.na(topf$P.Value),]
write.table(topf,file="20250723_ALL_Limma_EC.txt", sep="\t")
