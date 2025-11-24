#Metagenomics Alternative Differential Analysis####
#R version 4.4.0 (2024-04-24)
#ALDEx2

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ALDEx2")

#Loading Data ####
rm(list=ls())
setwd("")
library("ALDEx2")
library("nlme")
library("reshape2")
library("dplyr")
library("phyloseq")
library("ggplot2")

load("./240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq.RData")

#Filtro per profondità di lettura
ps <- prune_samples(rowSums(otu_table(ps)) > 20000, ps)
# abbiamo tolto i campioni 1514431congelatoDNAestrattodaloro, 20274congelatoDNAestrattodaloro, bianco, IRLS024
#PelletbufferRLT1514431estrattodanoi, Pelletcong1514431estrattodanoi; rimangono 90 campioni
str(ps)

# Filtro per annotazione a regno batetri Remove NA Phyla and keep only Bacteria and Archaea
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

#GENUS Analysis ####
ps_GENUS <- tax_glom(physeq_sub_filt,"Genus",NArm = T)
## Use absolute counts ####
counts <- data.frame(t(otu_table(ps_GENUS)))


#iRLS VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)
  
## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command GLM####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_GENUS))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_GENUS_iRLSvCTRL_Aldex2.txt",sep="\t",row.names=F)

#iRLS VS insonne ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","insonne"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation  ####
taxa <- data.frame(tax_table(ps_GENUS))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_B <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_GENUS_iRLSvinsonne_Aldex2.txt",sep="\t",row.names=F)

#insonne VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("insonne","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_GENUS))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_C <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_GENUS_insonnevCTRL_Aldex2.txt",sep="\t",row.names=F)

#FAMILY Analysis ####

ps_FAMILY <- tax_glom(physeq_sub_filt,"Family",NArm = T)
## Use absolute counts ####
counts <- data.frame(t(otu_table(ps_FAMILY)))


#iRLS VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command GLM####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_FAMILY))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_FAMILY_iRLSvCTRL_Aldex2.txt",sep="\t",row.names=F)

#iRLS VS insonne ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","insonne"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_FAMILY))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_B <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_FAMILY_iRLSvinsonne_Aldex2.txt",sep="\t",row.names=F)

#insonne VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("insonne","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_FAMILY))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_C <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_FAMILY_insonnevCTRL_Aldex2.txt",sep="\t",row.names=F)

#PHYLUM Analysis####
ps_PHYLUM <- tax_glom(physeq_sub_filt,"Phylum",NArm = T)
## Use absolute counts ####
counts <- data.frame(t(otu_table(ps_PHYLUM)))


#iRLS VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command GLM####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_PHYLUM))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_PHYLUM_iRLSvCTRL_Aldex2.txt",sep="\t",row.names=F)

#iRLS VS insonne ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("iRLS","insonne"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_PHYLUM))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_B <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_PHYLUM_iRLSvinsonne_Aldex2.txt",sep="\t",row.names=F)

#insonne VS CTRL ####
sampledata_temp<-sampledata_filt[sampledata_filt$Group3%in%c("insonne","CTRL"),]
counts_temp <- counts[,colnames(counts)%in%sampledata_temp$SampleID]
table(colnames(counts_temp)==sampledata_temp$SampleID)

## Set model matrix ####
mm_temp <- model.matrix(~ Group3 + Age + Sesso + BMI1 + Psichiatrico + Bach, sampledata_temp)

## Launch command ####
x.glm <- aldex.clr(counts_temp, mm_temp, mc.samples=128, denom="all", verbose=T)
glm.test <- aldex.glm(x.glm, mm_temp, fdr.method='BH')
glm.eff<- aldex.glm.effect(x.glm)
glm.eff.df <- data.frame(glm.eff)

## Results Annotation ####
taxa <- data.frame(tax_table(ps_PHYLUM))
final <- merge(taxa,glm.test, by="row.names")
rownames(final) <- final$Row.names
final <- final[,-1]
final_effect <- merge(final,glm.eff.df,by="row.names")

sig_final_C <- final_effect[final_effect$`Group3iRLS:pval`<0.05,]
write.table(final_effect,file="iRLS2_Metagenomics_PHYLUM_insonnevCTRL_Aldex2.txt",sep="\t",row.names=F)

