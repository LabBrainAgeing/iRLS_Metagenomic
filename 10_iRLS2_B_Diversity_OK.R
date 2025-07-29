rm(list=ls())
setwd("")
library("nlme")
library("reshape2")
library("dplyr")
library("phyloseq")
library("ggplot2")
library("vegan")

load("./240617_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq_PhyloTree.RData")


ps <- prune_samples(rowSums(otu_table(ps)) > 20000, ps)

str(ps)

# Group_missing <- sample_data(ps)[is.na(sample_data(ps)$Group)]
# Group_missing
# ps.temp <- prune_samples(!is.na(sample_data(ps)$Group), ps)
# str(ps.temp)
# setdiff(rownames(sample_data(ps)),rownames(sample_data(ps.temp)))

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
sampledata <- data.frame(sample_data(ps))
sampledata_filt<-sampledata[sampledata$Group3%in%c("CTRL","iRLS","insonne"),]
sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, ps2)


# BETA DIVERSITY
# - PCoA on Bray-Curtis dissimilarity
# - DPCoA
# - PCoA on  Weighted/Unweighted Unifrac distance
# - PERMANOVA Analysis of Group distances
###################################
# PCoA with Bray Curtis
pslog <- transform_sample_counts(physeq_sub_filt, function(x) log(1 + x))
pslog_sampledata<-data.frame(sample_data(pslog))


out.pcoa.log <- ordinate(pslog,  method="MDS", distance="bray")
evals <- out.pcoa.log$values[,1]
pdf("iRLS2_PCoA_BrayCurtis_bDiversity_OTU_Group3.pdf")
plot_ordination(pslog, out.pcoa.log, color="Group3", shape="Sesso") +
  labs(col = "Group3", shape = "Sesso")+
  stat_ellipse(aes(group=pslog_sampledata$Group3),type="norm")+
  coord_fixed(sqrt(evals[2] / evals[1]))
dev.off()

##BRAY CURTIS###

#data(enterotype)
#dist_methods <- unlist(distanceMethodList)

bray_dist<-phyloseq::distance(pslog, method="bray")
labels<-attr(bray_dist,"Labels")

table(pslog_sampledata$SampleID==labels)

permanova_bray<-adonis2(bray_dist~Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data=pslog_sampledata)
print(permanova_bray)

groups<-unique(pslog_sampledata$Group3)
pairmanova_res<-list()

for (i in 1:(length(groups)-1)){
  for (j in (i+1):length(groups)) {
  group_i<-groups[i]
  group_j<-groups[j]
  
  subset_indices<-pslog_sampledata$Group3 %in% c(group_i,group_j)
  subset_data<-pslog_sampledata[subset_indices,]
  subset_dist<-as.dist(as.matrix(bray_dist)[subset_indices,subset_indices])
  subset_grouping<-pslog_sampledata$Group3[subset_indices]
  
  result<-adonis2(subset_dist~subset_grouping+Age+Sesso+BMI1+Psichiatrico+Bach, data=subset_data)
  pairmanova_res[[paste(group_i,"vs",group_j)]]<-result
  }
}

pairmanova_res$ALL<-permanova_bray

sink("iRLSC_PCoA_BrayCurtis_bDiversity_OTU_PCoA_PERMANOVA_ALL_Pairwise_Group3.txt")
print(pairmanova_res)
sink()

## PCoA PERMANOVA Analysis

wunifrac_dist<-UniFrac(pslog, weighted=T)
labels<-attr(wunifrac_dist,"Labels")

table(pslog_sampledata$SampleID==labels)

permanova_wunifrac<-adonis2(wunifrac_dist~Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data=pslog_sampledata)
print(permanova_wunifrac)

groups<-unique(pslog_sampledata$Group3)
pairmanova_res<-list()

for (i in 1:(length(groups)-1)){
  for (j in (i+1):length(groups)) {
  group_i<-groups[i]
  group_j<-groups[j]
  
  subset_indices<-pslog_sampledata$Group3 %in% c(group_i,group_j)
  subset_data<-pslog_sampledata[subset_indices,]
  subset_dist<-as.dist(as.matrix(wunifrac_dist)[subset_indices,subset_indices])
  subset_grouping<-pslog_sampledata$Group3[subset_indices]
  
  result<-adonis2(subset_dist~subset_grouping+Age+Sesso+BMI1+Psichiatrico+Bach, data=subset_data)
  pairmanova_res[[paste(group_i,"vs",group_j)]]<-result
  }
}

pairmanova_res$ALL<-permanova_wunifrac

sink("iRLSC_PCoA_bDiversity_OTU_PCoA_wUnifrac_PERMANOVA_ALL_Pairwise_Group3.txt")
print(pairmanova_res)
sink()


# PCoA with unweighted Unifrac ####
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="unifrac")
evals <- out.wuf.log$values$Eigenvalues
pdf("iRLSC_PCoA_bDiversity_OTU_Unifrac_OTU_Group3.pdf")
pcoa_plot<-plot_ordination(pslog, out.wuf.log, color = "Group3",
                shape = "Sesso") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col="Group3", shape="Sesso") + 
  stat_ellipse(aes(group=pslog_sampledata$Group3),type="norm")
  print(pcoa_plot)
dev.off()

## PCoA PERMANOVA Analysis

unifrac_dist<-UniFrac(pslog, weighted=F)
labels<-attr(unifrac_dist,"Labels")

table(pslog_sampledata$SampleID==labels)

permanova_unifrac<-adonis2(unifrac_dist~Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data=pslog_sampledata)
print(permanova_unifrac)

groups<-unique(pslog_sampledata$Group3)
pairmanova_res<-list()

for (i in 1:(length(groups)-1)){
  for (j in (i+1):length(groups)) {
  group_i<-groups[i]
  group_j<-groups[j]
  
  subset_indices<-pslog_sampledata$Group3 %in% c(group_i,group_j)
  subset_data<-pslog_sampledata[subset_indices,]
  subset_dist<-as.dist(as.matrix(unifrac_dist)[subset_indices,subset_indices])
  subset_grouping<-pslog_sampledata$Group3[subset_indices]
  
  result<-adonis2(subset_dist~subset_grouping+Age+Sesso+BMI1+Psichiatrico+Bach, data=subset_data)
  pairmanova_res[[paste(group_i,"vs",group_j)]]<-result
  }
}

pairmanova_res$ALL<-permanova_unifrac

sink("iRLSC_PCoA_bDiversity_OTU_Unifrac_PERMANOVA_ALL_Pairwise_Group3.txt")
print(pairmanova_res)
sink()






