rm(list=ls())
library(dada2)

# Set number of threads for parallel computing:
# setThreadOptions(numThreads=10)

experiment_folder <-  ""
output_folder <- ""
output_file <- "240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO"
seqtab1 <- readRDS("./iRLS_iSNSIBO_seqtab.rds")
seqtab2 <- readRDS("./240513_MMISEQ_Metagenomica_iRLS_iSNSIBO_seqtab.rds")


setwd(experiment_folder)

rownames(seqtab2)[which(rownames(seqtab2)%in%rownames(seqtab1))]<-paste0(rownames(seqtab2)[which(rownames(seqtab2)%in%rownames(seqtab1))],"b2")

#BIG DATA####
##    Just run it once  #
st.all <- mergeSequenceTables(seqtab1,seqtab2)
# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
# Assign taxonomy

fastaRef <- "./silva_nr99_v138.1_train_set.fa.gz"
speciesRef <- "./silva_species_assignment_v138.1.fa.gz"


#### IMPORTANTE: per preparare il SampleSheet, usa l'ordine dei campioni nel file *_Reads_discarded.txt
#samdf <- read.table("/standard/users/m.bacalini/METAGENOMICS/IRLS2_Analysis/SampleSheet_IRLS_Full_BigData.txt",sep="\t",header=TRUE,quote="")
samdf <- read.table("./SampleSheet_IRLS_Full_BigData_20240724.txt",sep="\t",header=TRUE,quote="")
ls()

#	INTEGRATE SPECIES	####
#		Just run it once				#
#########################	
tax <- assignTaxonomy(seqtab, fastaRef, multithread=TRUE)
tax.plus <- addSpecies(tax, speciesRef, verbose=TRUE)
#
taxa.print_df <- data.frame(tax.plus, stringsAsFactors = T) # Removing sequence rownames for display only
rownames(taxa.print_df) <- NULL
head(taxa.print_df)
summary(taxa.print_df)

# Write to disk
saveRDS(seqtab, paste0(output_folder, output_file, "_Assign_Species_seqtab.rds")) # CHANGE ME to where you want sequence table saved
saveRDS(tax.plus, paste0(output_folder, output_file, "_Assign_Species_tax.rds")) # CHANGE ME ...


# ANALYSIS START ####

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

suppressMessages(library("knitr"))
suppressMessages(library("BiocStyle"))
suppressMessages(library("ggplot2"))
suppressMessages(library(gridExtra))
#suppressMessages(library(dada2))
suppressMessages(library(phyloseq))
suppressMessages(library(DECIPHER))
suppressMessages(library(phangorn))
suppressMessages(library(RcppParallel))


# %%
#####################################
# Check SampleSheet:
#####################################
str(samdf)
dim(samdf)
colnames(samdf)[1] <- "SampleID"
rownames(samdf) <- samdf$SampleID
table(rownames(seqtabNoC)==samdf$SampleID)
samdf<-samdf[match(rownames(seqtabNoC),samdf$SampleID),]
table(rownames(seqtabNoC)==samdf$SampleID)

# %%
##########################
# Here starts the analysis:
##########################
# Create phyloseq object:
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxTab))#,phy_tree(fitGTR$tree)) #PROPAV A TPGLIERE fitGTR
save(ps, file=paste0(output_folder, output_file,"_ps_silvaspecies_phyloseq.RData"))
load("./240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq.RData")

str(ps)

# Filtering of the samples with a small number of reads
min(rowSums(otu_table(ps))) 
max(rowSums(otu_table(ps))) 
data.frame(rownames(seqtabNoC),rowSums(otu_table(ps)))
pdf(paste0(output_file,"_Barplot_ReadsSums.pdf"))
barplot(rowSums(otu_table(ps)),las=2,names.arg=rownames(ps), cex.names=0.4)
dev.off()

ps <- prune_samples(rowSums(otu_table(ps)) > 20000, ps)
str(ps)

##################################
# TAXONOMIC AND PREVALENCE FILTERING
# Explore feature prevalence in the dataset: number of samples in which a taxon appears at least once.
# PREVALENCE = number of samples in which a taxon appears at least once
#################################

rank_names(ps)
table(tax_table(ps)[, "Kingdom"], exclude = NULL) 


table(tax_table(ps)[, "Phylum"], exclude = NULL)


# Remove NA Phyla and keep only Bacteria and Archaea
ps <- subset_taxa(ps, !is.na(Phylum) & Kingdom %in% c("Archaea","Bacteria"))
str(ps)

table(tax_table(ps)[, "Kingdom"], exclude = NULL) 
# %%

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

# Print mean prevalence and total prevalence of reads classified in each phyla
# In how many samples it is present (on average because it is phylum level)
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$TotalAbundance),sum(df1$TotalAbundance))})



#	Define phyla to filter (for now I don't do this)
#	Supervised Prevalence filtering
# ps1 = subset_taxa(ps, !Phylum %in% c("p__Fusobacteria", "p__TM7"))

## Prevalence filtering:
ps1 = ps
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
nrow(prevdf)
nrow(prevdf1)

# %%
###############################
# Define prevalence threshold: keep taxa if present in >= 5% of the samples
# Unsupervised Prevalence Filtering
###############################
prevalenceThreshold =  10
#prevalenceThreshold =  0.10 * nsamples(ps)
#keep all 
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)

# Abundance filtering
#ps_abundance_filter  = transform_sample_counts(ps2, function(x) x / sum(x) )
#ps_abundance_filter = filter_taxa(ps_abundance_filter, function(x) mean(x) > 1e-5, TRUE)
#ps_abundance_filter = filter_taxa(ps2, function(x) max(x) > 9, TRUE)
#keepTaxa = colnames(otu_table(ps_abundance_filter))
#ps2 <- prune_taxa(keepTaxa, ps2)
#str(ps2)

# How many genera would be present before filtering?
length(get_taxa_unique(ps, taxonomic.rank = "Genus"))

# How many genera would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Genus"))
#
# %%
# How many families would be present before filtering?
length(get_taxa_unique(ps, taxonomic.rank = "Family"))

# How many families would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Family"))
#

# How many species would be present before filtering?
length(get_taxa_unique(ps, taxonomic.rank = "Species"))

# How many species would be present after filtering?
length(get_taxa_unique(ps2, taxonomic.rank = "Species"))
# 

# How many OTU would be present before filtering?
dim(tax_table(ps))

# How many OTU would be present after filtering?
dim(tax_table(ps2))
#
# %%


#################################################################################################################################################################################################
# GENUS
#################################################################################################################################################################################################

####################################
# Aggregate
####################################
# AGGREATE at GENUS level
ps_GENUS = tax_glom(ps2, "Genus", NArm = TRUE)

# Group3 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic Group2s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "Group3"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("CTRL","iRLS","insonne"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)

  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES

# Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~Group3+Sesso+Age+BMI1+Psichiatrico+Bach) 
  
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="Group3_iRLS_vs_CTRL", type="apeglm")
  
  #risultati_anova[rownames(res_anova),paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  Group3_means e i Group3_sd
  lista_Group3s <- as.character(unique(sample_data(physeq_sub_filt)$Group3)) # I do this also for numeric columns
  Group3_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  Group3_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_Group3s_i in lista_Group3s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_Group3_i <- subset_samples(physeq_sub_filt, Group3==lista_Group3s_i)
    Group3_means[,as.character(lista_Group3s_i)] <- colMeans(data.frame(otu_table(physeq_sub_Group3_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Group3_i)))))
    Group3_sd[,as.character(lista_Group3s_i)] <- apply(data.frame(otu_table(physeq_sub_Group3_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Group3_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_Group3_i)))
  }
  rownames(Group3_means) <- Group3_means$taxa
  rownames(Group3_sd) <- Group3_sd$taxa
  Group3_means <- Group3_means[,colnames(Group3_means)!="taxa"]
  Group3_sd <- Group3_sd[,colnames(Group3_sd)!="taxa"]
  Group3_means <- t(Group3_means)
  Group3_sd <- t(Group3_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_Group3s <- unique(sampledata_filt[,label_col])
    for (Group3_i in lista_Group3s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", Group3_i)] <- paste0(round(t(Group3_means[as.character(Group3_i), rownames(risultati_anova_sign)]),3),"(",round(t(Group3_sd[as.character(Group3_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$Group3<-factor(sampledata_filt$Group3,levels=c("CTRL","iRLS","insonne"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#7fc97f","#fdc086","#beaed4"),ylab="log10(Abundance)",xlab="",names=c("CTRL","RLS","IN"))
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    for (n in 1:length(levels(factor(sampledata_filt$Group3)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Group3==levels(factor(sampledata_filt$Group3))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
}
  dev.off()
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_Group3_ANGELICA.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$Group3<-factor(sampledata_filt$Group3,levels=c("iRLS","insonne","CTRL"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#fdc086","#beaed4","#7fc97f"),ylab="log10(Abundance)",xlab="",names=c("RLS","INS","CTRL"),ylim=c(0,7))
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    for (n in 1:length(levels(factor(sampledata_filt$Group3)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Group3==levels(factor(sampledata_filt$Group3))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
}
  dev.off()
  
  #Vertical Barplot
 pdf(paste0(output_folder,output_file,"_GENUS_sign_barplorLogFC_Group3.pdf"))
 par(xpd=T)
 tobarplot<-risultati_anova_sign[,c("Log2_FC_Group3_iRLS_vs_CTRL","Log2_FC_Group3_CTRL_vs_insonne")]
 rownames(tobarplot)<-risultati_anova_sign$Genus
 barplot(t(tobarplot), col=c("steelblue","orange"),horiz=T, beside=T, xlim=c(-4,1),yaxt="n")
 text(-3, seq(2, by=3,length=nrow(tobarplot)), labels=rownames(tobarplot),las=2,pos=2)
 legend("topleft",legend=c("Log2_FC_Group3_iRLS_vs_CTRL","Log2_FC_Group3_CTRL_vs_insonne"),col=c("steelblue","orange"),pch=5)
 dev.off() 

 
 rel<-apply(counts_taxa[,9:ncol(counts_taxa)],2,FUN=function(x) x/sum(x))
counts_taxa_rel<-data.frame(counts_taxa[,1:8],rel)

pdf(paste0(output_folder,output_file,"_GENUS_sign_relative_boxplots_Group3.pdf"))
  sel <- counts_taxa_rel[counts_taxa_rel$Genus %in% sel_for_plots,]
  #sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(as.numeric(as.character(sel[f,9:ncol(sel)]))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#998ec3","#f7f7f7","#f1a340"),ylab="log10(Abundance)",xlab="")
    stripchart(as.numeric(as.character(sel[f,9:ncol(sel)]))~ sampledata_filt$Group3,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    for (n in 1:length(levels(factor(sampledata_filt$Group3)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Group3==levels(factor(sampledata_filt$Group3))[n]],na.rm=T)
 box_w<-0.4
 segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="cyan",cex=1.3)
 }
}
  dev.off()
  
sel_for_plots <- as.character(risultati_anova$Genus)

pdf(paste0(output_folder,output_file,"_GENUS_ALL_boxplots_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova[rownames(risultati_anova)==sel_OTU,grep("Wald",colnames(risultati_anova))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#998ec3","#f7f7f7","#f1a340"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    for (n in 1:length(levels(factor(sampledata_filt$Group)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Group==levels(factor(sampledata_filt$Group))[n]],na.rm=T)
 box_w<-0.4
 segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="cyan",cex=1.3)
 }
}
  dev.off()

###AbundanceDistribution_Barplot ####
library(pals)
pdf("AbundanceDistribution_GENUS_Group3.pdf",width=10,height=6)
par(mfrow=c(1,3))
counts_taxa$Genus<-factor(counts_taxa$Genus)
r1<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-counts_taxa$Genus
barplot(as.matrix(toplot),col=rainbow(factor(rownames(toplot))))

counts_taxa_rel<-apply(counts_taxa[,9:ncol(counts_taxa)],2,FUN=function(x) x/sum(x))
r1<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-counts_taxa$Genus
barplot(as.matrix(toplot),col=rainbow(factor(rownames(toplot))))

plot(x=NULL, xlim=c(0,100),ylim=c(0,100))
legend("topleft",legend=factor(rownames(toplot)),pch=13,col=rainbow(factor(rownames(toplot))),cex=.5,bty="n")
dev.off()

pdf("Top20_AbundanceDistribution_GENUS_Group3.pdf",width=10,height=6)
par(mfrow=c(1,3))
counts_taxa$Genus<-factor(counts_taxa$Genus)

absum<-rowSums(counts_taxa[,9:ncol(counts_taxa)])
count_temp<-data.frame(counts_taxa[,1:8],absum)
count_temp<-count_temp[order(count_temp$absum,decreasing=T),]
top20genera<-count_temp$Genus[1:20]

count_filt<-counts_taxa[counts_taxa$Genus%in%top20genera,]
r1<-rowMeans(count_filt[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(count_filt[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(count_filt[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-count_filt$Genus
barplot(as.matrix(toplot),col=polychrome(nrow(toplot)))

counts_filt_rel<-apply(count_filt[,9:ncol(count_filt)],2,FUN=function(x) x/sum(x))
r1<-rowMeans(counts_filt_rel[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(counts_filt_rel[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(counts_filt_rel[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-count_filt$Genus
barplot(as.matrix(toplot),col=polychrome(nrow(toplot)))

plot(x=NULL, xlim=c(0,100),ylim=c(0,100),frame.plot=F, xaxt='n', yaxt='n')
legend("topleft",legend=factor(rownames(toplot)),pch=15,col=polychrome(nrow(toplot)),cex=1,bty="n")
dev.off()


###Disbiosis ####
library(dplyr)
library(car)
counts_normalized <- counts(diagdds_anova, normalized=TRUE)
counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
colnames(counts_normalized)[1] <- "OTU"
taxa <- as.data.frame(physeq_sub@tax_table@.Data)
taxa <- data.frame(rownames(taxa),taxa)
colnames(taxa)[1] <- "OTU"
counts_taxa <- merge(taxa,counts_normalized,by="OTU")
dim(counts_taxa)
sampledata_d <- data.frame(sample_data(physeq_sub_filt))

Phylum<-counts_taxa[,"Phylum"]
countemp<-data.frame(Phylum,counts_taxa[,c(9:ncol(counts_taxa))])
count_p<- countemp %>% group_by(Phylum)%>% summarize(across(everything(),sum, na.rm=TRUE))
count_p<-data.frame(count_p)
table(sampledata_d$SampleID==colnames(count_p)[2:ncol(count_p)])

count_norm<-count_p
for (c in 2:ncol(count_norm)){
count_norm[,c]<-count_norm[,c]/colSums(count_norm[,2:ncol(count_norm)],na.rm=T)[c-1]
}

ratio<-c()
for (s in 2:ncol(count_norm)){
ratio[s-1]<-count_norm[count_norm$Phylum=="Firmicutes",s]/count_norm[count_norm$Phylum=="Bacteroidota",s]
}

stats<-data.frame(matrix(nrow=2,ncol=3))
colnames(stats)<-c("iRLSvsCTRL","iRLSvsinsonne","insonnevsCTRL")
rownames(stats)<-c("Ave","pVal")

fit1<-Anova(lm(ratio[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL"]~sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL","Group3"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL","Sesso"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL","Age"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL","BMI1"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="CTRL","Psichiatrico"]),type="II")
ave1<-round(mean(ratio[sampledata_d$Group3=="iRLS"])-mean(ratio[sampledata_d$Group3=="CTRL"]),3)
stats[2,1]<-fit1$"Pr(>F)"[1]
stats[1,1]<-paste0(ave1," (",round(mean(ratio[sampledata_d$Group3=="iRLS"]),3)," - ",round(mean(ratio[sampledata_d$Group3=="CTRL"]),3),")")
fit2<-Anova(lm(ratio[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne"]~sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne","Group3"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne","Sesso"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne","Age"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne","BMI1"]+sampledata_d[sampledata_d$Group3=="iRLS"|sampledata_d$Group3=="insonne","Psichiatrico"]),type="II")
ave2<-round(mean(ratio[sampledata_d$Group3=="iRLS"])-mean(ratio[sampledata_d$Group3=="insonne"]),3)
stats[2,2]<-fit2$"Pr(>F)"[1]
stats[1,2]<-paste0(ave2," (",round(mean(ratio[sampledata_d$Group3=="iRLS"]),3)," - ",round(mean(ratio[sampledata_d$Group3=="insonne"]),3),")")
fit3<-Anova(lm(ratio[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL"]~sampledata_d[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL","Group3"]+sampledata_d[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL","Sesso"]+sampledata_d[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL","Age"]+sampledata_d[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL","BMI1"]+sampledata_d[sampledata_d$Group3=="insonne"|sampledata_d$Group3=="CTRL","Psichiatrico"]),type="II")
ave3<-round(mean(ratio[sampledata_d$Group3=="insonne"])-mean(ratio[sampledata_d$Group3=="CTRL"]),3)
stats[2,3]<-fit3$"Pr(>F)"[1]
stats[1,3]<-paste0(ave3," (",round(mean(ratio[sampledata_d$Group3=="insonne"]),3)," - ",round(mean(ratio[sampledata_d$Group3=="CTRL"]),3),")")

write.table(stats,file=paste0(output_folder,output_file,"_Dysbiosis_Group3.txt"),sep="\t")

pdf(paste0(output_folder,output_file,"_Dysbiosis_Group3.pdf"))
boxplot(ratio~sampledata_d$Group3,col=c("#998ec3","#f7f7f7","#f1a340"))
stripchart(ratio~sampledata_d$Group3, add=T, vertical=T, method="jitter",pch=16,cex=.7)
dev.off()
}

rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group_i,group_means,group_sd,is_categorical,label_col,lista_groups,lista_groups_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_group_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa,sampledata_d,sampledata_rep,selb1,selb2,rep_dds_diagdds_rep,physeq_sub_rep,stats,fit1,fit2,fit3,ave1,ave2,ave3,count_norm)

# Group3 iRLS Psichiatrico 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic Psichiatricos are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "Psichiatrico"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
  #1 iRLS filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~Psichiatrico+Sesso+Age+BMI1) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="Psichiatrico_SI_vs_NO", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  Psichiatrico_means e i Psichiatrico_sd
  lista_Psichiatricos <- as.character(unique(sample_data(physeq_sub_filt)$Psichiatrico)) # I do this also for numeric columns
  Psichiatrico_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  Psichiatrico_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_Psichiatricos_i in lista_Psichiatricos){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_Psichiatrico_i <- subset_samples(physeq_sub_filt, Psichiatrico==lista_Psichiatricos_i)
    Psichiatrico_means[,as.character(lista_Psichiatricos_i)] <- colMeans(data.frame(otu_table(physeq_sub_Psichiatrico_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Psichiatrico_i)))))
    Psichiatrico_sd[,as.character(lista_Psichiatricos_i)] <- apply(data.frame(otu_table(physeq_sub_Psichiatrico_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Psichiatrico_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_Psichiatrico_i)))
  }
  rownames(Psichiatrico_means) <- Psichiatrico_means$taxa
  rownames(Psichiatrico_sd) <- Psichiatrico_sd$taxa
  Psichiatrico_means <- Psichiatrico_means[,colnames(Psichiatrico_means)!="taxa"]
  Psichiatrico_sd <- Psichiatrico_sd[,colnames(Psichiatrico_sd)!="taxa"]
  Psichiatrico_means <- t(Psichiatrico_means)
  Psichiatrico_sd <- t(Psichiatrico_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_Psichiatricos <- unique(sampledata_filt[,label_col])
    for (Psichiatrico_i in lista_Psichiatricos){
      risultati_anova_sign[,paste0("Mean (st.err.) ", Psichiatrico_i)] <- paste0(round(t(Psichiatrico_means[as.character(Psichiatrico_i), rownames(risultati_anova_sign)]),3),"(",round(t(Psichiatrico_sd[as.character(Psichiatrico_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_Psichiatrico_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_Psichiatrico_iRLS_Group3.csv" ))
  
if (min(risultati_anova$Wald_pval_Psichiatrico_NO_vs_SI,na.rm=T)<alpha) {
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_Psichiatrico_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Psichiatrico,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4),"\nLRT:",round(as.numeric(unlist(LRT_main))[2],4)," ",round(as.numeric(unlist(LRT_main))[3],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Psichiatrico,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    
for (n in 1:length(levels(factor(sampledata_filt$Psichiatrico)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Psichiatrico==levels(factor(sampledata_filt$Psichiatrico))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }    
}
  dev.off()
} else {}

# Group3 insonne Psichiatrico 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic Psichiatricos are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "Psichiatrico"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
 #1 insonne filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("insonne"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~Psichiatrico+Sesso+Age+BMI1) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="Psichiatrico_SI_vs_NO", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  Psichiatrico_means e i Psichiatrico_sd
  lista_Psichiatricos <- as.character(unique(sample_data(physeq_sub_filt)$Psichiatrico)) # I do this also for numeric columns
  Psichiatrico_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  Psichiatrico_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_Psichiatricos_i in lista_Psichiatricos){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_Psichiatrico_i <- subset_samples(physeq_sub_filt, Psichiatrico==lista_Psichiatricos_i)
    Psichiatrico_means[,as.character(lista_Psichiatricos_i)] <- colMeans(data.frame(otu_table(physeq_sub_Psichiatrico_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Psichiatrico_i)))))
    Psichiatrico_sd[,as.character(lista_Psichiatricos_i)] <- apply(data.frame(otu_table(physeq_sub_Psichiatrico_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Psichiatrico_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_Psichiatrico_i)))
  }
  rownames(Psichiatrico_means) <- Psichiatrico_means$taxa
  rownames(Psichiatrico_sd) <- Psichiatrico_sd$taxa
  Psichiatrico_means <- Psichiatrico_means[,colnames(Psichiatrico_means)!="taxa"]
  Psichiatrico_sd <- Psichiatrico_sd[,colnames(Psichiatrico_sd)!="taxa"]
  Psichiatrico_means <- t(Psichiatrico_means)
  Psichiatrico_sd <- t(Psichiatrico_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_Psichiatricos <- unique(sampledata_filt[,label_col])
    for (Psichiatrico_i in lista_Psichiatricos){
      risultati_anova_sign[,paste0("Mean (st.err.) ", Psichiatrico_i)] <- paste0(round(t(Psichiatrico_means[as.character(Psichiatrico_i), rownames(risultati_anova_sign)]),3),"(",round(t(Psichiatrico_sd[as.character(Psichiatrico_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_Psichiatrico_insonne_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_Psichiatrico_insonne_Group3.csv" ))
  
if (min(risultati_anova$Wald_pval_Psichiatrico_SI_vs_NO,na.rm=T)<alpha) {
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_Psichiatrico_insonne_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Psichiatrico,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4),"\nLRT:",round(as.numeric(unlist(LRT_main))[2],4)," ",round(as.numeric(unlist(LRT_main))[3],4)),col=c("#fed98e","#d95f0e"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Psichiatrico,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    
for (n in 1:length(levels(factor(sampledata_filt$Psichiatrico)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Psichiatrico==levels(factor(sampledata_filt$Psichiatrico))[n]],na.rm=T)
 box_w<-0.4
 segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="cyan",cex=1.3)
 }    
}
  dev.off()
} else {}

}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group_i,group_means,group_sd,is_categorical,label_col,lista_groups,lista_groups_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_group_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)

# Group3 iRLS EarlyLateonset2 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic EarlyLateonset2s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "EarlyLateonset2"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~EarlyLateonset2+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="EarlyLateonset2_late_vs_early", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  EarlyLateonset2_means e i EarlyLateonset2_sd
  lista_EarlyLateonset2s <- as.character(unique(sample_data(physeq_sub_filt)$EarlyLateonset2)) # I do this also for numeric columns
  EarlyLateonset2_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  EarlyLateonset2_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_EarlyLateonset2s_i in lista_EarlyLateonset2s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_EarlyLateonset2_i <- subset_samples(physeq_sub_filt, EarlyLateonset2==lista_EarlyLateonset2s_i)
    EarlyLateonset2_means[,as.character(lista_EarlyLateonset2s_i)] <- colMeans(data.frame(otu_table(physeq_sub_EarlyLateonset2_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_EarlyLateonset2_i)))))
    EarlyLateonset2_sd[,as.character(lista_EarlyLateonset2s_i)] <- apply(data.frame(otu_table(physeq_sub_EarlyLateonset2_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_EarlyLateonset2_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_EarlyLateonset2_i)))
  }
  rownames(EarlyLateonset2_means) <- EarlyLateonset2_means$taxa
  rownames(EarlyLateonset2_sd) <- EarlyLateonset2_sd$taxa
  EarlyLateonset2_means <- EarlyLateonset2_means[,colnames(EarlyLateonset2_means)!="taxa"]
  EarlyLateonset2_sd <- EarlyLateonset2_sd[,colnames(EarlyLateonset2_sd)!="taxa"]
  EarlyLateonset2_means <- t(EarlyLateonset2_means)
  EarlyLateonset2_sd <- t(EarlyLateonset2_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_EarlyLateonset2s <- unique(sampledata_filt[,label_col])
    for (EarlyLateonset2_i in lista_EarlyLateonset2s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", EarlyLateonset2_i)] <- paste0(round(t(EarlyLateonset2_means[as.character(EarlyLateonset2_i), rownames(risultati_anova_sign)]),3),"(",round(t(EarlyLateonset2_sd[as.character(EarlyLateonset2_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_EarlyLateonset2_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_EarlyLateonset2_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_EarlyLateonset2_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$EarlyLateonset2<-factor(sampledata_filt$EarlyLateonset2,levels=c("late","early"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$EarlyLateonset2,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$EarlyLateonset2,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
for (n in 1:length(levels(factor(sampledata_filt$EarlyLateonset2)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$EarlyLateonset2==levels(factor(sampledata_filt$EarlyLateonset2))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4,cex=1.3,lwd=2,col="dodgerblue4")
 }
}
  
  dev.off()
  
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group_i,group_means,group_sd,is_categorical,label_col,lista_groups,lista_groups_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_group_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)


# Group3 iRLS IRLSS1 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic IRLSS1s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "IRLSS1"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$IRLSS1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~IRLSS1+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name=label_col)
    risultati_anova[rownames(res_treat), paste0("Log2_FC_", label_col)] <- res_treat$log2FoldChange
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="IRLSS1", type="ashr")
  
 #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
 #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  IRLSS1_means e i IRLSS1_sd
  lista_IRLSS1s <- as.character(unique(sample_data(physeq_sub_filt)$IRLSS1)) # I do this also for numeric columns
  IRLSS1_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  IRLSS1_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_IRLSS1s_i in lista_IRLSS1s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_IRLSS1_i <- subset_samples(physeq_sub_filt, IRLSS1==lista_IRLSS1s_i)
    IRLSS1_means[,as.character(lista_IRLSS1s_i)] <- colMeans(data.frame(otu_table(physeq_sub_IRLSS1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_IRLSS1_i)))))
    IRLSS1_sd[,as.character(lista_IRLSS1s_i)] <- apply(data.frame(otu_table(physeq_sub_IRLSS1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_IRLSS1_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_IRLSS1_i)))
  }
  rownames(IRLSS1_means) <- IRLSS1_means$taxa
  rownames(IRLSS1_sd) <- IRLSS1_sd$taxa
  IRLSS1_means <- IRLSS1_means[,colnames(IRLSS1_means)!="taxa"]
  IRLSS1_sd <- IRLSS1_sd[,colnames(IRLSS1_sd)!="taxa"]
  IRLSS1_means <- t(IRLSS1_means)
  IRLSS1_sd <- t(IRLSS1_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_IRLSS1s <- unique(sampledata_filt[,label_col])
    for (IRLSS1_i in lista_IRLSS1s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", IRLSS1_i)] <- paste0(round(t(IRLSS1_means[as.character(IRLSS1_i), rownames(risultati_anova_sign)]),3),"(",round(t(IRLSS1_sd[as.character(IRLSS1_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_IRLSS1_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_IRLSS1_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_IRLSS1_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$IRLSS1,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="IRLS1_score")
    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$IRLSS1))
  }
  
  dev.off()
  
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group3_i,Group3_means,Group3_sd,is_categorical,label_col,lista_Group3s,lista_Group3s_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_Group3_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)


# Group3 iRLS ISI1 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic ISI1s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "ISI1"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$ISI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~ISI1+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name=label_col)
    risultati_anova[rownames(res_treat), paste0("Log2_FC_", label_col)] <- res_treat$log2FoldChange
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="ISI1", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  ISI1_means e i ISI1_sd
  lista_ISI1s <- as.character(unique(sample_data(physeq_sub_filt)$ISI1)) # I do this also for numeric columns
  ISI1_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  ISI1_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_ISI1s_i in lista_ISI1s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_ISI1_i <- subset_samples(physeq_sub_filt, ISI1==lista_ISI1s_i)
    ISI1_means[,as.character(lista_ISI1s_i)] <- colMeans(data.frame(otu_table(physeq_sub_ISI1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_ISI1_i)))))
    ISI1_sd[,as.character(lista_ISI1s_i)] <- apply(data.frame(otu_table(physeq_sub_ISI1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_ISI1_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_ISI1_i)))
  }
  rownames(ISI1_means) <- ISI1_means$taxa
  rownames(ISI1_sd) <- ISI1_sd$taxa
  ISI1_means <- ISI1_means[,colnames(ISI1_means)!="taxa"]
  ISI1_sd <- ISI1_sd[,colnames(ISI1_sd)!="taxa"]
  ISI1_means <- t(ISI1_means)
  ISI1_sd <- t(ISI1_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_ISI1s <- unique(sampledata_filt[,label_col])
    for (ISI1_i in lista_ISI1s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", ISI1_i)] <- paste0(round(t(ISI1_means[as.character(ISI1_i), rownames(risultati_anova_sign)]),3),"(",round(t(ISI1_sd[as.character(ISI1_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_ISI1_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_ISI1_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_ISI1_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$ISI1,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="ISI1_score")
    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$ISI1))
  }
  
  dev.off()
  
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group3_i,Group3_means,Group3_sd,is_categorical,label_col,lista_Group3s,lista_Group3s_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_Group3_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)

# Group3 insonne ISI1 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic ISI1s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "ISI1"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("insonne"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$ISI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~ISI1+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name=label_col)
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="ISI1", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  ISI1_means e i ISI1_sd
  lista_ISI1s <- as.character(unique(sample_data(physeq_sub_filt)$ISI1)) # I do this also for numeric columns
  ISI1_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  ISI1_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_ISI1s_i in lista_ISI1s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_ISI1_i <- subset_samples(physeq_sub_filt, ISI1==lista_ISI1s_i)
    ISI1_means[,as.character(lista_ISI1s_i)] <- colMeans(data.frame(otu_table(physeq_sub_ISI1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_ISI1_i)))))
    ISI1_sd[,as.character(lista_ISI1s_i)] <- apply(data.frame(otu_table(physeq_sub_ISI1_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_ISI1_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_ISI1_i)))
  }
  rownames(ISI1_means) <- ISI1_means$taxa
  rownames(ISI1_sd) <- ISI1_sd$taxa
  ISI1_means <- ISI1_means[,colnames(ISI1_means)!="taxa"]
  ISI1_sd <- ISI1_sd[,colnames(ISI1_sd)!="taxa"]
  ISI1_means <- t(ISI1_means)
  ISI1_sd <- t(ISI1_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_ISI1s <- unique(sampledata_filt[,label_col])
    for (ISI1_i in lista_ISI1s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", ISI1_i)] <- paste0(round(t(ISI1_means[as.character(ISI1_i), rownames(risultati_anova_sign)]),3),"(",round(t(ISI1_sd[as.character(ISI1_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_ISI1_insonne_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_ISI1_insonne_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_ISI1_insonne_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$ISI1,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="ISI1_score")
    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$ISI1))
  }
  dev.off()
  
  #sel_for_plots <- as.character(risultati_anova[!is.na(risultati_anova$Wald_pval_ISI1),"Genus"])
  #pdf(paste0(output_folder,output_file,"_GENUS_ALL_boxplots_ISI1_insonne_Group3.pdf"))
  #sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  #sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  #for (f in 1:nrow(sel)) {
  #  sel_OTU <- sel[f,"OTU"]
  #  Wald_main <- risultati_anova[rownames(risultati_anova)==sel_OTU,grep("Wald",colnames(risultati_anova))]
  #  #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
  #  plot(sampledata_filt$ISI1,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="ISI1_score")
  # abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$ISI1))
  #}
  #dev.off()
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group3_i,Group3_means,Group3_sd,is_categorical,label_col,lista_Group3s,lista_Group3s_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_Group3_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)


# Group3 iRLS PSQI 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic PSQIs are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "PSQI"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$PSQI),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~PSQI+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name=label_col)
    risultati_anova[rownames(res_treat), paste0("Log2_FC_", label_col)] <- res_treat$log2FoldChange
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #_anova <- lfcShrink(diagdds_anova, coef="PSQI", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  PSQI_means e i PSQI_sd
  lista_PSQIs <- as.character(unique(sample_data(physeq_sub_filt)$PSQI)) # I do this also for numeric columns
  PSQI_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  PSQI_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_PSQIs_i in lista_PSQIs){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_PSQI_i <- subset_samples(physeq_sub_filt, PSQI==lista_PSQIs_i)
    PSQI_means[,as.character(lista_PSQIs_i)] <- colMeans(data.frame(otu_table(physeq_sub_PSQI_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_PSQI_i)))))
    PSQI_sd[,as.character(lista_PSQIs_i)] <- apply(data.frame(otu_table(physeq_sub_PSQI_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_PSQI_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_PSQI_i)))
  }
  rownames(PSQI_means) <- PSQI_means$taxa
  rownames(PSQI_sd) <- PSQI_sd$taxa
  PSQI_means <- PSQI_means[,colnames(PSQI_means)!="taxa"]
  PSQI_sd <- PSQI_sd[,colnames(PSQI_sd)!="taxa"]
  PSQI_means <- t(PSQI_means)
  PSQI_sd <- t(PSQI_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_PSQIs <- unique(sampledata_filt[,label_col])
    for (PSQI_i in lista_PSQIs){
      risultati_anova_sign[,paste0("Mean (st.err.) ", PSQI_i)] <- paste0(round(t(PSQI_means[as.character(PSQI_i), rownames(risultati_anova_sign)]),3),"(",round(t(PSQI_sd[as.character(PSQI_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_PSQI_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_PSQI_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_PSQI_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$PSQI,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="PSQI_score")
    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$PSQI))
  }
  
  dev.off()
  
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group3_i,Group3_means,Group3_sd,is_categorical,label_col,lista_Group3s,lista_Group3s_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_Group3_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)


 # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "PSQI"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("insonne"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$PSQI),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~PSQI+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name=label_col)
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #_anova <- lfcShrink(diagdds_anova, coef="PSQI", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  PSQI_means e i PSQI_sd
  lista_PSQIs <- as.character(unique(sample_data(physeq_sub_filt)$PSQI)) # I do this also for numeric columns
  PSQI_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  PSQI_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_PSQIs_i in lista_PSQIs){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_PSQI_i <- subset_samples(physeq_sub_filt, PSQI==lista_PSQIs_i)
    PSQI_means[,as.character(lista_PSQIs_i)] <- colMeans(data.frame(otu_table(physeq_sub_PSQI_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_PSQI_i)))))
    PSQI_sd[,as.character(lista_PSQIs_i)] <- apply(data.frame(otu_table(physeq_sub_PSQI_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_PSQI_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_PSQI_i)))
  }
  rownames(PSQI_means) <- PSQI_means$taxa
  rownames(PSQI_sd) <- PSQI_sd$taxa
  PSQI_means <- PSQI_means[,colnames(PSQI_means)!="taxa"]
  PSQI_sd <- PSQI_sd[,colnames(PSQI_sd)!="taxa"]
  PSQI_means <- t(PSQI_means)
  PSQI_sd <- t(PSQI_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_PSQIs <- unique(sampledata_filt[,label_col])
    for (PSQI_i in lista_PSQIs){
      risultati_anova_sign[,paste0("Mean (st.err.) ", PSQI_i)] <- paste0(round(t(PSQI_means[as.character(PSQI_i), rownames(risultati_anova_sign)]),3),"(",round(t(PSQI_sd[as.character(PSQI_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_PSQI_insonne_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_PSQI_insonne_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_PSQI_insonne_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$PSQI,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="PSQI_score")
    #abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$PSQI))
  }
  
  dev.off()

#sel_for_plots <- as.character(risultati_anova$Genus)
#  
#  pdf(paste0(output_folder,output_file,"_GENUS_ALL_boxplots_PSQI_insonne_Group3.pdf"))
#  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
#  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
#  for (f in 1:nrow(sel)) {
#    sel_OTU <- sel[f,"OTU"]
#    Wald_main <- risultati_anova[rownames(risultati_anova)==sel_OTU,grep("Wald",colnames(risultati_anova))]
#    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
#    plot(sampledata_filt$PSQI,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as##.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="PSQI_score")
#    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$PSQI))
#  }
#  
#  dev.off()

}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group3_i,Group3_means,Group3_sd,is_categorical,label_col,lista_Group3s,lista_Group3s_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_Group3_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)

# Group3 iRLS Durata2 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic Durata2s are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "Durata2"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~Durata2+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="Durata2_short_vs_long", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  Durata2_means e i Durata2_sd
  lista_Durata2s <- as.character(unique(sample_data(physeq_sub_filt)$Durata2)) # I do this also for numeric columns
  Durata2_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  Durata2_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_Durata2s_i in lista_Durata2s){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_Durata2_i <- subset_samples(physeq_sub_filt, Durata2==lista_Durata2s_i)
    Durata2_means[,as.character(lista_Durata2s_i)] <- colMeans(data.frame(otu_table(physeq_sub_Durata2_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Durata2_i)))))
    Durata2_sd[,as.character(lista_Durata2s_i)] <- apply(data.frame(otu_table(physeq_sub_Durata2_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_Durata2_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_Durata2_i)))
  }
  rownames(Durata2_means) <- Durata2_means$taxa
  rownames(Durata2_sd) <- Durata2_sd$taxa
  Durata2_means <- Durata2_means[,colnames(Durata2_means)!="taxa"]
  Durata2_sd <- Durata2_sd[,colnames(Durata2_sd)!="taxa"]
  Durata2_means <- t(Durata2_means)
  Durata2_sd <- t(Durata2_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_Durata2s <- unique(sampledata_filt[,label_col])
    for (Durata2_i in lista_Durata2s){
      risultati_anova_sign[,paste0("Mean (st.err.) ", Durata2_i)] <- paste0(round(t(Durata2_means[as.character(Durata2_i), rownames(risultati_anova_sign)]),3),"(",round(t(Durata2_sd[as.character(Durata2_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_Durata2_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_Durata2_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_Durata2_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$Durata2<-factor(sampledata_filt$Durata2,levels=c("short","long"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
   # LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Durata2,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Durata2,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
for (n in 1:length(levels(factor(sampledata_filt$Durata2)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$Durata2==levels(factor(sampledata_filt$Durata2))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
}
  
  
  dev.off()
  
}
rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group_i,group_means,group_sd,is_categorical,label_col,lista_groups,lista_groups_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_group_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)

# Group3 iRLS SM 
{
  # %%
  ##################################
  # DIFFERENTIAL ANALYSIS
  # Hierarchical multiple testing: taxonomic SMs are only tested if higher levels are found to be be associated
  # Normalization: we consider a variance stabilizing transformation available in the DESeq2 package.
  # Results are similar to those with log-transform data.
  ##################################
  # Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "SM"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~SM+Sesso+Age+BMI1+Psichiatrico) 
  
  # Wald test (pairwise comparisons):
  diagdds = DESeq(ps_dds, test="Wald", fitType="parametric",sfType="poscounts")
  
  if (is_categorical){
    if (length(unique(sampledata_filt[,label_col]))>2){
      pairs <- combn(unique(sampledata_filt[,label_col]), 2)
    }else{
      pairs <- as.matrix(unique(sampledata_filt[,label_col]))
    }
    
    for (pair_i in 1:ncol(pairs)){
      col_i <- paste0(label_col,'_', pairs[1,pair_i],'_vs_', pairs[2,pair_i])
      if ((sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[1,pair_i])) > 1) &(sum(sample_data(physeq_sub_filt)[,label_col]==as.character(pairs[2,pair_i])) > 1)){
        res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, contrast=c(label_col,as.character(pairs[1,pair_i]),as.character(pairs[2,pair_i])))
        risultati_anova[rownames(res_treat), paste0("Log2_FC_", col_i)] <- res_treat$log2FoldChange
        risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", col_i)] <- res_treat$lfcSE
        risultati_anova[rownames(res_treat), paste0("FC_", col_i)] <- 2**(res_treat$log2FoldChange)
        risultati_anova[rownames(res_treat), paste0("FC_SE_", col_i)] <- 2**(res_treat$lfcSE)
        risultati_anova[rownames(res_treat), paste0("Wald_pval_", col_i)] <- res_treat$pvalue
        risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", col_i)] <- res_treat$padj
      }else{
        risultati_anova[, paste0("Log2_FC_", col_i)] <- NA
        risultati_anova[, paste0("Log2_FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("FC_", col_i)] <- NA
        risultati_anova[, paste0("FC_SE_", col_i)] <- NA
        risultati_anova[, paste0("Wald_pval_", col_i)] <- NA
        risultati_anova[, paste0("Wald_AdjPval_", col_i)] <- NA
      }
    }
  }else{
    res_treat <- results(diagdds, pAdjustMethod="BH", cooksCutoff=FALSE, name="label_col")
    risultati_anova[rownames(res_treat), paste0("Log2_FC_SE_", label_col)] <- res_treat$lfcSE
    risultati_anova[rownames(res_treat), paste0("FC_", label_col)] <- 2**(res_treat$log2FoldChange)
    risultati_anova[rownames(res_treat), paste0("FC_SE_", label_col)] <- 2**(res_treat$lfcSE)
    risultati_anova[rownames(res_treat), paste0("Wald_pval_", label_col)] <- res_treat$pvalue
    risultati_anova[rownames(res_treat), paste0("Wald_AdjPval_", label_col)] <- res_treat$padj
  }
  
  ## LTR test (1 way ANOVA)
  diagdds_anova = DESeq(ps_dds, test="LRT",reduced=~1, fitType="parametric",sfType="poscounts") # reduced= ~ somma delle covariate (che non interessano)
  
  #res_anova <- results(diagdds_anova, pAdjustMethod="BH",cooksCutoff = FALSE)
  
  resultsNames(diagdds_anova)
  
  #res_anova <- lfcShrink(diagdds_anova, coef="SM_sm_vs_m", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  SM_means e i SM_sd
  lista_SMs <- as.character(unique(sample_data(physeq_sub_filt)$SM)) # I do this also for numeric columns
  SM_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  SM_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_SMs_i in lista_SMs){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_SM_i <- subset_samples(physeq_sub_filt, SM==lista_SMs_i)
    SM_means[,as.character(lista_SMs_i)] <- colMeans(data.frame(otu_table(physeq_sub_SM_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_SM_i)))))
    SM_sd[,as.character(lista_SMs_i)] <- apply(data.frame(otu_table(physeq_sub_SM_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_SM_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_SM_i)))
  }
  rownames(SM_means) <- SM_means$taxa
  rownames(SM_sd) <- SM_sd$taxa
  SM_means <- SM_means[,colnames(SM_means)!="taxa"]
  SM_sd <- SM_sd[,colnames(SM_sd)!="taxa"]
  SM_means <- t(SM_means)
  SM_sd <- t(SM_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_SMs <- unique(sampledata_filt[,label_col])
    for (SM_i in lista_SMs){
      risultati_anova_sign[,paste0("Mean (st.err.) ", SM_i)] <- paste0(round(t(SM_means[as.character(SM_i), rownames(risultati_anova_sign)]),3),"(",round(t(SM_sd[as.character(SM_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_SM_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_SM_iRLS_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_SM_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$SM,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$SM,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    
for (n in 1:length(levels(factor(sampledata_filt$SM)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$SM==levels(factor(sampledata_filt$SM))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
}
dev.off()
  
}

rm(alpha,col_i,counts_normalized,counts_taxa,diagdds,diagdds_anova,f,Group_i,group_means,group_sd,is_categorical,label_col,lista_groups,lista_groups_i,pair_i,pairs,physeq_sub,physeq_sub_filt,physeq_sub_group_i,ps_dds,res_anova,res_treat,risultati_anova,risultati_anova_sign,sampledata,sampledata_filt,sel,sel_for_plots,taxa)
