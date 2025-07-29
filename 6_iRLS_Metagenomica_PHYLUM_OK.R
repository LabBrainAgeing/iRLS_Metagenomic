ps_PHYLUM = tax_glom(ps2, "Phylum", NArm = TRUE)

#Group3 ####
library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_PHYLUM
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
write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_PHYLUM_Group3.csv"))
write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_PHYLUM_Group3.csv" )) 
  
counts_normalized <- counts(diagdds_anova, normalized=TRUE)
counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
colnames(counts_normalized)[1] <- "OTU"
taxa <- as.data.frame(physeq_sub@tax_table@.Data)
taxa <- data.frame(rownames(taxa),taxa)
colnames(taxa)[1] <- "OTU"
counts_taxa <- merge(taxa,counts_normalized,by="OTU")
dim(counts_taxa)

#library(RColorBrewer)
#base_palette<-brewer.pal(11,"RdBu")
#counts_taxa$Phylum<-factor(counts_taxa$Phylum)
#palettef<-colorRampPalette(base_palette)(nlevels(counts_taxa$Phylum))

library(pals)
pdf("AbundanceDistribution_PHYLUM_Group3.pdf",width=10,height=6)
par(mfrow=c(1,3))
counts_taxa$Phylum<-factor(counts_taxa$Phylum)
r1<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(counts_taxa[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-counts_taxa$Phylum
barplot(as.matrix(toplot),col=polychrome(nrow(toplot)))

counts_taxa_rel<-apply(counts_taxa[,9:ncol(counts_taxa)],2,FUN=function(x) x/sum(x))
r1<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="CTRL","SampleID"]])
r2<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="insonne","SampleID"]])
r3<-rowMeans(counts_taxa_rel[,sampledata_filt[sampledata_filt$Group3=="iRLS","SampleID"]])
toplot<-data.frame(CTRL=r1,insonne=r2,iRLS=r3)
rownames(toplot)<-counts_taxa$Phylum
barplot(as.matrix(toplot),col=polychrome(nrow(toplot)))

plot(x=NULL, xlim=c(0,100),ylim=c(0,100),frame.plot=F, xaxt='n', yaxt='n')
legend("topleft",legend=factor(rownames(toplot)),pch=15,col=polychrome(nrow(toplot)),cex=1,bty="n")
dev.off()




