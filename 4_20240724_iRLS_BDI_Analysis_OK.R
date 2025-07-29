# BDI.II

# Differential analysis:
  library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "BDI.II"
  is_categorical <- F
  sampledata <- data.frame(sample_data(ps))
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3 %in% c("iRLS"),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BDI.II),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES
  
  # Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~BDI.II+Sesso+Age+BMI1) 
  
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
  
  #_anova <- lfcShrink(diagdds_anova, coef="BDI.II", type="ashr")
  
  #risultati_anova[rownames(res_anova), paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #######
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_BDI.II_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_BDI.II__Group3.csv" ))
  
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
  
  sampledata_insonne<-sampledata[sampledata$Group3 %in% c("insonne"),]
  sampledata_insonne<-sampledata_insonne[!is.na(sampledata_insonne$BDI.II),]
  physeq_insonne<-prune_samples(sampledata_insonne$SampleID, physeq_sub)
  ps_dds_insonne <- phyloseq_to_deseq2(physeq_insonne, design = ~BDI.II+Sesso+Age+BMI1) 
  diagdds_insonne = DESeq(ps_dds_insonne, test="Wald", fitType="parametric",sfType="poscounts")
  counts_normalized_insonne <- counts(diagdds_insonne, normalized=TRUE)
  counts_normalized_insonne <- data.frame(rownames(counts_normalized_insonne), counts_normalized_insonne)
  colnames(counts_normalized_insonne)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa_insonne <- merge(taxa,counts_normalized_insonne,by="OTU")
  
  table(counts_taxa$Genus==counts_taxa_insonne$Genus)
  
  pdf(paste0(output_folder,output_file,"_GENUS_sign_iRLS_insonne_BDI.II_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel_insonne<-counts_taxa_insonne[counts_taxa_insonne$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  sel_insonne[] <- lapply(sel_insonne, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    plot(sampledata_filt$BDI.II,log10(as.numeric(as.character(sel[f,9:ncol(sel)]))),outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),ylab="log10(Abundance)",xlab="BDI.II_score")
    abline(lm(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~sampledata_filt$BDI.II))
    points(sampledata_insonne$BDI.II, log10(as.numeric(as.character(sel_insonne[f,9:ncol(sel_insonne)]))),col="blue")
    abline(lm(log10(as.numeric(as.character(sel_insonne[f,9:ncol(sel_insonne)])))~sampledata_insonne$BDI.II),col="blue")
  }
  
  dev.off()
