#Group3 Treatments ####
#Gabaergic Treatment Effect ####

library("reshape2")
  library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "gabaergici"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  sampledata$gabaergici<-ifelse(sampledata$gabaergici==1,"Treat","Untreat")
  
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3=="iRLS",]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  sample_data(physeq_sub_filt)$gabaergici<-ifelse(sample_data(physeq_sub_filt)$gabaergici==1,"Treat","Untreat")

  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES

# Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~gabaergici+Sesso+Age+BMI1+Psichiatrico+Bach) 
  
  
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
  
  #res_anova <- lfcShrink(diagdds_anova, coef="gabaergici_iRLS_vs_CTRL", type="apeglm")
  
  #risultati_anova[rownames(res_anova),paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  gabaergici_means e i gabaergici_sd
  lista_gabaergicis <- as.character(unique(sample_data(physeq_sub_filt)$gabaergici)) # I do this also for numeric columns
  gabaergici_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  gabaergici_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_gabaergicis_i in lista_gabaergicis){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_gabaergici_i <- subset_samples(physeq_sub_filt, gabaergici==lista_gabaergicis_i)
    gabaergici_means[,as.character(lista_gabaergicis_i)] <- colMeans(data.frame(otu_table(physeq_sub_gabaergici_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_gabaergici_i)))))
    gabaergici_sd[,as.character(lista_gabaergicis_i)] <- apply(data.frame(otu_table(physeq_sub_gabaergici_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_gabaergici_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_gabaergici_i)))
  }
  rownames(gabaergici_means) <- gabaergici_means$taxa
  rownames(gabaergici_sd) <- gabaergici_sd$taxa
  gabaergici_means <- gabaergici_means[,colnames(gabaergici_means)!="taxa"]
  gabaergici_sd <- gabaergici_sd[,colnames(gabaergici_sd)!="taxa"]
  gabaergici_means <- t(gabaergici_means)
  gabaergici_sd <- t(gabaergici_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_gabaergicis <- unique(sampledata_filt[,label_col])
    for (gabaergici_i in lista_gabaergicis){
      risultati_anova_sign[,paste0("Mean (st.err.) ", gabaergici_i)] <- paste0(round(t(gabaergici_means[as.character(gabaergici_i), rownames(risultati_anova_sign)]),3),"(",round(t(gabaergici_sd[as.character(gabaergici_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_gabaergici_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_gabaergici_iRLS_Group3.csv" ))

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
  #sel_for_plots<-"Lachnoclostridium"
 
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_gabaergici_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$gabaergici<-factor(sampledata_filt$gabaergici,levels=c("Untreat","Treat"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$gabaergici,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4),"\nLRT:",round(as.numeric(unlist(LRT_main))[2],4)," ",round(as.numeric(unlist(LRT_main))[3],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$gabaergici,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    
for (n in 1:length(levels(factor(sampledata_filt$gabaergici)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$gabaergici==levels(factor(sampledata_filt$gabaergici))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }    
}
  dev.off()
  

#Dopa-Agonist Effect####

library("reshape2")
library("DESeq2")
  
  alpha <- 0.05
  physeq_sub <- ps_GENUS
  label_col <- "dopaagonisti"
  is_categorical <- T
  sampledata <- data.frame(sample_data(ps))
  sampledata$dopaagonisti<-ifelse(sampledata$dopaagonisti==1,"Treat","Untreat")
  
  
  #filter for samples of interest
  sampledata_filt<-sampledata[sampledata$Group3=="iRLS",]
  sampledata_filt<-sampledata_filt[!is.na(sampledata_filt$BMI1),]
  physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, physeq_sub)
  sample_data(physeq_sub_filt)$dopaagonisti<-ifelse(sample_data(physeq_sub_filt)$dopaagonisti==1,"Treat","Untreat")

  
  risultati_anova <- data.frame(tax_table(physeq_sub)) # NON FUNZIONA QUANDO AGGREGHIAMO PER LE FAMIGLIE (Continua a non andare)
  risultati_anova  <- as.data.frame(physeq_sub_filt@tax_table@.Data) #L'ORDINE DEI TAXA è LO STESSO CHE SI OTTIENE CON TAXA_NAMES

# Normalize counts:
  ps_dds <- phyloseq_to_deseq2(physeq_sub_filt, design = ~dopaagonisti+Sesso+Age+BMI1+Psichiatrico+Bach) 
  
  
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
  
  #res_anova <- lfcShrink(diagdds_anova, coef="dopaagonisti_iRLS_vs_CTRL", type="apeglm")
  
  #risultati_anova[rownames(res_anova),paste0("LRT_Log2_FC")] <- res_anova$log2FoldChange
  #risultati_anova[rownames(res_anova),paste0("LRT_pval")] <- res_anova$pvalue
  #risultati_anova[rownames(res_anova),paste0("LRT_AdjPval")] <- res_anova$padj
  #colnames(risultati_anova) <- gsub("label_col",label_col,colnames(risultati_anova))
  # Aggiungo Mean (st.err)
  risultati_anova_sign <- risultati_anova[which(apply(risultati_anova[,grep("Wald_",colnames(risultati_anova))],MARGIN=1,FUN=function(x) min(x,na.rm=TRUE)) < alpha),] # alpha va fissato prima (per es. alpha=0.05)
  
  #####  dopaagonisti_means e i dopaagonisti_sd
  lista_dopaagonistis <- as.character(unique(sample_data(physeq_sub_filt)$dopaagonisti)) # I do this also for numeric columns
  dopaagonisti_means <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  dopaagonisti_sd <- data.frame(taxa=colnames(otu_table(physeq_sub_filt)))
  
  for (lista_dopaagonistis_i in lista_dopaagonistis){ # aggiungere % media e SD per ogni gruppo
    physeq_sub_dopaagonisti_i <- subset_samples(physeq_sub_filt, dopaagonisti==lista_dopaagonistis_i)
    dopaagonisti_means[,as.character(lista_dopaagonistis_i)] <- colMeans(data.frame(otu_table(physeq_sub_dopaagonisti_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_dopaagonisti_i)))))
    dopaagonisti_sd[,as.character(lista_dopaagonistis_i)] <- apply(data.frame(otu_table(physeq_sub_dopaagonisti_i))*100/as.vector(rowSums(data.frame(otu_table(physeq_sub_dopaagonisti_i)))),2,FUN=sd)/sqrt(nrow(sample_data(physeq_sub_dopaagonisti_i)))
  }
  rownames(dopaagonisti_means) <- dopaagonisti_means$taxa
  rownames(dopaagonisti_sd) <- dopaagonisti_sd$taxa
  dopaagonisti_means <- dopaagonisti_means[,colnames(dopaagonisti_means)!="taxa"]
  dopaagonisti_sd <- dopaagonisti_sd[,colnames(dopaagonisti_sd)!="taxa"]
  dopaagonisti_means <- t(dopaagonisti_means)
  dopaagonisti_sd <- t(dopaagonisti_sd)
  
  if (nrow(risultati_anova_sign)>=1){
    lista_dopaagonistis <- unique(sampledata_filt[,label_col])
    for (dopaagonisti_i in lista_dopaagonistis){
      risultati_anova_sign[,paste0("Mean (st.err.) ", dopaagonisti_i)] <- paste0(round(t(dopaagonisti_means[as.character(dopaagonisti_i), rownames(risultati_anova_sign)]),3),"(",round(t(dopaagonisti_sd[as.character(dopaagonisti_i),rownames(risultati_anova_sign)]),3),")")
    }
  }
  
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_GENUS_dopaagonisti_iRLS_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_GENUS_dopaagonisti_iRLS_Group3.csv" ))

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
  #sel_for_plots<-"Lachnoclostridium"
 
  pdf(paste0(output_folder,output_file,"_GENUS_sign_boxplots_dopaagonisti_iRLS_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Genus %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$dopaagonisti<-factor(sampledata_filt$dopaagonisti,levels=c("Untreat","Treat"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$dopaagonisti,outline=T,main=paste0(sel[f,]$Genus,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4),"\nLRT:",round(as.numeric(unlist(LRT_main))[2],4)," ",round(as.numeric(unlist(LRT_main))[3],4)),col=c("#fee090","#f46d43"),ylab="log10(Abundance)",xlab="")
    stripchart(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$dopaagonisti,col= "black", vertical = TRUE, method = "jitter", add = TRUE, pch = 16,cex=0.7)
    
for (n in 1:length(levels(factor(sampledata_filt$dopaagonisti)))) {
 box_x<-n
 temp<-log10(as.numeric(as.character(sel[f,9:ncol(sel)])))
 box_y<-mean(temp[sampledata_filt$dopaagonisti==levels(factor(sampledata_filt$dopaagonisti))[n]],na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue",cex=1.3,lwd=2)
 }    
}
  dev.off()