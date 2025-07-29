ps_FAMILY = tax_glom(ps2, "Family", NArm = TRUE)

# Group3 
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
  physeq_sub <- ps_FAMILY
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
  write.csv(risultati_anova, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Output_AssignSpecies_FAMILY_Group3.csv"))
  write.csv(risultati_anova_sign, file=paste0(output_folder,output_file,"_DifferentialAnalisis_Sign_Output_AssignSpecies_FAMILY_Group3.csv" ))
  
  counts_normalized <- counts(diagdds_anova, normalized=TRUE)
  counts_normalized <- data.frame(rownames(counts_normalized), counts_normalized)
  colnames(counts_normalized)[1] <- "OTU"
  taxa <- as.data.frame(physeq_sub@tax_table@.Data)
  taxa <- data.frame(rownames(taxa),taxa)
  colnames(taxa)[1] <- "OTU"
  counts_taxa <- merge(taxa,counts_normalized,by="OTU")
  dim(counts_taxa)
  sampledata <- data.frame(sample_data(ps2))
  sel_for_plots <- as.character(risultati_anova_sign$Family)
  
  pdf(paste0(output_folder,output_file,"_FAMILY_sign_boxplots_Group3.pdf"))
  sel <- counts_taxa[counts_taxa$Family %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$Group3<-factor(sampledata_filt$Group3,levels=c("CTRL","iRLS","insonne"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Family,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#7fc97f","#fdc086","#beaed4"),ylab="log10(Abundance)",xlab="", names=c("CTRL","RLS","IN"))
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
  
   pdf(paste0(output_folder,output_file,"_FAMILY_sign_boxplots_Group3_ANGELICA.pdf"))
  sel <- counts_taxa[counts_taxa$Family %in% sel_for_plots,]
  sel[] <- lapply(sel, function(x) (gsub("\\<0\\>", NA, x)))
  for (f in 1:nrow(sel)) {
    sel_OTU <- sel[f,"OTU"]
    sampledata_filt$Group3<-factor(sampledata_filt$Group3,levels=c("iRLS","insonne","CTRL"))
    Wald_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("Wald",colnames(risultati_anova_sign))]
    #LRT_main <- risultati_anova_sign[rownames(risultati_anova_sign)==sel_OTU,grep("LRT",colnames(risultati_anova_sign))]
    boxplot(log10(as.numeric(as.character(sel[f,9:ncol(sel)])))~ sampledata_filt$Group3,outline=T,main=paste0(sel[f,]$Family,"\nWald:",round(as.numeric(unlist(Wald_main))[1],4)," ",round(as.numeric(unlist(Wald_main))[2],4)),col=c("#fdc086","#beaed4","#7fc97f"),ylab="log10(Abundance)",xlab="", names=c("RLS","INS","CTRL"),ylim=c(0,7))
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


###Grafico x Cami ####
library(pals)
pdf("AbundanceDistribution_FAMILY_Group3.pdf",width=10,height=6)
par(mfrow=c(1,3))
counts_taxa$Genus<-factor(counts_taxa$Family)
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

pdf("Top20_AbundanceDistribution_FAMILY_Group3.pdf",width=10,height=6)
par(mfrow=c(1,3))
counts_taxa$Genus<-factor(counts_taxa$Family)

absum<-rowSums(counts_taxa[,9:ncol(counts_taxa)])
count_temp<-data.frame(counts_taxa[,1:8],absum)
count_temp<-count_temp[order(count_temp$absum,decreasing=T),]
top20genera<-count_temp$Family[1:20]

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
