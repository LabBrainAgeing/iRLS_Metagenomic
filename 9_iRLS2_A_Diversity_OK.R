##############################-
# ALPHA DIVERSITY
##############################
rm(list=ls())
setwd("")
library("nlme")
library("reshape2")
library("dplyr")
library("phyloseq")
library("ggplot2")
library(car)
library(lm.beta)

load("./240513_MMISEQ_Metagenomica_iRLSFull_iSNSIBO_ps_silvaspecies_phyloseq.RData")

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
sampledata <- data.frame(sample_data(ps2))
sampledata_filt<-sampledata[sampledata$Group3%in%c("CTRL","iRLS","insonne"),]

physeq_sub_filt<-prune_samples(sampledata_filt$SampleID, ps2)

#physeq_sub_filt = tax_glom(physeq_sub_filt, "Genus", NArm = TRUE)

index<-c("Shannon","Simpson","Chao1","Observed")
alphastat<-data.frame(matrix(ncol=length(index),nrow=9))
colnames(alphastat)<-c("Shannon","Simpson","Chao1","Observed")
rownames(alphastat)<-c("iRLSvCTR.pval","iRLSvCTR.ave","iRLSvins.pval","iRLSvins.ave","insvCTR.pval","insvCTR.ave","iRLSvCTR.est","iRLSvins.est","insvCTR.est")
# Compute Shannon diversity for each sample:
for (i in index){
ps_alpha_div <- estimate_richness(physeq_sub_filt, split = TRUE, measure = i)
ps_alpha_div$SampleID <- as.factor(rownames(ps_alpha_div))
ps_samp <- sample_data(physeq_sub_filt) %>%
  unclass() %>%
  data.frame() %>%
  left_join(ps_alpha_div, by = "SampleID") %>%
  melt(measure.vars = i,
       variable.name = "diversity_measure",
       value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ps_samp %>%
  group_by(Group3) %>%	#time or Group
  summarise(mean_div = mean(alpha_diversity)) %>%
  arrange(mean_div)

ps_samp$Group3 <- factor(ps_samp$Group3)
ps_samp$Sesso<-factor(ps_samp$Sesso)

# Compare alpha diversities:
# You can do a linear model, a Mixed effect model or what you prefer

# alpha_div_model <- lme(fixed = alpha_diversity ~ Group, data = ps_samp,random = ~ 1 | Time)

mod1<-ps_samp[ps_samp$Group3=="iRLS"|ps_samp$Group3=="CTRL",]
mod1<-droplevels(mod1)
alpha_div_model1 <- lm(alpha_diversity ~ Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data = mod1)
alphastat[1,i]<-round(Anova(alpha_div_model1, type="II",test.statistic="Wald")[[4]][1],3)
ave1<-round(mean(mod1[mod1$Group3=="iRLS","alpha_diversity"])-mean(mod1[mod1$Group3=="CTRL","alpha_diversity"]),3)
alphastat[2,i]<-paste0(ave1," (",round(mean(mod1[mod1$Group3=="iRLS","alpha_diversity"]),3)," - ",round(mean(mod1[mod1$Group3=="CTRL","alpha_diversity"]),3),")")
alphastat[7,i]<-round(summary(lm.beta(alpha_div_model1))[[4]][2],3)

mod2<-ps_samp[ps_samp$Group3=="iRLS"|ps_samp$Group3=="insonne",]
mod2<-droplevels(mod2)
alpha_div_model2 <- lm(alpha_diversity ~ Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data = mod2)
alphastat[3,i]<-round(Anova(alpha_div_model2, type="II",test.statistic="Wald")[[4]][1],3)
ave2<-round(mean(mod2[mod2$Group3=="iRLS","alpha_diversity"])-mean(mod2[mod2$Group3=="insonne","alpha_diversity"]),3)
alphastat[4,i]<-paste0(ave2," (",round(mean(mod2[mod2$Group3=="iRLS","alpha_diversity"]),3)," - ",round(mean(mod2[mod2$Group3=="insonne","alpha_diversity"]),3),")")
alphastat[8,i]<-round(summary(lm.beta(alpha_div_model2))[[4]][2],3)

mod3<-ps_samp[ps_samp$Group3=="insonne"|ps_samp$Group3=="CTRL",]
mod3<-droplevels(mod3)
alpha_div_model3 <- lm(alpha_diversity ~ Group3+Age+Sesso+BMI1+Psichiatrico+Bach, data = mod3)
alphastat[5,i]<-round(Anova(alpha_div_model3, type="II",test.statistic="Wald")[[4]][1],3)
ave3<-round(mean(mod3[mod3$Group3=="insonne","alpha_diversity"])-mean(mod3[mod3$Group3=="CTRL","alpha_diversity"]),3)
alphastat[6,i]<-paste0(ave3," (",round(mean(mod3[mod3$Group3=="insonne","alpha_diversity"]),3)," - ",round(mean(mod3[mod3$Group3=="CTRL","alpha_diversity"]),3),")")
alphastat[9,i]<-round(summary(lm.beta(alpha_div_model3))[[4]][2],3)

# fitted values, with error bars
# alpha_ps <- ps_samp %>% left_join(new_data)

ps_samp$Group<-factor(ps_samp$Group,levels=c("iRLS","insonne","CTRL"))
pdf(paste0("iRLS2_iSNSIBO_OTU_alpha_Group3_",i,".pdf"))
boxplot(ps_samp$alpha_diversity ~ ps_samp$Group, ylab=i, col=c("#fdc086","#beaed4","#7fc97f"), outline=T,xlab="",names=c("RLS","INS","CTRL"),ylim=c(2,6))  
stripchart(ps_samp$alpha_diversity ~ ps_samp$Group, add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
for (n in 1:length(levels(factor(ps_samp$Group)))) {
 box_x<-n
 temp<-ps_samp$alpha_diversity
 box_y<-mean(temp[ps_samp$Group==levels(factor(ps_samp$Group))[n]], na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=3, col="black")
 points(box_x, box_y, pch=4,cex=1.3, lwd=2, col="dodgerblue4")
 }
dev.off()

# Association with clinical features: ####
iRLS<-ps_samp[ps_samp$Group3=="iRLS",]

pdf(paste0("iRLS2_iSNSIBO_OTU_PSQI_",i,".pdf"))
plot(iRLS$PSQI,iRLS$alpha_diversity,xlab="PSQI_score",ylab=i)
fit<-lm(iRLS$alpha_diversity ~ iRLS$PSQI+iRLS$Age+iRLS$Sesso+iRLS$BMI1+iRLS$Psichiatrico+iRLS$Bach)
fit1<-lm(iRLS$alpha_diversity ~ iRLS$PSQI)
abline(fit1, col="cyan")
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," beta: ",round(summary(fit)[[4]][2],3)), side=3)
dev.off()

pdf(paste0("iRLS2_iSNSIBO_OTU_IRLSS1_",i,".pdf"))
plot(iRLS$IRLSS1,iRLS$alpha_diversity,xlab="IRLSS1_score",ylab=i)
fit<-lm(iRLS$alpha_diversity ~ iRLS$IRLSS1+iRLS$Age+iRLS$Sesso+iRLS$BMI1+iRLS$Psichiatrico+iRLS$Bach)
fit1<-lm(iRLS$alpha_diversity ~ iRLS$IRLSS1)
abline(fit1, col="cyan")
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," beta: ",round(summary(fit)[[4]][2],3)), side=3)
dev.off()

pdf(paste0("iRLS2_iSNSIBO_OTU_ISI1_",i,".pdf"))
plot(iRLS$ISI1,iRLS$alpha_diversity,xlab="ISI1_score",ylab=i)
fit<-lm(iRLS$alpha_diversity ~ iRLS$ISI1+iRLS$Age+iRLS$Sesso+iRLS$BMI1+iRLS$Psichiatrico+iRLS$Bach)
fit1<-lm(iRLS$alpha_diversity ~ iRLS$ISI1)
abline(fit1,col="cyan")
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," beta: ",round(summary(fit)[[4]][2],3)), side=3)
dev.off()

pdf(paste0("iRLS2_iSNSIBO_OTU_Durata1_",i,".pdf"))
plot(iRLS$Durata1,iRLS$alpha_diversity,xlab="Disease duration [y]",ylab=i)
points(iRLS[iRLS$EarlyLateonset2=="early","Durata1"],iRLS[iRLS$EarlyLateonset2=="early","alpha_diversity"], pch=16,col="#f46d43")
points(iRLS[iRLS$EarlyLateonset2=="late","Durata1"],iRLS[iRLS$EarlyLateonset2=="late","alpha_diversity"], pch=16,col="#fee090")
fit<-lm(iRLS$alpha_diversity ~ iRLS$Durata1+iRLS$Age+iRLS$Sesso+iRLS$BMI1+iRLS$Psichiatrico+iRLS$Bach)
lm.beta(fit)
fit1<-lm(iRLS$alpha_diversity ~ iRLS$Durata1)
abline(fit1, col="black")
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," Beta:",round(summary(lm.beta(fit))[[4]][2],3)), side=3)
legend("bottomright", legend=c("early","late"),pch=16, col=c("#f46d43","#fee090"))
dev.off()

pdf(paste0("ANGELICA_iRLS2_iSNSIBO_OTU_Durata1_",i,".pdf"))
plot(iRLS$Durata1,iRLS$alpha_diversity,xlab="Disease duration [y]",ylab=i)
points(iRLS[iRLS$SM=="m","Durata1"],iRLS[iRLS$SM=="m","alpha_diversity"], pch=16,col="#f46d43")
points(iRLS[iRLS$SM=="sm","Durata1"],iRLS[iRLS$SM=="sm","alpha_diversity"], pch=16,col="#fee090")
fit<-lm(iRLS$alpha_diversity ~ iRLS$Durata1+iRLS$Age+iRLS$Sesso+iRLS$BMI1+iRLS$Psichiatrico+iRLS$Bach)
lm.beta(fit)
fit1<-lm(iRLS$alpha_diversity ~ iRLS$Durata1)
abline(fit1, col="black")
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," Beta:",round(summary(lm.beta(fit))[[4]][2],3)), side=3)
legend("bottomright", legend=c("m","sm"),pch=16, col=c("#f46d43","#fee090"))
dev.off()


#pdf(paste0("iRLS2_iSNSIBO_OTU_Psichiatrico_",i,".pdf"))
#fit<-lm(alpha_diversity ~ Psichiatrico+Age+Sesso+BMI1, data=iRLS)
#stat<- Anova(fit, type="II",test.statistic="Wald")
#boxplot(iRLS$alpha_diversity ~ iRLS$Psichiatrico, ylab=i, outline=T,xlab="",,col=c("#fed98e","#d95f0e"))
#stripchart(iRLS$alpha_diversity ~ iRLS$Psichiatrico, add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
#mtext(paste0("iRSL_NoVSiRLS_SI pVal: ",stat[[4]][1]),side=3)
#for (n in 1:length(levels(factor(iRLS$Durata2)))) {
# box_x<-n
# temp<-iRLS$alpha_diversity
# box_y<-mean(temp[iRLS$Psichiatrico==levels(factor(iRLS$Psichiatrico))[n]], na.rm=T)
# box_w<-0.4
# segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
# points(box_x, box_y, pch=4, col="cyan",cex=1.3)
# }
#dev.off()

pdf(paste0("iRLS2_iSNSIBO_OTU_Durata2_",i,".pdf"))
iRLS$Durata2<-factor(iRLS$Durata2,levels=c("short","long"))
fit<-lm(alpha_diversity ~ Durata2+Age+Sesso+BMI1+Psichiatrico+Bach, data=iRLS)
stat<- Anova(fit, type="II",test.statistic="Wald")
#iRLS_t<-iRLS[iRLS$Durata2!="",]
#iRLS_t<-droplevels(iRLS_t)
boxplot(iRLS$alpha_diversity ~ iRLS$Durata2, ylab=i, outline=T,xlab="",,col=c("#fee090","#f46d43"),ylim=c(2,6))
stripchart(iRLS$alpha_diversity ~ iRLS$Durata2, add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," Beta:",round(summary(lm.beta(fit))[[4]][2],3)), side=3)
for (n in 1:length(levels(factor(iRLS$Durata2)))) {
 box_x<-n
 temp<-iRLS$alpha_diversity
 box_y<-mean(temp[iRLS$Durata2==levels(factor(iRLS$Durata2))[n]], na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
dev.off()

pdf(paste0("iRLS2_iSNSIBO_OTU_SM_",i,".pdf"))
fit<-lm(alpha_diversity ~ SM+Age+Sesso+BMI1+Psichiatrico+Bach, data=iRLS)
stat<- Anova(fit, type="II",test.statistic="Wald")
boxplot(iRLS$alpha_diversity ~ factor(iRLS$SM,levels=c("sm","m")), ylab=i, outline=T,xlab="",col=c("#fee090","#f46d43"),names=c("sensory-motor","motor"),ylim=c(2,6))
stripchart(iRLS$alpha_diversity ~ factor(iRLS$SM,levels=c("sm","m")), add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
mtext(paste0("pVal:",round(summary(fit)[[4]][23],3)," Beta:",round(summary(lm.beta(fit))[[4]][2],3)), side=3)
for (n in 1:length(levels(factor(iRLS$SM)))) {
 box_x<-n
 temp<-iRLS$alpha_diversity
 box_y<-mean(temp[iRLS$SM==levels(factor(iRLS$SM,levels=c("sm","m")))[n]], na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
dev.off()


pdf(paste0("iRLS2_iSNSIBO_OTU_Onset_",i,".pdf"))
iRLS$EarlyLateonset2<-factor(iRLS$EarlyLateonset2,levels=c("late","early"))
fit<-lm(alpha_diversity ~ EarlyLateonset2+Age+Sesso+BMI1+Psichiatrico+Bach, data=iRLS)
stat<- Anova(fit, type="II",test.statistic="Wald")
boxplot(iRLS$alpha_diversity ~ iRLS$EarlyLateonset2, ylab=i, outline=T,xlab="",,col=c("#fee090","#f46d43"),ylim=c(2,6))
stripchart(iRLS$alpha_diversity ~ iRLS$EarlyLateonset2, add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
mtext(paste0("iRSL_lateVSiRLS_early pVal: ",round(stat[[4]][1],3)," Beta:",round(summary(lm.beta(fit))[[4]][2],3)),side=3)
for (n in 1:length(levels(factor(iRLS$EarlyLateonset2)))) {
 box_x<-n
 temp<-iRLS$alpha_diversity
 box_y<-mean(temp[iRLS$EarlyLateonset2==levels(factor(iRLS$EarlyLateonset2))[n]], na.rm=T)
 box_w<-0.4
 #segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=3, col="black")
 points(box_x, box_y, pch=4, col="dodgerblue4",cex=1.3,lwd=2)
 }
dev.off()

insonne<-ps_samp[ps_samp$Group3=="insonne",]

pdf(paste0("Insonne_iSNSIBO_OTU_PSQI_",i,".pdf"))
plot(insonne$PSQI,insonne$alpha_diversity,xlab="PSQI_score",ylab=i)
fit<-lm(insonne$alpha_diversity ~ insonne$PSQI+insonne$Age+insonne$Sesso+insonne$BMI1+insonne$Psichiatrico+insonne$Bach)
fit1<-lm(insonne$alpha_diversity ~ insonne$PSQI)
abline(fit1, col="cyan")
mtext(paste0("pVal:",summary(fit)[[4]][23], side=3))
dev.off()

#pdf(paste0("Insonne_iSNSIBO_OTU_Psichiatrico_",i,".pdf"))
#fit<-lm(alpha_diversity ~ Psichiatrico+Age+Sesso+BMI1, data=insonne)
#stat<- Anova(fit, type="II",test.statistic="Wald")
#boxplot(insonne$alpha_diversity ~ insonne$Psichiatrico, ylab=i, outline=T,xlab="",,col=c("#fed98e","#d95f0e"))
#stripchart(insonne$alpha_diversity ~ insonne$Psichiatrico, add=T,vertical=T,pch=16, method="jitter",jitter=0.3,cex=.7)
#mtext(paste0("insonne_NoVSiRLS_SI pVal: ",stat[[4]][1]),side=3)
#for (n in 1:length(levels(factor(insonne$Durata2)))) {
# box_x<-n
# temp<-insonne$alpha_diversity
# box_y<-mean(temp[insonne$Psichiatrico==levels(factor(insonne$Psichiatrico))[n]], na.rm=T)
# box_w<-0.4
# segments(box_x-box_w, box_y, box_x+box_w, box_y, lty=2, lwd=2, col="cyan")
# points(box_x, box_y, pch=4, col="cyan",cex=1.3)
# }
#dev.off()

pdf(paste0("Insonne_iSNSIBO_OTU_ISI1_",i,".pdf"))
plot(insonne$ISI1,insonne$alpha_diversity,xlab="ISI1_score",ylab=i)
fit<-lm(insonne$alpha_diversity ~ insonne$ISI1+insonne$Age+insonne$Sesso+insonne$BMI1+insonne$Psichiatrico+insonne$Bach)
fit1<-lm(insonne$alpha_diversity ~ insonne$ISI1)
abline(fit1,col="cyan")
mtext(paste0("pVal:",summary(fit)[[4]][23], side=3))
dev.off()

}
write.table(alphastat,file="iRLS2_alpha_OTU_stats_Group3.txt",sep="\t",row.names=T)
