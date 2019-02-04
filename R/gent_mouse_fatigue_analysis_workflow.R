setwd("/Users/David/Desktop/multiomics/oren_tz/fatigue/")
source("~/Desktop/repos/multiomics/R/aux_functions.R")

# Prepare the workspace, preprocess and reshape the data.

# Read the abundance datasets
rnaseq = read.table("rnaseq.txt",header = T,row.names = 1,sep="\t",stringsAsFactors = F)
rnaseq = rnaseq[,-1]
metab_csf = read.delim("metabolomics_targeted_csf.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
metab_plasma = read.delim("metabolomics_targeted_plasma.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
metab_csf = t(metab_csf)
metab_plasma = t(metab_plasma)

# As of Jan 2019 we exclude "anti" columns
metab_csf = metab_csf[,!grepl("ant",colnames(metab_csf),ignore.case = T)]
metab_plasma = metab_plasma[,!grepl("ant",colnames(metab_plasma),ignore.case = T)]
rnaseq = rnaseq[,!grepl("ant",colnames(rnaseq),ignore.case = T)]
boxplot(log(rnaseq[,1:23])) # Jan 2019 - data seems reasonably normalized

# Split into data matrices, each with a case-control dataset in a tissue
split_dataset<-function(x,ctrl_regex = "Ctrl",tissues=c("")){
  ds = list()
  for(tt in tissues){
    currx = as.matrix(x[,grepl(tt,colnames(x))])
    mode(currx)="numeric"
    curry = grepl(ctrl_regex,colnames(currx),ignore.case = T)
    ds[[tt]] = list(currx,curry)
  }
  return(ds)
}

# Put the filtered abundance datasets in one list so that it will be easier to work with
datasets = list()
datasets[["metabolites;plasma"]] = split_dataset(metab_plasma,"dss")[[1]]
datasets[["metabolites;csf"]] = split_dataset(metab_csf,"dss")[[1]]
datasets[["rnaseq;colon"]] = split_dataset(rnaseq,"dss",tissues="colon")[[1]]
datasets[["rnaseq;CP"]] = split_dataset(rnaseq,"dss",tissues="CP")[[1]]
sapply(datasets,function(x)table(x[[2]]))
sapply(datasets,function(x)dim(x[[1]]))

# Normalize the count data to be log scale - useful for comparisons and linear models
for(j in 1:length(datasets)){
  datasets[[j]][[1]] = log(datasets[[j]][[1]]+1,base=2)
}

# Exclude zero SD analytes
dataset_analyte_sd = lapply(datasets,function(x)apply(x[[1]],1,sd))
sapply(dataset_analyte_sd,function(x)table(x==0))
for(nn in names(datasets)){
  datasets[[nn]][[1]] = datasets[[nn]][[1]][dataset_analyte_sd[[nn]]>0,]
}

# Boxplot of all matrices: make sure that transcriptomics data are still well normalized
# after the filters above.
par(mfrow=c(2,2))
for(nn in names(datasets)){
	boxplot(datasets[[nn]][[1]],main=nn,ylab="Log abundance",xlab="Samples",names=F)
}
# RNA-Seq boxplots only
par(mfrow=c(1,2))
for(nn in names(datasets)[3:4]){
  boxplot(datasets[[nn]][[1]],main=nn,ylab="Log abundance",xlab="Samples",names=F)
}

# Keep a version of the data matrices with gene names instead of transcripts
datasets_genes = list()
for(nn in names(datasets)){
  curr_names = rownames(datasets[[nn]][[1]])
  curr_names = sapply(curr_names,function(x)strsplit(x,split="\\.")[[1]][1])
  if(any(!is.na(gene_names[curr_names]))){curr_names = gene_names[curr_names]}
  to_keep = !is.na(curr_names)
  datasets_genes[[nn]] = datasets[[nn]][[1]][to_keep,]
  rownames(datasets_genes[[nn]]) = curr_names[to_keep]
}
sapply(datasets_genes,dim)


# Load the RNASeq biospecimen metadata

# For analysis of colon
rnaseq_colon_bios = read.table("rnaseq_colon_biospec_metadata.txt",sep="\t",header=T)
# match to the column names in our dataset, look at colnames(datasets[[3]][[1]])
rnaseq_colon_bios = rnaseq_colon_bios[!grepl("anti", rnaseq_colon_bios[,2]),]
xx = datasets[["rnaseq;colon"]][[1]]
rnaseq_colon_bios = rnaseq_colon_bios[c(1:10,12:13,15:16),] # based on manual examination of the sample names
rnaseq_colon_bios[,2] = colnames(xx)

# For analysis of CP
rnaseq_cp_bios = read.table("rnaseq_cp_biospec_metadata.txt",sep="\t",header=T)
rnaseq_cp_bios = rnaseq_cp_bios[order(rnaseq_cp_bios[,1]),]
# match to the column names in our dataset, look at colnames(datasets[[3]][[1]])
rnaseq_cp_bios = rnaseq_cp_bios[!grepl("anti", rnaseq_cp_bios[,1]),]
rnaseq_cp_bios = rnaseq_cp_bios[c(1:3,5:10),] # based on manual examination
xx = datasets[["rnaseq;CP"]][[1]] # manual check: the colnames of xx fit the metadata above
colnames(xx)

# PCA of the CP matrix
xx = t(datasets_genes$`rnaseq;CP`)
yy1 = as.factor(rnaseq_cp_bios$Batch)
yy2 = as.factor(grepl("DSS",rnaseq_cp_bios$Sample.ID))
xx_pca = prcomp(xx,retx=T)$x[,1:2]
dev.off()
plot(xx_pca[,1],xx_pca[,2],pch=as.numeric(yy1),col=yy2,lwd=4,
     ylab="PC2",xlab="PC1",cex.lab=1.2)
legend("bottomleft",as.character(unique(yy1)),pch=unique(as.numeric(yy1)),cex=1.3)
legend("topright",as.character(unique(yy2)),fill=unique(yy2),cex=1.3)

# PCA of the colon matrix
xx = t(datasets_genes$`rnaseq;colon`)
yy1 = as.factor(rnaseq_colon_bios$Batch)
yy2 = as.factor(grepl("DSS",rnaseq_colon_bios$Name))
xx_pca = prcomp(xx,retx=T)$x[,1:2]
dev.off()
plot(xx_pca[,1],xx_pca[,2],pch=as.numeric(yy1),col=yy2,lwd=4,
     ylab="PC2",xlab="PC1",cex.lab=1.2)
legend("bottom",as.character(unique(yy1)),pch=unique(as.numeric(yy1)),cex=1.3)
legend("topright",as.character(unique(yy2)),fill=unique(yy2),cex=1.3)

library(lme4)
# Simple differential abundance analysis without adjusting for 
# technical covariates
# Result: mean difference when compared to controls and the significance
# of the statistical test
simple_diff_analysis<-function(x,y,ctrl=FALSE,func=t.test,...){
  x1 = x[y==ctrl]
  x2 = x[y!=ctrl]
  mean_diff = mean(x2)-mean(x1)
  if(sd(x1)==0 || sd(x2)==0){return(c(0,1))}
  pval = func(x1,x2,...)$p.value
  return(c(mean_diff=mean_diff,pval=pval))
}

# # Deprecated for now. We used fixed-effects analysis to adjust for covariates
# # g1 = a factor for random effects
# # g2 = a covariate to add
# simple_diff_analysis_rand_effects<-function(x,y,g1=NULL,g2=NULL,...){
# 	m0 = NULL
#   if(is.null(g2)){
#   	g1 = as.factor(g1)
#   	d = data.frame(x,y,g1)
#   	m = lmer(x~y + (1|g1),data=d,REML=F)
#   	m0 = lmer(x~ (1|g1),data=d, REML =F) 
#   } 
#   if(is.null(g1)){
#   	d = data.frame(x,y,g2)
#   	m = lm(x~y + g2,data=d)
#   }
#   if(!is.null(g2) & !is.null(g1)){
#   	g1 = as.factor(g1)
#   	d = data.frame(x,y,g1,g2)
#   	m = lmer(x~y + g2 + (1|g1),data=d, REML =F)
#   	m0 = lmer(x~ g2 + (1|g1),data=d, REML =F)
#   }
#   if(!is.null(m0)){
#   	an = anova(m,m0)
#   	pval = an[2,8]
#   	fc = summary(m)$coefficients[2,1]
#   	return(c(fc,pval))
#   }
#   m = summary(m)
#   return(m$coefficients[2,c(1,4)])
# }

# A simple function that adjusts for having one or two covariates.
# Current implementation is fine but not generalizable for more covariates,
# which will be implemented later if necessary.
simple_lm_analysis<-function(x,y,g1=NULL,g2=NULL,...){
  d = data.frame(x,y)
  if(is.null(g2) && !is.null(g2)){
  	d = data.frame(x,y,g1)
  } 
  if(is.null(g1)&&!is.null(g1)){
  	d = data.frame(x,y,g2)
  }
  if(!is.null(g2) & !is.null(g1)){
  	d = data.frame(x,y,g1,g2)
  }
  m = summary(lm(x~.,data=d))
  return(m$coefficients[2,c(1,4)])
}

# Simple pairwise tests
diff_res_ttest = list()
for(nn in names(datasets)){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  res = t(apply(curr_x,1,simple_diff_analysis,y=curr_y))
  diff_res_ttest[[nn]] = res
}
diff_res_wilcox = list()
for(nn in names(datasets)){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  res = t(apply(curr_x,1,simple_diff_analysis,y=curr_y,func=wilcox.test))
  diff_res_wilcox[[nn]] = res
}

# Plot the results above. These are baseline for comparisons against
# more "advanced" analyses (such as those that adjust for batch)
l = diff_res_ttest
par(mfrow=c(2,2))
for(nn in names(l)){
	hist(l[[nn]][,2],main=nn,xlab="P-value",xlim=c(0,1))
}
par(mfrow=c(2,2))
for(nn in names(l)){
	qqplot(y=-log(l[[nn]][,2],10),x=-log(runif(10000),10),main=nn,ylab="Sample quantiles",xlab="Theoretical quantiles")
	abline(0,1,lty=2,lwd=2,col="red")
}

# Here we perform the main analysis: a naive t-test vs.
# regression-based adjustment for technical variables such as
# batch and concentration
diff_res = list()
# Non-rnaseq: use simple t-tests
for(nn in names(datasets)[1:2]){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  res = t(apply(curr_x,1,simple_diff_analysis,y=curr_y))
  diff_res[[nn]] = res
}
# rnaseq: use FIXED effects analysis with the metadata
for(nn in names(datasets)[3:4]){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  curr_m = rnaseq_cp_bios
  if(grepl("colon",nn,ignore.case=T)){curr_m = rnaseq_colon_bios}
  g1 = curr_m$Batch
  g2 = curr_m[,ncol(curr_m)]
  res1 = t(apply(curr_x,1,simple_diff_analysis,y=curr_y))
  res2 = t(apply(curr_x,1, simple_lm_analysis,y=curr_y,g1=g1,g2=g2))
  diff_res[[paste(nn,"_naive",sep="")]] = res1
  diff_res[[paste(nn,"_batch_conc_corrected",sep="")]] = res2
}

# Plot the results side by side: useful to show if the adjustment worked
l = diff_res[3:6]
par(mfrow=c(2,2))
for(nn in names(l)){
  curr_title = "P-value histogram: t-test"
  if(grepl("corrected",nn)){
    curr_title = "P-value histogram: adjusted lm"
  }
  if(grepl("colon",nn)){
    curr_title = paste(curr_title,"colon",sep=",")
  }
  else{
    curr_title = paste(curr_title,"CP",sep=",")
  }
	hist(l[[nn]][,2],main=curr_title,xlab="P-value",xlim=c(0,1))
}
# # In case we want qq plots
# par(mfrow=c(2,2))
# for(nn in names(l)){
# 	qqplot(y=-log(l[[nn]][,2],10),x=-log(runif(10000),10),
# 	       main=nn,ylab="Sample quantiles",xlab="Theoretical quantiles")
# 	abline(0,1,lty=2,lwd=2,col="red")
# }

# Specific plots for the CP data alone
par(mfrow=c(1,2))
pvals1 = diff_res[[5]][,2]
pvals2 = diff_res[[6]][,2]
hist(pvals1,main="P-value histogram: t-test",xlab="P-value",xlim=c(0,1))
# qqplot(y=-log(pvals1,10),
#        x=-log(runif(50000),10),main="QQ-plot: ttest",ylim=c(0,6),
#        ylab="Sample quantiles",xlab="Theoretical quantiles",pch=20,cex=0.5)
# abline(0,1,lty=2,lwd=2,col="red")
hist(pvals2,main="P-value histogram: adjusted lm",xlab="P-value",xlim=c(0,1))
# qqplot(y=-log(pvals2,10),
#        x=-log(runif(50000),10),main="QQ-plot: adjusted lm",ylim=c(0,6),
#        ylab="Sample quantiles",xlab="Theoretical quantiles",pch=20,cex=0.5)
# abline(0,1,lty=2,lwd=2,col="red")

# Compare t-test to adjusted results
l = diff_res[3:6]

# Correlations: CP dataset
y1 = -log(l[[3]][,2])
y2 = -log(l[[4]][,2])
plot(y1,y2);abline(0,1,lty=2,lwd=2)
cor(y1,y2,method = "spearman")
cor(y1,y2)
y1 = l[[3]][,1]
y2 = l[[4]][,1]
plot(y1,y2);abline(0,1,lty=2,lwd=2)
cor(y1,y2,method = "spearman")
cor(y1,y2)
# Correlations: colon dataset
y1 = -log(l[[1]][,2])
y2 = -log(l[[2]][,2])
plot(y1,y2);abline(0,1,lty=2,lwd=2)
cor(y1,y2,method = "spearman")
y1 = l[[1]][,1]
y2 = l[[2]][,1]
plot(y1,y2);abline(0,1,lty=2,lwd=2)
cor(y1,y2,method = "spearman")
cor(y1,y2)

# Adjust and select genes based on FDR
FDR_threshold = 0.1
l = diff_res[c(3:6)]
all_ps = unlist(c(sapply(l,function(x)x[,2])))
thr = max(all_ps[p.adjust(all_ps,method="fdr") < FDR_threshold])
sapply(l,function(x,y)sum(x[,2]<y),y=thr)
selected_results_adjusted = sapply(l,function(x,y)x[x[,2]<y,],y=thr)

# # Simple overlap-based enrichment analysis
# all_selected_results = c(selected_results_adjusted)
# bg = union(rownames(datasets[[3]][[1]]),rownames(datasets[[4]][[1]]))
# bg = sapply(bg,function(x)strsplit(x,split="\\.")[[1]][1])
# bg = gene_names[bg]
# bg = bg[!is.na(bg)]
# enrichment_results = c()
# for(cc in names(all_selected_results)){
#   set1 = rownames(all_selected_results[[cc]])
#   set1 = sapply(set1,function(x)strsplit(x,split="\\.")[[1]][1])
#   if(is.element(cc,set=names(selected_results_naive))){
#     selected_results_naive[[cc]] = list(scores=selected_results_naive[[cc]],gene_names=set1)
#   }
#   if(is.element(cc,set=names(selected_results_adjusted))){
#     selected_results_adjusted[[cc]] = list(scores=selected_results_adjusted[[cc]],gene_names=set1)
#   }
#   set1 = set1[!is.na(set1)]
#   enrichment_res = sapply(mm_pathway,simple_enrichment_analysis,set2=set1,bg=bg)  
#   tb = enrichment_list_to_table(enrichment_res)
#   curr_cl = cc
#   tb = cbind(rep(curr_cl,nrow(tb)),tb)
#   enrichment_results = rbind(enrichment_results,tb)
# }
# # compare enrichment results before and after adjustment
# m = enrichment_results[grepl("naive", enrichment_results[,1]),]
# all_ps = as.numeric(m[,3])
# all_qs = p.adjust(all_ps,method='fdr')
# naive_analysis_results = m[all_qs < 0.1,]
# m = enrichment_results[!grepl("naive", enrichment_results[,1]),]
# all_ps = as.numeric(m[,3])
# all_qs = p.adjust(all_ps,method='fdr')
# adj_analysis_results = m[all_qs < 0.2,]

# GSEA: use p-values or fold changes, compare before and after adjustment
fgsea_enrichment_results = list()
fgsea_ranks = list()
for(cc in names(diff_res)[3:6]){
  curr_scores = -log(diff_res[[cc]][,2])
  names(curr_scores) = sapply(names(curr_scores),function(x)strsplit(x,split="\\.")[[1]][1])
  names(curr_scores) = gene_names[names(curr_scores)]
  fgsea_enrichment_results[[paste(cc,"pval_ranking",sep=";")]] = 
    fgsea(mm_pathway,curr_scores,nperm = 50000,minSize = 5,maxSize = 200)
  fgsea_ranks[[paste(cc,"pval_ranking",sep=";")]] = curr_scores
  
  curr_scores = diff_res[[cc]][,1]
  names(curr_scores) = sapply(names(curr_scores),function(x)strsplit(x,split="\\.")[[1]][1])
  names(curr_scores) = gene_names[names(curr_scores)]
  fgsea_enrichment_results[[paste(cc,"fchange_ranking",sep=";")]] = 
    fgsea(mm_pathway,curr_scores,nperm = 50000,minSize = 5,maxSize = 200)  
  fgsea_ranks[[paste(cc,"fchange_ranking",sep=";")]] = curr_scores
}

# Distribution and statistics of the GSEA results
gsea_res = fgsea_enrichment_results$`rnaseq;CP_naive;fchange_ranking`
plot(gsea_res$NES,gsea_res$pval,pch=20,xlab = "Normalized GSEA score",
     ylab = "Empirical p-value",cex.axis=1.2,cex.lab=1.3,main="CP",
     xlim = c(-4,4))
# Same for colon
gsea_res = fgsea_enrichment_results$`rnaseq;colon_naive;fchange_ranking`
plot(gsea_res$NES,gsea_res$pval,pch=20,xlab = "Normalized GSEA score",
     ylab = "Empirical p-value",cex.axis=1.2,cex.lab=1.3)

# Interpret the GSEA results
# Transform data frames to matricex
fgsea_matrix_enrichment_results = lapply(fgsea_enrichment_results,as.matrix)
adjust_fgsea_res<-function(x,thr=0.01){
  res = x[x[,3] <= thr,]
  return(res)
}
# Select adjusted results
fgsea_matrix_enrichment_results = lapply(fgsea_matrix_enrichment_results,adjust_fgsea_res)
sapply(fgsea_matrix_enrichment_results,dim)
# Save the analysis
save(diff_res,fgsea_ranks,
     selected_results_adjusted,selected_results_naive,
     enrichment_results,fgsea_enrichment_results,fgsea_matrix_enrichment_results,
     file="diff_abundance_analysis_results.RData")

# Compare CP results befor and after the analysis
x1 = unlist(fgsea_matrix_enrichment_results[["rnaseq;CP_batch_conc_corrected;fchange_ranking"]][,1])
x2 = unlist(fgsea_matrix_enrichment_results[["rnaseq;CP_naive;fchange_ranking"]][,1])
length(intersect(x1,x2))
length(x1);length(x2)
x1[grepl("GLUTA",x1)]

# Compare colon results befor and after the analysis
x1 = unlist(fgsea_matrix_enrichment_results[["rnaseq;colon_batch_conc_corrected;fchange_ranking"]][,1])
x2 = unlist(fgsea_matrix_enrichment_results[["rnaseq;colon_naive;fchange_ranking"]][,1])
length(intersect(x1,x2))
length(x1);length(x2)
x1[grepl("GLUTA",x1)]

# Plotting of selected GSEA results
library(data.table);library(fgsea);library(ggplot2);library(gplots)
specific_pathway = "PANTHER_MM_THYROTROPIN-RELEASING_HORMONE_RECEPTOR_SIGNALING_PATHWAY" 
specific_pathway = "WIKIPATHWAYS_MM_MAPK_SIGNALING_PATHWAY-WP382"  
specific_pathway = "PANTHER_MM_DOPAMINE_RECEPTOR_MEDIATED_SIGNALING_PATHWAY" 
specific_pathway = "BIOCARTA_MM_REGULATION_OF_CK1_CDK5_BY_TYPE_1_GLUTAMATE_RECEPTORS" 
specific_pathway = "PANTHER_MM_IONOTROPIC_GLUTAMATE_RECEPTOR_PATHWAY" 

# Further display items: CP
ranks_name = names(diff_res)[5]
specific_pathway = "PANTHER_MM_IONOTROPIC_GLUTAMATE_RECEPTOR_PATHWAY" 
curr_ranks_name = names(fgsea_ranks)[grepl(ranks_name,names(fgsea_ranks)) &
                                  grepl("fchange",names(fgsea_ranks))]
curr_ranks = fgsea_ranks[[curr_ranks_name]]
plotEnrichment(mm_pathway[[specific_pathway]],
               curr_ranks) + labs(title=specific_pathway)
all_gsea_res = as.matrix(fgsea_enrichment_results[[curr_ranks_name]])
rownames(all_gsea_res) = all_gsea_res[,1]
curr_genes = all_gsea_res[specific_pathway,8][[1]]
curr_ranks[curr_genes]
m = datasets_genes[[4]][curr_genes,]
m_y = datasets[[4]][[2]]
curr_meta = rnaseq_cp_bios
colnames(m) = rep("control",ncol(m))
colnames(m)[m_y] = "case"
heatmap.2(m,trace = "none",scale = "row",col="bluered")
# Print the selected GSEA results table in a nice table format
write.table(fgsea_matrix_enrichment_results[["rnaseq;CP_batch_conc_corrected;fchange_ranking"]],
            sep="\t", file = "CP_GSEA_results.txt",row.names = F,col.names = T)

# Further display items: colon
ranks_name = names(diff_res)[3]
specific_pathway = "INOH_MM_GLUTAMATE_GLUTAMINE_METABOLISM"
curr_ranks_name = names(fgsea_ranks)[grepl(ranks_name,names(fgsea_ranks)) &
                                       grepl("fchange",names(fgsea_ranks))]
curr_ranks = fgsea_ranks[[curr_ranks_name]]
plotEnrichment(mm_pathway[[specific_pathway]],
               curr_ranks) + labs(title=specific_pathway)
all_gsea_res = as.matrix(fgsea_enrichment_results[[curr_ranks_name]])
rownames(all_gsea_res) = all_gsea_res[,1]
curr_genes = all_gsea_res[specific_pathway,8][[1]]
curr_ranks[curr_genes]
m = datasets_genes[[3]][curr_genes,]
m_y = datasets[[3]][[2]]
curr_meta = rnaseq_cp_bios
colnames(m) = rep("control",ncol(m))
colnames(m)[m_y] = "case"
heatmap.2(m,trace = "none",scale = "row",col="bluered",mar=c(8,8))
# Print the selected GSEA results table in a nice table format
write.table(fgsea_matrix_enrichment_results[["rnaseq;colon_naive;fchange_ranking"]],
            sep="\t", file = "colon_GSEA_results.txt",row.names = F,col.names = T,quote=F)


