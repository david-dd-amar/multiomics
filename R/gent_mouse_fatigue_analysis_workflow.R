setwd("/Users/David/Desktop/multiomics/oren_tz/fatigue/")
source("~/Desktop/repos/multiomics/R/aux_functions.R")
list.files()

rnaseq = read.table("rnaseq.txt",header = T,row.names = 1,sep="\t",stringsAsFactors = F)
rnaseq = rnaseq[,-1]
metab_csf = read.delim("metabolomics_targeted_csf.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
metab_plasma = read.delim("metabolomics_targeted_plasma.txt",header = T,sep="\t",stringsAsFactors = F,row.names = 1)
metab_csf = t(metab_csf)
metab_plasma = t(metab_plasma)

# Exclude "anti" columns
metab_csf = metab_csf[,!grepl("ant",colnames(metab_csf),ignore.case = T)]
metab_plasma = metab_plasma[,!grepl("ant",colnames(metab_plasma),ignore.case = T)]
rnaseq = rnaseq[,!grepl("ant",colnames(rnaseq),ignore.case = T)]
colnames(rnaseq)
colnames(metab_plasma)
colnames(metab_csf)
boxplot(rnaseq[,1:23])

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
datasets = list()
datasets[["metabolites;plasma"]] = split_dataset(metab_plasma,"dss")[[1]]
datasets[["metabolites;csf"]] = split_dataset(metab_csf,"dss")[[1]]
datasets[["rnaseq;colon"]] = split_dataset(rnaseq,"dss",tissues="colon")[[1]]
datasets[["rnaseq;CP"]] = split_dataset(rnaseq,"dss",tissues="CP")[[1]]
sapply(datasets,function(x)table(x[[2]]))
sapply(datasets,function(x)dim(x[[1]]))
for(j in 1:length(datasets)){
  datasets[[j]][[1]] = log(datasets[[j]][[1]]+1,base=2)
}

# Exclude zero SD analytes
dataset_analyte_sd = lapply(datasets,function(x)apply(x[[1]],1,sd))
sapply(dataset_analyte_sd,function(x)table(x==0))
for(nn in names(datasets)){
  datasets[[nn]][[1]] = datasets[[nn]][[1]][dataset_analyte_sd[[nn]]>0,]
}

par(mfrow=c(2,2))
sapply(datasets,function(x)boxplot(x[[1]],names = NA))
sapply(datasets,function(x)table(x[[1]]==0))
sapply(datasets,function(x)nrow(x[[1]]==0))
sapply(datasets,function(x)x[[2]])

# Load the RNASeq biospecimen metadata
# colon
rnaseq_colon_bios = read.table("rnaseq_colon_biospec_metadata.txt",sep="\t",header=T)
# match to the column names in our dataset, look at colnames(datasets[[3]][[1]])
rnaseq_colon_bios = rnaseq_colon_bios[!grepl("anti", rnaseq_colon_bios[,2]),]
# boxplot(Concentration~Batch,data=rnaseq_colon_bios)
# summary(lm(Concentration~Batch,data=rnaseq_colon_bios))
xx = datasets[["rnaseq;colon"]][[1]]
rnaseq_colon_bios = rnaseq_colon_bios[c(1:10,12:13,15:16),] # based on manual examination of the sample names
rnaseq_colon_bios[,2] = colnames(xx)
# CP
rnaseq_cp_bios = read.table("rnaseq_cp_biospec_metadata.txt",sep="\t",header=T)
rnaseq_cp_bios = rnaseq_cp_bios[order(rnaseq_cp_bios[,1]),]
# match to the column names in our dataset, look at colnames(datasets[[3]][[1]])
rnaseq_cp_bios = rnaseq_cp_bios[!grepl("anti", rnaseq_cp_bios[,1]),]
rnaseq_cp_bios = rnaseq_cp_bios[c(1:3,5:10),]
xx = datasets[["rnaseq;CP"]][[1]] # manual check: the colnames of xx fit the metadata above
colnames(xx)

library(lme4)
simple_diff_analysis<-function(x,y,ctrl=FALSE,func=t.test,...){
  x1 = x[y==ctrl]
  x2 = x[y!=ctrl]
  mean_diff = mean(x2)-mean(x1)
  if(sd(x1)==0 || sd(x2)==0){return(c(0,1))}
  pval = func(x1,x2,...)$p.value
  return(c(mean_diff=mean_diff,pval=pval))
}

# g1 = a factor for random effects
# g2 = a covariate to add
simple_diff_analysis_rand_effects<-function(x,y,g1=NULL,g2=NULL,...){
	m0 = NULL
  if(is.null(g2)){
  	g1 = as.factor(g1)
  	d = data.frame(x,y,g1)
  	m = lmer(x~y + (1|g1),data=d,REML=F)
  	m0 = lmer(x~ (1|g1),data=d, REML =F) 
  } 
  if(is.null(g1)){
  	d = data.frame(x,y,g2)
  	m = lm(x~y + g2,data=d)
  }
  if(!is.null(g2) & !is.null(g1)){
  	g1 = as.factor(g1)
  	d = data.frame(x,y,g1,g2)
  	m = lmer(x~y + g2 + (1|g1),data=d, REML =F)
  	m0 = lmer(x~ g2 + (1|g1),data=d, REML =F)
  }
  if(!is.null(m0)){
  	an = anova(m,m0)
  	pval = an[2,8]
  	fc = summary(m)$coefficients[2,1]
  	return(c(fc,pval))
  }
  m = summary(m)
  return(m$coefficients[2,c(1,4)])
}

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


diff_res = list()
# Non-rnaseq: use simple t-tests
for(nn in names(datasets)[1:2]){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  res = t(apply(curr_x,1,simple_diff_analysis,y=curr_y))
  diff_res[[nn]] = res
}
# rnaseq: use rand effects analysis with the metadata
for(nn in names(datasets)[3:4]){
  curr_x = datasets[[nn]][[1]]
  curr_y = datasets[[nn]][[2]]
  curr_m = rnaseq_cp_bios
  if(grepl("colon",nn,ignore.case=T)){curr_m = rnaseq_colon_bios}
  g1 = curr_m$Batch
  g2 = curr_m[,ncol(curr_m)]
  #res1 = t(apply(curr_x,1, simple_diff_analysis_rand_effects,y=curr_y,g1=g1))
  #res2 = t(apply(curr_x,1, simple_diff_analysis_rand_effects,y=curr_y,g2=g2))
  #res3 = t(apply(curr_x,1, simple_diff_analysis_rand_effects,y=curr_y,g1=g1,g2=g2))
  res4 = t(apply(curr_x,1,simple_diff_analysis,y=curr_y))
  res5 = t(apply(curr_x,1, simple_lm_analysis,y=curr_y,g1=g1,g2=g2))
  #diff_res[[paste(nn,"re",sep="")]] = res1
  #diff_res[[paste(nn,"fe",sep="")]] = res2
  #diff_res[[paste(nn,"re,fe",sep="")]] = res3
  diff_res[[paste(nn,"naive",sep="")]] = res4
  diff_res[[paste(nn,"batch_conc_corrected",sep="")]] = res5
}

# For QA and testing
gene = "ENSMUSG00000000159.15"
res5[gene,]
res4[gene,]
x = datasets[[3]][[1]][gene,]
y = datasets[[3]][[2]]
g1 = rnaseq_colon_bios$Batch
g2 = rnaseq_colon_bios[,4]
d = data.frame(x,y,g1,g2)
summary(lm(x~.,data=d))
par(mfrow=c(1,3))
hist(res4[,2]);hist(res5[,2])
qqplot(res4[,2],res5[,2]);abline(0,1,lwd=2,lty=2)
cor(res4[,2],res5[,2],method="spearman")
cor(res4[,1],res5[,1])
table(res4[,2]<1e-4,res5[,2]<1e-4)

par(mfrow=c(2,2))
for(nn in names(datasets)){
  hist(diff_res[[nn]][,2],main=nn,xlab = "t-test p-value",xlim=c(0,1))
}
dev.off()
par(mfrow=c(2,2))
for(nn in names(datasets)){
  v1 = -log(diff_res[[nn]][,2])
  v2 = -log(runif(10000))
  qqplot(x=v2,y=v1)
  abline(0,1,lwd=2,lty=2)
}
dev.off()

# compare tissues
x1 = diff_res[[3]][,1]
x2 = diff_res[[4]][,1]
inds = intersect(names(x1),names(x2))
x1 = x1[inds];x2=x2[inds]
plot(x1,x2)

all_pvals = unlist(c(sapply(diff_res,function(x)x[,2])))
hist(all_pvals)
all_qvals = p.adjust(all_pvals,method="BY")
thr = max(all_pvals[all_qvals < 0.1])
passed_results = lapply(diff_res,function(x,y)x[x[,2]<=y,],y=thr)
sapply(passed_results,dim)




