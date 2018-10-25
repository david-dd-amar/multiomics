setwd("/Users/David/Desktop/multiomics/oren_tz/mouse_als/")
list.files()
# Files in this dir
# combined.gene.sf.tpm - gene expression matrix
# MetabolomicRAW data OD6-9_Bart_Report_MCF2332.xlsx- targeted metabolomics, samples x metabolites
# Metabolomics.csv - targeted metabolomics, metabolites x samples
# MiceGeneConvertion.R - a script that uses biomart for id conversion
# Proteomics_massaged.csv - proteomics proteins x samples
# RNA-hFUS_Mice_stranded.csv - sample metadata, samples x features
# SamplesToMatch_ER.xlsx - omics-specific metadata, need to examine manually
# tx2gene.csv - ensembl transcript to gene mapping

abundance_data = list()
abundance_data[["ge"]] = read.table("combined.gene.sf.tpm",header = T,row.names = 1,stringsAsFactors = F)
abundance_data[["metab"]] = read.csv("Metabolomics.csv",header = T,row.names = 1,stringsAsFactors = F)
abundance_data[["prot"]] = read.csv("Proteomics_massaged.csv",header = T,row.names = 2,stringsAsFactors = F)
metadata = list()
metadata[["ge"]] = read.csv("RNA-hFUS_Mice_stranded.csv",row.names = 1,header = T)
metadata[["metab"]] = read.table("meta_metabolomics.txt",row.names = 1,header = T,sep="\t")
metadata[["prot"]] = read.table("meta_proteomics.txt",row.names = 1,header = T,sep="\t")

# Some stats
sapply(metadata,dim)
lapply(metadata, function(x)apply(x,2,table))

# format the abundance data into numeric matrices or data frames
# and put the data in log scale
abundance_data[["prot"]] = abundance_data[["prot"]][,-c(1:2)]
abundance_data[["prot"]][,1:10] = log(abundance_data[["prot"]][,1:10]+1,base=2)
abundance_data[["ge"]] = log(abundance_data[["ge"]]+1,base=2)

# Some basic stats
sapply(abundance_data,max)
sapply(abundance_data,min)
sapply(abundance_data,function(x)sum(x==0))
sapply(abundance_data,function(x)quantile(cor(t(x[1:2000,]),method = "spearman"),na.rm = T))
sapply(abundance_data,dim)
sapply(abundance_data,function(x)quantile(x[,1]))

par(mfrow = c(2,2))
boxplot(abundance_data[["ge"]],las=2,main = "Transcriptomics",names=NULL)
boxplot(abundance_data[["prot"]][,11:19],las=2,main="Proteomics, LFQ",names=NULL)
boxplot(abundance_data[["prot"]][,2:10],main="Proteomics, iBAQ",names=NULL)
boxplot(abundance_data[["metab"]],las=2,main = "Metabolomics",names=NULL)
colnames(abundance_data[["prot"]])
metadata[["prot"]]
library(corrplot)
corrplot(cor(abundance_data[["prot"]][,11:19],abundance_data[["prot"]][,2:10],method = "spearman"))
dev.off()

exclude_analyte_by_thrs<-function(x,intensity_thr = 3, percent_thr = 0.8){
  return(sum(x<intensity_thr)/length(x) >= percent_thr)
}
# Filter out unreliable analytes from GE and proteomics
to_rem = apply(abundance_data[["ge"]],1,exclude_analyte_by_thrs,intensity_thr = 2) |
  apply(abundance_data[["ge"]]==0,1,any)
table(to_rem)
abundance_data[["ge"]] = abundance_data[["ge"]][!to_rem,]
to_rem = apply(abundance_data[["prot"]][,11:19],1,exclude_analyte_by_thrs,intensity_thr = 24)
table(to_rem)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Helper functions
select_analyte_pval_fchange<-function(x,p,f){
  n = length(x)
  fs = x[1:(n/2)]
  ps = x[((n/2)+1):n]
  return(any(abs(fs)>=f & ps <= p))
}
get_diff_abundance_raw_results<-function(lm_diff_res,selected_analytes){
  m = c()
  for(nn in names(lm_diff_res)){
    curr_res = lm_diff_res[[nn]][,selected_analytes[[nn]]]
    curr_res = t(curr_res)
    curr_res = cbind(curr_res,rep(nn,nrow(curr_res)))
    m = rbind(m,curr_res)
  }
  dim(m)
  m = cbind(rownames(m),m)
  colnames(m)[1] = "analyte"
  return(m)
}
# Assumes that many to many mappings is not a big issue
transform_raw_results_to_gene_based<-function(m,uniprot2ensemble_pro){
  curr_prot = m[m[,6]=="prot",1]
  curr_prot2ensembl_prot = uniprot2ensemble_pro[is.element(uniprot2ensemble_pro[,1],set=curr_prot),]
  curr_prot2ensembl_prot[curr_prot2ensembl_prot[,1]=="P97372",]
  curr_prot2ensembl_gene = uniprot2ensemble_gene[is.element(uniprot2ensemble_gene[,1],set=curr_prot),]
  curr_prot2ensembl_gene = curr_prot2ensembl_gene[is.element(curr_prot2ensembl_gene[,3],set=rownames(abundance_data[["ge"]])),]
  # curr_prot2ensembl_gene[curr_prot2ensembl_gene[,1]=="Q80WG5",]
  # lm_diff_res[["ge"]][,"ENSMUSG00000007476"]
  rough_intersection_pval = phyper(31,500,3000,116,lower.tail = F)
  all_our_genes = union(m[m[,6]=="ge",1],curr_prot2ensembl_gene[,3])
  new_m = c()
  count = 0
  for(g in all_our_genes){
    v1 = rep(NA,nrow(lm_diff_res[["ge"]]))
    if(is.element(g,set=colnames(lm_diff_res[["ge"]]))){
      v1 = lm_diff_res[["ge"]][,g] 
    }
    curr_prot = all_prot2ensembl_gene[all_prot2ensembl_gene[,3]==g,1]
    count = count + length(curr_prot)
    v2 = rep(NA,nrow(lm_diff_res[["ge"]]))
    if(length(curr_prot)>0){
      v2 = lm_diff_res[["prot"]][,curr_prot[1]] 
    }
    new_m = rbind(new_m,c(v1,v2,curr_prot[1]))
    rownames(new_m)[nrow(new_m)] = g
  }
  new_m = cbind(rownames(new_m),new_m)
  colnames(new_m)[1] = "analyte"
  return(list(new_m,rough_intersection_pval))
}
# Read protein mapping info: we want uniprot to ensembl
# Downloaded from uniprot
raw_proteomics = read.delim("MOUSE_10090_idmapping.dat",header = F,stringsAsFactors = F)
uniprot2ensemble_pro = raw_proteomics[raw_proteomics[,2] == "Ensembl_PRO",]
uniprot2ensemble_gene = raw_proteomics[raw_proteomics[,2] == "Ensembl",]
all_prot2ensembl_gene = uniprot2ensemble_gene[is.element(uniprot2ensemble_gene[,1],
                                                         set=rownames(abundance_data[["prot"]])),]
# clustering, enrichment, and intepretation of results
library('gskb')
data(mm_pathway)
gene_info = read.table("Mus_musculus.GRCm38.94.gtf_genes_only.txt",stringsAsFactors = F)
gene_names = gene_info[,16]
names(gene_names) = gene_info[,10]
gene_names[1:10]
is.element("Apoe",set=gene_names)
gene_names = toupper(gene_names)
bg = gene_names[rownames(abundance_data[["ge"]])]
simple_enrichment_analysis<-function(set1,set2,bg){
  set1 = intersect(set1,bg)
  set2 = intersect(set2,bg)
  if(length(set1)==0 || length(set2)==0){return(1)}
  x1 = is.element(bg,set=set1)
  x2 = is.element(bg,set=set2)
  tb = table(x1,x2)
  p = fisher.test(tb,alternative = "g")$p.value
  genes = intersect(set1,set2)
  return(list(p=p,genes=genes))
}
enrichment_list_to_table<-function(l){
  pvsl = sapply(l,function(x)x[[1]])
  genes = sapply(l,function(x)paste(x[[2]],collapse=","))
  return(cbind(names(l),pvsl,genes))
}

# Differential abundance using simple linear regression
lm_get_effects_and_pvals<-function(y,x){
  m = lm(y~.,data=data.frame(y,x))
  m = summary(m)$coefficients
  v = c(m[-1,1],m[-1,4])
  names(v)[1:(nrow(m)-1)] = paste("effect_",names(v)[1:(nrow(m)-1)],sep="")
  names(v)[nrow(m):length(v)] = paste("pval_", names(v)[nrow(m):length(v)],sep="")
  return(v)
}

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


lm_diff_res = list()
# metabolomics
# metadata rows fit the abundance data column: manual examination
ys = as.matrix(abundance_data[["metab"]])
x = metadata[["metab"]][,2:3]
lm_diff_res[["metab"]] = apply(ys,1,lm_get_effects_and_pvals,x=x)
# gene expression
ys = as.matrix(abundance_data[["ge"]])
x = metadata[["ge"]][,2:3]
all(colnames(ys)==rownames(x))
all(is.element(colnames(ys),set=rownames(x)))
x = x[colnames(ys),] # correct metadata order
lm_diff_res[["ge"]] = apply(ys,1,lm_get_effects_and_pvals,x=x)
dim(lm_diff_res[["ge"]])
hist(lm_diff_res[["ge"]][2,])
# Proteomics
ys = as.matrix(abundance_data[["prot"]][,2:10])
x = metadata[["prot"]][,1:2]
col_ys_info = sapply(colnames(ys),function(x)strsplit(x,split="\\.")[[1]])
colnames(ys) = col_ys_info[2,]
x = x[colnames(ys),] # correct metadata order
lm_diff_res[["prot"]] = apply(ys,1,lm_get_effects_and_pvals,x=x)
par(mfrow=c(1,2))
hist(lm_diff_res[["prot"]][2,])
hist(lm_diff_res[["prot"]][1,])
dev.off()

# look at p-values for treatments
par(mfrow=c(3,1))
hist(lm_diff_res[["ge"]][4,],main="Transcriptomics")
hist(lm_diff_res[["metab"]][4,],main="Metabolomics")
hist(lm_diff_res[["prot"]][4,],main = "Proteomics")
# qqplot(-log(runif(1000)),-log(lm_diff_res[["ge"]][4,]));abline(0,1)
# qqplot(-log(runif(1000)),-log(lm_diff_res[["prot"]][4,]));abline(0,1)
# qqplot(-log(runif(1000)),-log(lm_diff_res[["metab"]][4,]));abline(0,1)
dev.off()

# Some plots of fold change data
x1 = lm_diff_res[["ge"]][1,]
x2 = -lm_diff_res[["ge"]][2,]
plot(x=x1,y=x2,main=paste("Transcriptomics, rho=",round(cor(x1,x2),digits = 4)),
     ylab = "Effect of treatment (log fchange)",
     xlab = "Effect of genotype (log fchange)",pch=20,col="blue")
m = lm(x2~x1)
abline(m[[1]][1],m[[1]][2],lwd=2,lty=2)

x1 = lm_diff_res[["prot"]][1,]
x2 = -lm_diff_res[["prot"]][2,]
plot(x=x1,y=x2,main=paste("Proteomics, rho=",round(cor(x1,x2),digits = 4)),
     ylab = "Effect of treatment (log fchange)",
     xlab = "Effect of genotype (log fchange)",pch=20,col="blue")
m = lm(x2~x1)
abline(m[[1]][1],m[[1]][2],lwd=2,lty=2)

x1 = lm_diff_res[["metab"]][1,]
x2 = -lm_diff_res[["metab"]][2,]
plot(x=x1,y=x2,main=paste("Metabolomics, rho=",round(cor(x1,x2),digits = 4)),
     ylab = "Effect of treatment (log fchange)",
     xlab = "Effect of genotype (log fchange)",pch=20,col="blue")
m = lm(x2~x1)
abline(m[[1]][1],m[[1]][2],lwd=2,lty=2)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# Selection and interpretation

sapply(lm_diff_res,dim)
ps = sapply(lm_diff_res,function(x)c(x[3:4,]))
ps_v = unlist(ps)
table(p.adjust(ps_v,method="fdr")<0.01)
qs_v = p.adjust(ps_v,method="fdr")
pval_thr = max(ps_v[qs_v < 0.05]) 

# Analysis for a specific fchange analysis
selected_analytes = lapply(lm_diff_res,function(x,...)
  apply(x,2,select_analyte_pval_fchange,...),f=0.5,p=pval_thr)
sapply(selected_analytes,table)
# print the output
m = get_diff_abundance_raw_results(lm_diff_res,selected_analytes)
write.table(m,file="analyte_diff_analysis_results.txt",sep="\t",quote = T,col.names = T,row.names = F)
new_m = transform_raw_results_to_gene_based(m,uniprot2ensemble_pro)
write.table(new_m[[1]],file="analyte_diff_analysis_results_by_gene.txt",sep="\t",quote = T,col.names = T,row.names = F)

# Check proteomics-transcriptomics agreement
new_m = new_m[[1]]
i1 = 3
i2 = 7
inds = !is.na(new_m[,i1]) & !is.na(new_m[,i2])
x1 = as.numeric(new_m[inds,i1])
x2 = as.numeric(new_m[inds,i2])
plot(x=x1,y=x2,main=paste("Treatment effect, rho=",round(cor(x1,x2,method="spearman"),digits = 4)),
     ylab = "Transcriptomics (log fchange)",
     xlab = "Proteomics (log fchange)",pch=20,col="blue")
m = lm(x2~x1)
abline(m[[1]][1],m[[1]][2],lwd=2,lty=2)
i1 = 2
i2 = 6
inds = !is.na(new_m[,i1]) & !is.na(new_m[,i2])
x1 = as.numeric(new_m[inds,i1])
x2 = as.numeric(new_m[inds,i2])
plot(x=x1,y=x2,main=paste("FUS effect, rho=",round(cor(x1,x2,method="spearman"),digits = 4)),
     ylab = "Transcriptomics (log fchange)",
     xlab = "Proteomics (log fchange)",pch=20,col="blue")
m = lm(x2~x1)
abline(m[[1]][1],m[[1]][2],lwd=2,lty=2)

# partition the genes into groups and run enrichment
new_m_effects = new_m[,c(2:3)]
mode(new_m_effects) = "numeric"

kmeans_sol = kmeans(new_m_effects,centers = 2)
table(kmeans_sol$cluster)
boxplot(new_m_effects[,1]~kmeans_sol$cluster)
boxplot(new_m_effects[,2]~kmeans_sol$cluster)

set1 = gene_names[rownames(new_m)[kmeans_sol$cluster==2]]
enrichment_res = sapply(mm_pathway,simple_enrichment_analysis,set2=set1,bg=bg)
enrichment_pvals = sapply(enrichment_res,function(x)x[1])
# write.table(enrichment_list_to_table(enrichment_res[(p.adjust(unlist(enrichment_pvals),method="fdr")<0.1)]),
#             file = "simple_enrichment_results.txt",sep="\t",row.names = F)
enrichment_list_to_table(enrichment_res[(p.adjust(unlist(enrichment_pvals),method="fdr")<0.05)])



