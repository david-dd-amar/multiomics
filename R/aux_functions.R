# PARAMs: external data used in the pipelines
mouse_genes_path = "/Users/David/Desktop/multiomics/oren_tz/mouse_als/Mus_musculus.GRCm38.94.gtf_genes_only.txt"
mouse_proteins_path = "/Users/David/Desktop/multiomics/oren_tz/mouse_als/MOUSE_10090_idmapping.dat"

# Read protein mapping info: we want uniprot to ensembl
# Data was downloaded from uniprot
raw_proteomics = read.delim(mouse_proteins_path,header = F,stringsAsFactors = F)
uniprot2ensemble_pro = raw_proteomics[raw_proteomics[,2] == "Ensembl_PRO",]
uniprot2ensemble_gene = raw_proteomics[raw_proteomics[,2] == "Ensembl",]
all_prot2ensembl_gene = uniprot2ensemble_gene[is.element(uniprot2ensemble_gene[,1],
                                                         set=rownames(abundance_data[["prot"]])),]
# for clustering, enrichment, and intepretation of results
library('gskb')
data(mm_pathway)
gene_info = read.table(mouse_genes_path,stringsAsFactors = F)
gene_names = gene_info[,16]
names(gene_names) = gene_info[,10]
gene_names = toupper(gene_names)
# bg = gene_names[rownames(abundance_data[["ge"]])]
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
  ps = c()
  genes = c()
  for(j in 1:length(l)){
    ps = c(ps,l[[j]][[1]])
    currg = ""
    if(length(l[[j]])>1){
      currg = paste(l[[j]][[2]],collapse=",")
    }
    genes = c(genes,currg)
  }
  return(cbind(names(l),ps,genes))
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