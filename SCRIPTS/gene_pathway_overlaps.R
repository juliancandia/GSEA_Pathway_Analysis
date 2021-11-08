#############################################################
# This script belongs to the GSEA PAthway Analysis Pipeline #
# written by Julián Candia, julian.candia@nih.gov (c) 2021  #
#############################################################

library(fgsea)

rm(list=ls())

PROJECT_DIR = "~/git/GSEA_Pathway_Analysis" # replace with your local path

# list of all gene sets available (do NOT modify!)
ref_label_all = c("Hallmarks","Canonical_Pathways","Chem_Genet_Perturb","miR_Targets","TF_Targets","GO_BP","GO_CC","GO_MF","Immune","Cell_Type")
ref_file_all = c("h.all","c2.cp","c2.cgp","c3.mir","c3.tft","c5.go.bp","c5.go.cc","c5.go.mf","c7.immunesigdb","c8.all")

##### HERE: provide info ######################
ref_label = c("Hallmarks","Canonical_Pathways") # gene sets to run (must be a subset of "ref_label_all")
#n_run = 1000 # number of GSEA runs
n_run = 10 # for testing purposes
frac_thres = 0.95 # overlap fraction for pathway removal (must be <=1, closer to 1 is more lenient - i.e. more pathways are kept in the final solution).
###############################################

ref_file = NULL
for (geneset in ref_label) {
    ref_file = c(ref_file,ref_file_all[ref_label_all==geneset])
}
n_ref = length(ref_file)
for (i_ref in 1:n_ref) {
    infile = file.path(PROJECT_DIR,"RESULTS",paste0(ref_label[i_ref],"_n",n_run,".txt"))
    data = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t"))
    pwy = data[,"pathway"]
    n_pwy = length(pwy)
    ref_pathways = fgsea::gmtPathways(file.path(PROJECT_DIR,"REF","MSigDB",paste0(ref_file[i_ref],".v7.4.symbols.gmt")))
    pwy_in_ref = names(ref_pathways)
    pwy_genes_this_ref = vector("list",n_pwy)
    for (i_pwy in 1:n_pwy) {
        index = which(pwy_in_ref==pwy[i_pwy])
        pwy_genes_this_ref[[i_pwy]] = ref_pathways[[index]]
    }
    if (i_ref==1) {
        merged_data = data
        pwy_genes = pwy_genes_this_ref
    } else {
        merged_data = rbind(merged_data,data)
        pwy_genes = append(pwy_genes,pwy_genes_this_ref)
    }
}

pwy = merged_data[,"pathway"]
n_pwy = length(pwy)
n_gene = rep(NA,n_pwy)
for (i_pwy in 1:n_pwy) {
    n_gene[i_pwy] = length(pwy_genes[[i_pwy]])
}

# we find pathway overlaps recursively
exit_code = 0
while (exit_code==0) {
    res = NULL
    for (i_pwy in 1:(n_pwy-1)) {
        for (j_pwy in (i_pwy+1):n_pwy) {
            n_overlap = sum(pwy_genes[[i_pwy]]%in%pwy_genes[[j_pwy]])
            if (n_overlap>0) {
                frac1 = n_overlap/n_gene[i_pwy]
                frac2 = n_overlap/n_gene[j_pwy]
                if (frac1>frac2) {
                    res = rbind(res,c(pwy[i_pwy],n_gene[i_pwy],frac1,pwy[j_pwy],n_gene[j_pwy],frac2))
                } else {
                    res = rbind(res,c(pwy[j_pwy],n_gene[j_pwy],frac2,pwy[i_pwy],n_gene[i_pwy],frac1))
                }
            }
        }
    }
    res_filt = res[as.numeric(res[,3])>=frac_thres,,drop=F]
    pwy_remove = names(which.max(table(res_filt[,1]))) # we remove the pathway with largest number of hits
    n_pwy_remove = length(pwy_remove)
    if (n_pwy_remove>=1) {
        if (n_pwy_remove>1) { # we resolve ties in this order: the pathway with largest maximum overlap fraction; the pathway with smallest size
            max.frac = rep(NA,n_pwy_remove)
            size = rep(NA,n_pwy_remove)
            for (i_pwy_remove in 1:n_pwy_remove) {
                max.frac[i_pwy_remove] = max(as.numeric(res_filt[res_filt[,1]==pwy_remove[i_pwy_remove],3]))
                size[i_pwy_remove] = as.numeric(res_filt[res_filt[,1]==pwy_remove[i_pwy_remove],2][1])
            }
            pwy_remove = order(-max.frac,size)
        }
        index_remove = which(pwy==pwy_remove[1]) # we only remove one pathway in each iteration.
        pwy = pwy[-index_remove]
        n_pwy = n_pwy - 1
        n_gene = n_gene[-index_remove]
        pwy_genes[[index_remove]] = NULL
    } else {
        exit_code = 1
    }
}

# we remove pathways that were filtered out during the overlap process
merged_data = merged_data[merged_data[,"pathway"]%in%pwy,]
output = rbind(colnames(merged_data),merged_data)
outfile = file.path(PROJECT_DIR,"RESULTS",paste0("merged_filt_",
paste(ref_label,collapse="_"),"_n",n_run,"_f",frac_thres,".txt"))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")

# we update the pathway lists for all the leading-edge genes
pwy_lead_genes = vector("list",n_pwy)
for (i_pwy in 1:n_pwy) {
    pwy_lead_genes[[i_pwy]] = as.character(unlist(strsplit(merged_data[i_pwy,"leading_edge_genes"],",")))
}
tmp = table(unlist(pwy_lead_genes))
gene = names(tmp[order(-tmp)])
n_gene = length(gene)
gene_pwy = matrix(rep(0,n_gene*n_pwy),ncol=n_pwy)
for (i_pwy in 1:n_pwy) {
  gene_pwy[gene%in%pwy_lead_genes[[i_pwy]],i_pwy] = 1
}
pwys_per_gene = rep("",n_gene)
n_pwys_per_gene = rep("",n_gene)
for (i_gene in 1:n_gene) {
  index = which(gene_pwy[i_gene,]==1)
  pwys_per_gene[i_gene] = paste0(pwy[index],collapse=",")
  n_pwys_per_gene[i_gene] = length(index)
}
outfile = file.path(PROJECT_DIR,"RESULTS",paste0("merged_filt_",
paste(ref_label,collapse="_"),"_n",n_run,"_f",frac_thres,"_pwys_per_gene.txt"))
output = rbind(c("leading_edge_gene","n_pathways","pathways"),cbind(gene,n_pwys_per_gene,pwys_per_gene))
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
