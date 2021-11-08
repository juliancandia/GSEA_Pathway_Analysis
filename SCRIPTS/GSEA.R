#############################################################
# This script belongs to the GSEA PAthway Analysis Pipeline #
# written by Julián Candia, julian.candia@nih.gov (c) 2021  #
#############################################################

library(fgsea)
library(dplyr)

rm(list=ls())

PROJECT_DIR = "~/git/GSEA_Pathway_Analysis" # replace with your local path

# list of all gene sets available (do NOT modify!)
ref_label_all = c("Hallmarks","Canonical_Pathways","Chem_Genet_Perturb","miR_Targets","TF_Targets","GO_BP","GO_CC","GO_MF","Immune","Cell_Type")
ref_file_all = c("h.all","c2.cp","c2.cgp","c3.mir","c3.tft","c5.go.bp","c5.go.cc","c5.go.mf","c7.immunesigdb","c8.all")

##### HERE: provide info ######################
input_filename = "DEG_topTable.txt" # input file name
ref_label = c("Hallmarks","Canonical_Pathways") # gene sets to run (must be a subset of "ref_label_all")

# This code generates n_run GSEA runs and finds consensus solutions (because GSEA is a stochastic process).
#n_run = 1000 # number of GSEA runs
n_run = 10 # for testing purposes

# Parameters to select pathways and leading edge genes (can be tweaked but you don't have to).
adj_pval_thres = 0.05 # adjusted p-value threshold for pathway selection
pwy_sel_thres = 80 # percent threshold across runs for pathway selection
gene_sel_thres = 80 # percent threshold across runs for gene selection
seed = 123 # fixed random number seed to ensure reproducibility
eps = 1e-10 # boundary for calculating pathway p-values used by fgseaMultilevel. Can be set to zero for better (but slower) estimation.
###############################################

infile = file.path(PROJECT_DIR,"REF","HGNC","genenames_map.txt")
map = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))

# Read the input list of genes from the case vs control differential expression analysis. Create a list ordered by significance based on -log10(p-value) such that the top genes are the most significantly overexpressed in cases, and the bottom genes are the most significantly overexpressed in controls.
infile = file.path(PROJECT_DIR,"DATA",input_filename)
DEG = as.matrix(read.table(infile,stringsAsFactors=F,header=T,sep="\t",quote="",comment.char=""))
# we remove duplicated names (if any)
if (any(duplicated(DEG[,"name"]))) {
  DEG = DEG[!duplicated(DEG[,"name"]),]
}
gene_list = (-log10(as.numeric(DEG[,"P.Value"])))*sign(as.numeric(DEG[,"logFC"]))
o = order(-gene_list)
DEG = DEG[o,]
gene_list = gene_list[o]

ref_file = NULL
for (geneset in ref_label) {
    ref_file = c(ref_file,ref_file_all[ref_label_all==geneset])
}
n_ref = length(ref_file)
for (i_ref in 1:n_ref) {
    ref_pathways = fgsea::gmtPathways(file.path(PROJECT_DIR,"REF","MSigDB",paste0(ref_file[i_ref],".v7.4.symbols.gmt")))
    # gene relabeled as needed
    ref_genes = NULL
    for (i_pwy in 1:length(ref_pathways)) {
        ref_genes = unique(c(ref_genes,ref_pathways[[i_pwy]]))
    }
    gene_names = DEG[,"name"]
    sel_missing = !gene_names%in%ref_genes
    gene_missing = DEG[sel_missing,c("ENSG","name")]
    n_missing_pre = nrow(gene_missing)
    for (i in 1:nrow(gene_missing)) {
        symbol = map[map[,"ENSG"]==gene_missing[i,"ENSG"],"Symbol"]
        sel = symbol%in%ref_genes
        if (sum(sel)>0) {
            gene_missing[i,"name"] = symbol[sel][1]
        }
    }
    gene_names[sel_missing] = gene_missing[,"name"]
    n_missing_post = sum(!gene_missing[,"name"]%in%ref_genes)
    cat(paste0("Number of gene names relabeled for ",ref_label[i_ref]," analysis:"),n_missing_pre-n_missing_post,"\n")
    names(gene_list) = gene_names
    
    set.seed(seed)
    fgRes = vector("list",n_run)
    for (i_run in 1:n_run) {
        fgRes[[i_run]] = fgseaMultilevel(pathways=ref_pathways,stats=gene_list,nPermSimple=10000,eps=eps) %>% as.data.frame() %>% dplyr::filter(padj<!!adj_pval_thres)
    }
    concat_pwy = NULL
    for (i_run in 1:n_run) {
        concat_pwy = c(concat_pwy,fgRes[[i_run]][,"pathway"])
    }
    table_pwy = 100*table(concat_pwy)/n_run
    table_pwy = table_pwy[table_pwy>pwy_sel_thres]
    pwy = names(table_pwy)
    n_pwy = length(pwy)
    pwy_genes = vector("list",n_pwy)
    pwy_pval_med = rep(NA,n_pwy)
    pwy_padj_med = rep(NA,n_pwy)
    pwy_ES_med = rep(NA,n_pwy)
    pwy_dir = rep(NA,n_pwy)
    pwy_size = rep(NA,n_pwy)
    for (i_pwy in 1:n_pwy) {
        concat_genes = NULL
        pwy_pval = NULL
        pwy_padj = NULL
        pwy_ES = NULL
        for (i_run in 1:n_run) {
            if (pwy[i_pwy]%in%fgRes[[i_run]][,"pathway"]) {
                pwy_index = which(fgRes[[i_run]][,"pathway"]==pwy[i_pwy])
                concat_genes = c(concat_genes,fgRes[[i_run]]$leadingEdge[pwy_index][[1]])
                pwy_pval = c(pwy_pval,fgRes[[i_run]][pwy_index,"pval"])
                pwy_padj = c(pwy_padj,fgRes[[i_run]][pwy_index,"padj"])
                pwy_ES = c(pwy_ES,fgRes[[i_run]][pwy_index,"ES"])
                pwy_size[i_pwy] = fgRes[[i_run]][pwy_index,"size"] # I can rewrite size across runs; it's a constant value
            }
        }
        table_genes = table(concat_genes)
        pwy_genes[[i_pwy]] = names(table_genes[table_genes>(gene_sel_thres*length(pwy_pval)/100)])
        pwy_pval_med[i_pwy] = median(pwy_pval)
        pwy_padj_med[i_pwy] = median(pwy_padj)
        pwy_ES_med[i_pwy] = median(pwy_ES)
        pwy_dir[i_pwy] = mean(sign(pwy_ES))
        if (abs(pwy_dir[i_pwy])<1) {
            pwy_dir[i_pwy] = 0 # remove pathways whose direction is not consistent
        }
    }
    o = order(pwy_padj_med)
    pwy_index = c(o[which(pwy_dir[o]==1)],o[which(pwy_dir[o]==(-1))]) # NOTE: pathways ordered by direction (first pwy_dir=1 i.e. cases, followed by pwy_dir=-1 i.e. controls), then by adjusted p-value within each group. If the direction is not always consistent, pwy_dir=0 and the pathway is not considered.
    
    # we save the significant pathways ordered by p-value within each enrichment group
    pwy_dir_label = rep("Cases",n_pwy)
    pwy_dir_label[pwy_dir==(-1)] = "Controls"
    leading_edge_genes = NULL
    for (i_pwy in pwy_index) {
      leading_edge_genes = rbind(leading_edge_genes,c(length(pwy_genes[[i_pwy]]),paste0(pwy_genes[[i_pwy]],collapse=",")))
    }
    output = rbind(c("pathway","n_measured_genes","enriched_class","perc_signif_padj_0.05","ES_median","pval_median","padj_median","n_leading_edge_genes","leading_edge_genes"),cbind(cbind(pwy,pwy_size,pwy_dir_label,table_pwy,pwy_ES_med,pwy_pval_med,pwy_padj_med)[pwy_index,],leading_edge_genes))
    outfile = file.path(PROJECT_DIR,"RESULTS",paste0(ref_label[i_ref],"_n",n_run,".txt"))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
    
    # we generate pathway lists for all the leading-edge genes
    tmp = table(unlist(pwy_genes))
    gene = names(tmp[order(-tmp)])
    n_gene = length(gene)
    gene_pwy = matrix(rep(0,n_gene*n_pwy),ncol=n_pwy)
    for (i_pwy in 1:n_pwy) {
      gene_pwy[gene%in%pwy_genes[[i_pwy]],i_pwy] = 1
    }
    pwys_per_gene = rep("",n_gene)
    n_pwys_per_gene = rep("",n_gene)
    for (i_gene in 1:n_gene) {
      index = which(gene_pwy[i_gene,]==1) 
      pwys_per_gene[i_gene] = paste0(pwy[index],collapse=",")
      n_pwys_per_gene[i_gene] = length(index)
    }
    outfile = file.path(PROJECT_DIR,"RESULTS",paste0(ref_label[i_ref],"_n",n_run,"_pwys_per_gene.txt"))
    output = rbind(c("leading_edge_gene","n_pathways","pathways"),cbind(gene,n_pwys_per_gene,pwys_per_gene))
    write(t(output),ncol=ncol(output),file=outfile,sep="\t")
}
