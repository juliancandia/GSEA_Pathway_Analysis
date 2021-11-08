rm(list=ls())

PROJECT_DIR = "~/git/GSEA_Pathway_Analysis" # replace with your local path

infile = file.path(PROJECT_DIR,"REF","HGNC","genenames.txt")
data = as.matrix(read.table(infile,comment.char="",quote="",header=T,sep="\t"))
data[is.na(data[,"Approved.symbol"]),"Approved.symbol"] = ""

res = c("HGNC.ID","ENSG","Symbol")
for (i in 1:nrow(data)) {
    uid = data[i,"HGNC.ID"] # unique identifier
    ensg1 = data[i,"Ensembl.gene.ID"]
    ensg2 = data[i,"Ensembl.ID.supplied.by.Ensembl."]
    if (ensg1!="") {
        ensg = ensg1
    } else {
        ensg = ensg2
    }
    symbol = data[i,"Approved.symbol"]
    symbol_prev = trimws(unlist(strsplit(data[i,"Previous.symbols"],",")))
    symbol_alias = trimws(unlist(strsplit(data[i,"Alias.symbols"],",")))
    symbol = unique(c(symbol,symbol_prev,symbol_alias))
    n_symbol = length(symbol)
    if (n_symbol==1) {
        res = rbind(res,c(uid,ensg,symbol))
    } else if (n_symbol>1) {
        res = rbind(res,cbind(rep(uid,n_symbol),rep(ensg,n_symbol),symbol))
    }
}

output = res
outfile = file.path(PROJECT_DIR,"REF","HGNC","genenames_map.txt")
write(t(output),ncol=ncol(output),file=outfile,sep="\t")
