#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

trans_df  = read.delim(args[1], row.names=1, sep="\t", check.names = F)
genes_expr= read.delim(args[2], sep="\t", header=F)
t2g_db   = read.delim(args[3])
output_fn = args[4]

transcript_filters <- function(input, genes_expr, t2g_db) {
    input = merge(t2g_db, input, by.x=1, by.y=0)
    colnames(input) =gsub("_1.fastq","",colnames(input))
    k=input[input$Gene %in% genes_expr$V1,]
    k1=k[c(3:ncol(k))]
    k2=t(t(k1)/(colSums(k1)/1000000))
    k2=round(k2,digits = 5)
    k2=cbind(k[c(1:2)],k2)
    s=aggregate(cbind(k2[c(3:ncol(k2))]), list(k2$Gene), FUN=sum)
    return(s)
}

filtered_df = transcript_filters(trans_df, genes_expr, t2g_db)
write.table(filtered_df, file=output_fn, quote = F, sep="\t", row.names = F)
