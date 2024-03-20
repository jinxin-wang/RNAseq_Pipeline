#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trans_df = read.delim(args[1], row.names=1, sep="\t", check.names = F)
t2g_db   = read.delim(args[2])
output_fn= args[3]

transcripts_to_genes <- function(input, t2g_db) {
    colnames(input) = gsub("_1.fastq", "", colnames(input))
    input = merge(t2g_db, input, by.x=1, by.y=0)
    res = aggregate(cbind(input[c(3:ncol(input))]),list(input$Gene), FUN = sum)
    return(res)
}


genes_df = transcripts_to_genes(trans_df, t2g_db)
write.table(genes_df, file=output_fn, quote = F, sep="\t", row.names = F)

