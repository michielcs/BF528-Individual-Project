ds = read.table("/projectnb/bf528/users/group6/project_2/analysis/cuffdiff_out/gene_exp.diff", head= TRUE, sep='\t')
ds = ds[order(ds$q_value),]

#get top 10 rows and only interested columns
top_diff=ds[1:10,]
top_diff=top_diff[c(3,8:10,12:13)]

histogram = hist(ds$log2.fold_change.,main="Histogram of fold change", xlab = "log2 fold change", breaks=20)

#only significant genes
sig = subset(ds, ds$significant=='yes')
histogram2 = hist(sig$log2.fold_change.,,main="Histogram of fold change for significant genes", xlab = "log2 fold change")

#filter up/down regulated genes
up_reg <- subset(sig$gene, sig$log2.fold_change.>1)
down_reg <- subset(sig$gene,sig$log2.fold_change.<1)

length(ds$gene)
length(up_reg)
length(down_reg)

#write up/down regulated genes into csv
write.csv(up_reg,"up_genes_project2.csv")
write.csv(down_reg,"down_genes_project2.csv")
