# for ballgown analysis
library(ballgown)
# for plot labeling
library(calibrate)
# for manipulating file names/strings
library(stringr) 

# make the ballgown object:
bg = ballgown (dataDir="./results/stringtie/" ,samplePattern="hcc", meas='all') 

## where did you merge from
bg@dirs

## when did you merge?
bg@mergedDate

## what did you import?
bg@meas

## ----getexpr-----
transcript_fpkm = texpr(bg, 'FPKM')
transcript_cov = texpr(bg, 'cov')
whole_tx_table = texpr(bg, 'all')
exon_mcov = eexpr(bg, 'mcov')
junction_rcount = iexpr(bg)
whole_intron_table = iexpr(bg, 'all')
gene_expression = gexpr(bg)


## ----struct-------
structure(bg)$exon
structure(bg)$intron
structure(bg)$trans

exon_transcript_table = indexes(bg)$e2t
transcript_gene_table = indexes(bg)$t2g

# View(head(transcript_gene_table))
# View(head(exon_transcript_table))

# how many exons are there?
length(rownames(exon_transcript_table))

# how many transcripts are there?
length(rownames(transcript_gene_table))

# how many genes are there?
length(unique(transcript_gene_table[,"g_id"])) #Unique Gene count

### transcript stats
## plot average transcript length
hist(whole_tx_table$length, breaks=50, xlab="Transcript length (bp)", main="Distribution of transcript lengths", col="steelblue")

# how many transcripts are there per gene? count the number of genes and count the number of transcripts pere gene and plot it.
counts=table(transcript_gene_table[,"g_id"])
# View(counts)

## some interesting stats
# genes with one transcript
c_one = length(which(counts == 1))

# genes with more than one transcript
c_more_than_one = length(which(counts > 1))

# what is the maximum number of transcript per gene
c_max = max(counts)

# plot above data
hist(counts, breaks=50, col="bisque4", xlab="Transcripts per gene", main="Distribution of transcript count per gene")

# add legend
legend_text = c(paste("Genes with one transcript =", c_one), paste("Genes with more than one transcript =", c_more_than_one), paste("Max transcripts for single gene = ", c_max))
legend("topright", legend_text)

## extract gene names and transcript names
gene_names=data.frame(SYMBOL=unique(rownames(gene_expression)))

#View(gene_names)
t_names=unique(whole_tx_table[,c(1,6)])

#View(whole_tx_table)

## sample meta data---
phenotype_table= data.frame(id=sampleNames(bg), group=rep(c("normal","tumor"), each=2))
pData(bg) =phenotype_table
phenotype_table

## plotTranscripts for a gene, for a sample
plotTranscripts(gene='NCF4', bg, samples='hcc1395_normal_rep1', meas='FPKM', colorby='transcript', main='transcripts from gene NCF4: hcc1395_normal_rep1, FPKM')

## plotTranscripts for a gene, for 3 samples
plotTranscripts('NCF4', bg,  samples=c('hcc1395_normal_rep1', 'hcc1395_normal_rep2', 'hcc1395_tumor_rep1', 'hcc1395_tumor_rep2'),  meas='FPKM', colorby='transcript')

## plot transcript means for all the samples, 
plotMeans('NCF4', bg, groupvar='group', meas='FPKM', colorby='transcript')

### boxplot with and without log transformation
par(mfrow=c(1,2))

boxplot(gene_expression, col=rainbow(4),  las=2, ylab="log2(FPKM)", main="Distribution of FPKMs for all 6 samples")

boxplot(log2(gene_expression+1), col=rainbow(6),  las=2, ylab="log2(FPKM)", main="log transformed distribution of FPKMs for all 6 samples")
# dev.off()

## differential transcript expression
results_txns = stattest(bg, feature='transcript', getFC = T, covariate='group',meas='FPKM' )

# Extract transcript names
t.ids=whole_tx_table[,c(1,6)]
head(results_txns)

# merge transcript results with transcript names
results_txns_merged = merge(results_txns,t.ids,by.x=c("id"),by.y=c("t_id"))
head(results_txns_merged)

# Differential gene expression
results_genes = stattest(bg, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
#View(head(results_genes))

## histogram of diffrentially expressed genes
# Log fold changes and store it in logfc column
results_genes[,"logfc"] = log2(results_genes[,"fc"])

# Fitler results by significant pvalue
sig=which(results_genes$qval<0.05)
#View(sig)
## View(results_genes[sig,])

# draw histogram
hist(results_genes[sig,"logfc"], breaks=50, col="seagreen", xlab="log2(Fold change) Tumor vs Normal", main="Distribution of differential expression values")

# Add vertical cut offs
abline(v=c(-2,2), col="black", lwd=2, lty=2)

# Add legend
legend("topright", "Fold change >2 and <-2", lwd=2, lty=2)

## correlation plot between tumor and normal samples. Average expression of Normal samples Vs Average expression of Tumor samples

# Convert the matrix to data
gene_expression=as.data.frame(gene_expression)    
## View(gene_expression)

# create normal means column
gene_expression$normal=rowMeans(gene_expression[,c(1:2)])

# create tumor means column
gene_expression$tumor=rowMeans(gene_expression[,c(3:4)])

#write.table(gene_expression, "gene_expression.txt", sep="\t")

# to avoid log 0, add 1 to log values. FPKM values are not normalized
x=log2(gene_expression[,"normal"]+1)
y=log2(gene_expression[,"tumor"]+1)
plot(x=x, y=y, pch=1, cex=2, xlab="Normal FPKM (log2)", ylab="Tumor (log2)", main="Tumor vs Normal FPKMs")
abline(a=0, b=1)

# qval significance
qsig=which(results_genes$qval<0.05)
xqsig=x[qsig]
yqsig=y[qsig]
points(x=xqsig, y=yqsig, col="green", pch=19, cex=2)

## fold change signiificance
fsig=which(abs(results_genes$logfc)>4)
xfsig=x[fsig]
yfsig=y[fsig]
points(x=xfsig, y=yfsig, col="red", pch=1, cex=2)

## legend
legend_text = c("Significant by Q value", "Significant by Fold change")
legend("topright", legend_text,bty="n",pch = c(19,19), col=c("green","red"))

# label the significant genes
textxy(xfsig,yfsig, cex=0.8, labs=row.names(gene_expression[fsig,]))


# add red line through 0
abline(v=0, col="red", lwd=3)
# add red line through fold change 4 (log2,2)
abline(v=c(4,-4), col="red", lwd=3)
abline(h=c(-4,4), col="red",lwd=3)

## volcano plot 
# filter by log fold change by 16 fold
fc_sig_results_genes=which(abs(results_genes$logfc)>4)

# Extract genes with fold change by 16 fold
fc_sig_results_genes_plot=results_genes[fc_sig_results_genes,]  

# plot
plot(results_genes$logfc,results_genes$qval, col="steelblue", pch=1) 

#abline
abline(v=c(2,-2), col="red", lwd=3)
abline(h=0.05, col="red",lwd=3)

# highlight the genes with color  
points(fc_sig_results_genes_plot$logfc,fc_sig_results_genes_plot$qval, col="green", pch=16) 

# label the significant genes
textxy(fc_sig_results_genes_plot$logfc,fc_sig_results_genes_plot$qval, labs=fc_sig_results_genes_plot$id, cex=1.2)

## density plot of differentially expressed genes
colors = colorRampPalette(c("white", "blue","red","green","yellow"))
par(mfrow=c(1,2))
plot(x,y, main="Scatter plot of DEGs")
smoothScatter(x,y, colramp = colors, main="Density plot of DEGs")

## write the results to the file
# Filter results_genes by p-value significance
sigpi = which(results_genes[,"pval"]<0.05)

# Extract p-significant genes in a separate object
sigp = results_genes[sigpi,]
## View(sigp)

# filter fc significant genes from p significant genes in a separate object
sigde = which(abs(sigp[,"logfc"]) >= 2)

# extract fc significant genes from p significant genes in a separate object
sig_tn_de = sigp[sigde,]
#View(sig_tn_de)

# Order by q value, followed by differential expression
o = order(sig_tn_de[,"qval"], -abs(sig_tn_de[,"logfc"]), decreasing=FALSE)

# write output to local disc with columns of desired output
output = sig_tn_de[o,c("id","fc","pval","qval","logfc")]

write.table(output, file="./results/ballgown/SigDE.txt", sep="\t", row.names=FALSE, quote=FALSE)


#View(gene_expression)
## heatmap
# Extract gene expression values using significant genes
dim(sig_tn_de)
#View(sig_tn_de)
dim(gene_expression)
# View(gene_expression)
length(sig_tn_de$id)
sig_gene_expression=gene_expression[rownames(gene_expression) %in% sig_tn_de$id,]
dim(sig_gene_expression)
#View(sig_gene_expression)

#remove tumor and normal columns
sig_gene_expression=sig_gene_expression[,-c(5:6)]
phenotype_table
# for pheatmap function, column names and row names of data and pdata mush be identical
# change the row names
rownames(phenotype_table)=phenotype_table[,1]

# remove the id column
phenotype_table=subset(phenotype_table, select = -c(id) )
# change the colnames to match with the sample names
colnames(sig_gene_expression)=row.names(phenotype_table)

# draw heatmap
library(pheatmap)
pheatmap(as.matrix(sig_gene_expression), scale = "row", clustering_distance_rows = "correlation", clustering_method = "complete",annotation_col = phenotype_table , main="Significant genes",fontsize_col=14, fontsize_row = 6 ,color = c("green","red"))

## Draw PCA plot
# transpose the data and compute principal components
pca_data=prcomp(t(sig_gene_expression))

# Calculate PCA component percentages
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)

print(pca_data_perc)

# Extract 1 and 2 principle components and create a data frame with sample names, first and second principal components and group information
df_pca_data = data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], sample = colnames(sig_gene_expression), condition = rep(c("Normal","Tumor"),each=2))
## View(df_pca_data)

## use ggplot package to draw
# color by sample
library(ggplot2)
library(ggrepel)
ggplot(df_pca_data, aes(PC1,PC2, color = sample))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")")) 

# color by condition/group
ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
  geom_text_repel(aes(label=sample),point.padding = 0.75)


# calculate the variance by each pc. Formula is variance/overall variance. Multiply this by 100 and round the result to 2 points


## gene wise pca
temp_data=prcomp(sig_gene_expression)
temp_data_df=data.frame(x=temp_data$x[,1], y=temp_data$x[,2])
ggplot(temp_data_df, aes(x,y))+geom_point()+geom_text(label=rownames(temp_data_df))

# save the workplace

# save.image("bg_08012018.rda")
# load("~/example_data/practice_rnaseq_data/ballgown_scripts/bg_08012018.rda")
 
