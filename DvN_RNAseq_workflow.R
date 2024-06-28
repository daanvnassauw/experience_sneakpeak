#Author:  Daan van Nassauw
#Course:  RNA-seq analysis - workflow


##Example of how to install a Bioconductior library...
##this is done only once per computer.
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")

##activate DESeq2, ggplot2 and pheatmap. 

library("DESeq2")  #for differential expression
library("ggplot2") #for plotting.
library("pheatmap") # for plotting pretty heatmaps
#############

###Read in count matrices.
#DeSeq expects input data that are a count of integer values (with the number of reads, fragments or pairs that are associated to each gene in each sample).
##In this case the data are in two different files, so they have to be combined...

#First transfrom the data into CSV files (this is the easiest thing to do....) Here I copied the files in the directory I am working...  I called themreadcount_May2018

#read in the read count files 
countData <- read.csv("read_count_file.csv", row.names=1)
#explore the object

###description of the experiments 
filenameSample= "sample_description.csv"
colData <- read.csv(filenameSample, row.names=1, stringsAsFac=TRUE)

#explore the object

#count matrix and column data need to be  consisent (same names in same order)
colnames(countData)==rownames(colData)  #check that they are in the same order.Everything should be "TRUE"



#read in gene annotation
annotation=read.table("gene_description.tab", sep="\t",header=TRUE, quote="#", row.names=1)
#explore the object



###Normalizing the data
####Create DESeq2 object.
dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ condition)  #"condition" is the colname in coldata



#normalize the data
dds <- DESeq(dds)



## DATA quality assessment. 

##Variance stabilization transformation
vsd <- vst(dds, blind=FALSE) # Variance stabilizing transformation 
vsdData <- assay(vsd)


#######Dedrograms of distance between samples

##Euclidean distance 
distanceEuclidean=dist(t(vsdData))                    
plot(hclust(distanceEuclidean),hang=-1, main="Distance between samples (Euclidean)", xlab=NA, sub=NA)  #Dendogram based on Euclidean Distance


##distance based on correlation
distanceCorr=as.dist(1-cor(vsdData))

#Dendogram based on Euclidean Distance


#########Heatmap of sample-to sample distance 
##use 
#based on correlation
pheatmap(as.matrix( cor((vsdData))))

###PCA with ALL Genes 
xx <- prcomp(t(as.matrix(assay(vsd))))
summary(xx)

##Ok-ish plot
rownames((xx$x))
mycolors=c("red", "red", "blue", "blue", "black","black")#make sure that they are in the same order as the rownames of xx$x
plot(xx$x[,1], xx$x[,2], col=mycolors, pch=16)
legend("topright",c("Control", "GN160", "GN800"),  col=c("red", "blue", "black"), pch=16)


##more elegant plot using the ggplot library.
z=summary(xx)
PC1 <- as.numeric(xx$x[,1]); PC2 <- as.numeric(xx$x[,2])
Sample <- rownames(xx$x)
Sample <- c("C_2h_1"  ,  "C_2h_2"  ,  "GN160_2h_1", "GN160_2h_2" ,"GN800_2h_1", "GN800_2h_2")
#Here you can put the description in the conditions... (as they are in sample) 
Condition <- colData[colnames(countData),"condition"]
Condition <- c("Control","Control", "GN160", "GN160","GN800", "GN800")
df <-data.frame(PC1, PC2, Sample, Condition)
myplot <- ggplot(df,aes( PC1,PC2, color=Condition), color='black') +
    geom_point(size=3,aes(shape=Condition )) +
    geom_text(data=df, mapping=aes(x=PC1+c(2.5,0,-0.5,-0.5,-5.5,-5.5),y=PC2+c(-0.5,-0.5,+0.5,-0.5, -0.5, +0.5), label=Sample), size=4, show_guide=FALSE) +  #in case you want to see the labels
    xlab(paste0("PC1 (",format(100* z$importance[2,1],digits=2),"%)")) +
    ylab(paste0("PC2 (",format(100* z$importance[2,2],digits=2),"%)"))   +
    theme_bw()+
    theme(axis.text.x = element_text(size=15), axis.text.y = element_text(size=15), text=element_text(size=15))    +
    theme(legend.title=element_blank())

myplot  #to see the plot

#to save the plot as a tiff
filename <- "PCA_all.tiff"  #changing the extension to pdf with create a pdf.
ggsave(myplot, filename=filename, width=7, height=7)


##########Differential expression 

##compare two conditions 
FCTh=1.5 #FC threshold to be used. 
padjTh=0.05  #FDR threshold to be used. 

##compare exposure GN160 vs control (control)
##here "condition" is as in the colData file (with sample information) and then what you want to compare. 
res160=results(dds, contrast=c("condition","GN160_2h", "Control"),alpha=padjTh)  
##res summarizes the genes that have passed the criteria

##get some statistics.
cat( "Total number of upregulated genes ", sum(res160$log2FoldChange> log(FCTh,2) & res160$padj<padjTh, na.rm=TRUE), "\n")



##get list of upregulated genes 
upDE160=rownames(res160)[which(res160$log2FoldChange> log(FCTh,2) & res160$padj<padjTh)]

##get list of downregulated genes 


#combine both lists 
AllDE160=c(upDE160, downDE160)

#save the list of genes data from the differential expression comparison.
filename="GN160_vsControl.csv"
write.csv(as.matrix(res160), file=filename )

##Add functional annotation 
xx=as.matrix(res160)
filename="GN160_vsControl_annotation.csv"
write.table(cbind(xx,annotation[rownames(xx),]), sep="\t", file=filename,row.names=TRUE, col.names=NA)



##Volcano plot
plot(res160$log2FoldChange, -log(res160$padj,10),xlab="log2(fold change)", ylab="-log10(Adj p-value)", pch=20, main="CN160 vs control")
    abline(v=c(-log(1.5,2), log(1.5,2)))
    abline(h=-log(0.05,10))
    up=which(res160$log2FoldChange> log(1.5,2) & res160$padj<0.05)
    down=which(res160$log2FoldChange < - log(1.5,2) & res160$padj<0.05)
    points(res160$log2FoldChange[up], -log(res160$padj[up],10), col="blue",pch=20)
    points(res160$log2FoldChange[down], -log(res160$padj[down],10), col="red", pch=20)




#compare exposure GN800 vs control (control)
##here "condition" is as in the colData file (with sample information) and then what you want to compare. 
res800=results(dds, contrast=c("condition","GN800_2h", "Control"),alpha=padjTh)  


###get the list of genes printed in the screen

cat(upDE800)



####OPTIONAL 
####Create DESeq2 object with the new design
dds2 <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ compound)  #"compound" is the colname in coldata


#normalize the data
dds2 <- DESeq(dds2)
##compare two conditions 
FCTh=1.5 #FC threshold to be used. 
padjTh=0.05  #FDR threshold to be used. 

##compare exposure  vs control (control)
##here "compound" is as in the colData file (with sample information) and then what you want to compare. 
res=results(dds2, contrast=c("compound","GN", "none"),alpha=padjTh)  
##res summarizes the genes that have passed the criteria

##get some statistics.
cat( "Total number of upregulated genes ", sum(res$log2FoldChange> log(FCTh,2) & res$padj<padjTh, na.rm=TRUE), "\n")
cat( "Total number of downregulated genes ", sum(res$log2FoldChange < -1* log(FCTh,2) & res$padj<padjTh, na.rm=TRUE), "\n")


