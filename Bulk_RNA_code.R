# This code receives raw matrix (df) with samples (columns) and genes (rows) and a metadata file (meta)
# output: differential expression tables and graphs.

library(reshape2)
library(ggplot2)
library(openxlsx)
library(RColorBrewer)
library(DESeq2)


dograph <- function(metadata,matr,name,genes){ 
         # input: 
         # metadata - an information metadat file, with sample names in rows.
         # matr - a dataframe of genes in columns, samples in rows
         # name - a name to add to the graph
         # genes - the genes to plot in the graph.
         # --------------------------------------
         # output: writes a graph to pdf file.
         
         df= matr[,intersect(genes,colnames(matr))]
         m = melt(cbind(metadata,df))
         if (!name=="mito"){
           p1 = ggplot(m, aes(x=health, y=value, fill=health)) +
           geom_boxplot()+ggtitle(paste(name))+facet_wrap(~variable,nrow=1, strip.position = "bottom",scale="free_y")+
           theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+xlab("Gene")+ylab("Normalized Expression")+
           geom_jitter(color="black", size=0.4, alpha=0.9) +
           scale_fill_manual(values=c("#F8766D", "#00BA38", "#00BFC4"))+
           theme_bw() + ylab("Normalized Expression") + ggtitle(paste("Healthy vs COVID N = ",length(rownames(metadata))))+
           theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20,angle = 60,hjust = 1))+
           theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=25),
           strip.text = element_text(size=25),
           axis.ticks.x = element_blank())#+ylim(0,3500)
         }else{
           p1 = ggplot(m, aes(x=health, y=value, fill=health)) +
           geom_boxplot()+ggtitle(paste(name))+facet_wrap(~variable,nrow=1, strip.position = "bottom")+
           theme(axis.text.x = element_text(angle = 60, vjust = 0.5, hjust=1))+xlab("mtDNA Gene")+ylab("Normalized Expression")+
           geom_jitter(color="black", size=0.4, alpha=0.9) +
           scale_fill_manual(values=c("#F8766D", "#00BA38", "#00BFC4"))+
           theme_bw() + ylab("Normalized Expression") + ggtitle(paste("Healthy vs COVID N = ",length(rownames(metadata))))+
           theme(text = element_text(size=20),axis.text.y = element_text(size=20),axis.text.x = element_text(size=20,angle = 60,hjust = 1))+
           theme(axis.text.x = element_blank(),
           axis.text.y = element_text(size=25),
           strip.text = element_text(size=20),
           axis.ticks.x = element_blank())#+ylim(0,75000)
         }
        if (name=="mito"){w=20}else{w=10}
        pdf(paste(name,".pdf"),width=w,height = 6)
        print(p1)
        dev.off()
}


getgenes<- function(DE, mygenes){
    # This function subset a set of genes from the original differential expression list
    # correct for multiple testing.
    # ----------------------------------------------------------------------------------
    # input:
    # DE - data frame of results(dds) deseq object.
    # mygenes - genes to test.
    DE <- DE[!is.na(DE$padj),]
	  de = DE[rownames(DE) %in% intersect(mygenes,rownames(DE)),]
    if (length(rownames(de))<1){return(de)}
	  de$p_val_adj_FDR = p.adjust(de$pvalue, method = "fdr") 
    de$p_val_adj_BH = p.adjust(de$pvalue, method = "BH") 
	  de = de[match(mygenes, de$gene),]
    return(de)
}


dode <- function(d,name,namestotest,listtotest){
 # This function writes the excel file differential expressed genes
 # ----------------------------------------------------------------
 # input:
 # d - data frame of results(dds) deseq object.
 # name - "COVID-19" patients
 # namestotest - pathway names
 # listtotest - a list of pathways' genes, respectively.
 
 ptable = getgenes(d,listtotest[[1]])
 ptable <- ptable %>% drop_na()
 ptable[is.na(ptable)] = 0
 ptable$group = unlist(lapply(ptable$log2FoldChange, function(x){if(x<0){x<-paste("healthy (",name," vs healthy)")}else{x<-name}}))
 ptable$pathway = namestotest[1]
 ptable$dataset = dataset_name
 for (i in 2:length(namestotest)){
     pathway = namestotest[i]
     mygenes = listtotest[[i]]
     ptabletmp = getgenes(d,listtotest[[i]])
     ptabletmp <- ptabletmp %>% drop_na()
     ptabletmp[is.na(ptabletmp)] = 0
     ptabletmp$pathway = namestotest[i]
     ptabletmp$dataset = dataset_name
     ptabletmp$group = unlist(lapply(as.numeric(ptabletmp$log2FoldChange), function(x){if(x<0){x<-paste("healthy (",name," vs healthy)")}else{x<-name}}))
     ptable = rbind(ptable,ptabletmp)
 }

 list_of_datasets <- list("All_DE_genes" = d,"mtDNA"=getgenes(d,mito.genes), "pathways" = ptable)
 write.xlsx(list_of_datasets, file = paste(name,"DE_groups.xlsx"))
}


# upload genes to test
mito.genes = c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ATP6", "MT-ATP8")  
pathways_names =unlist(lapply(list.files(ad,full.names = F),function(x){x<-strsplit(x,"[.]")[[1]][1]}))
pathways = lapply(list.files(ad,full.names = T),function(x){x<-as.character(read.table(x,header=T)$gene)})

# upload files
df<- read.table("raw_df.txt",header=TRUE)
info = read.csv("metadata.csv",header = T)
dataset_name = "" # fill in the name of the dataset
health =  unlist(lapply(as.character(colnames(df)),function(x){x<-as.character(info[as.character(info$sample) %in% x,]$health)}))
batch = unlist(lapply(as.character(colnames(df)),function(x){x<-as.character(info[as.character(info$sample) == x,]$batch)}))
metadata = as.data.frame(health)
metadata$batch = batch

if (!file.exists("dds.rds")){
  df <- df[index,]
  dds <- DESeqDataSetFromMatrix(countData = round(df), colData = metadata, design = ~ batch + health)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds, betaPrior = FALSE)
  res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
  saveRDS(dds, "dds.rds")
}else{dds = readRDS("dds.rds")}

res <- as.data.frame(results(dds))
res$gene = rownames(res)
write.csv(res,"degenes.csv")

# extract normalized matrix:
normdf = counts(dds, normalized=T)

# subset mito genes and plot graphs:
normdf = as.data.frame(t(normdf))
dograph(metadata,normdf,"mito",mito.genes) #bonferroni
dograph(metadata,normdf,"nuc",c("GAPDH","LDHA","LDHB"))
dograph(metadata,normdf,"nuc_reg",c("JUND","JUN","POLRMT"))

# test for differential expressed genes:
dode(res,"COVID-19",pathways_names,pathways)

