# This code receives normalized Seurat annotated object: with "severity", "infection" and "celltype" annotations.
# The code compares COVID-19 patients and healthy and all severity levels.
# output: cell type-specific graphs and differential expression tables.

library(Seurat)
library(reshape2)
library(ggplot2)
library(openxlsx)
library(RColorBrewer)
library(gridExtra)
library(stringr)
library(tidyr)
library(ggpubr)


### FUNCTIONS ###

filtermitozero<-function(ob,n,genes=mito.genes){
  # exclude cells with more than n zeros in mtdna genes/ other genes.
  
  df<- as.data.frame(ob@assays$RNA@data[genes,])
  exclude = names(df)[colSums(pmax(df == 0)) > n]
  ob = subset(ob, cells = colnames(ob)[!colnames(ob) %in% exclude])
  return(ob)
}


subset_obj <- function(ob,n=0,cells){
  # This function filters cells with more than n zeros in mtDNA genes and subset the object accordingly.
  
  ob = filtermitozero(ob,n)
  return(subset(ob,cells = rownames(ob@meta.data[ob@meta.data$celltype %in% cells,])))
}


do_plot<-function(table,metadata,myname,k,toadd,maingraphslist,celltype){
  # This function creates a graph comparing severity levels. 
  # output: a list with graphs.
  
  mylistofgenes = c(mito.genes)
  table = table[colnames(table) %in% c("severity",mylistofgenes)]
  ta = melt(table)
  colnames(ta) = c("severity","gene","value")
  ta$gene = ordered(ta$gene, levels = mylistofgenes)
  ta=  ta[ta$gene %in% mito.genes,] 
  maingraphslist[[myname]] = ggplot(ta, aes(x=severity, y=value, fill=severity)) +
    geom_boxplot()+#geom_violin()+ #geom_boxplot()
    facet_wrap(~gene, nrow=1, ncol=length(rownames(mito.genes)), strip.position = "bottom")+
    #geom_jitter( position=position_jitter(0.2))+ #shape=16,
    theme_bw() + ylab("Normalized Expression")+xlab("mtDNA genes") + 
    ggtitle(paste(celltype,", N = ",length(rownames(metadata))))+
    theme(text = element_text(size=20),axis.text.x = element_text(size=20,angle = 60,hjust = 1))+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),axis.text.y = element_text(size=25),
          strip.text = element_text(size=25))#hjust=1 #
  return(maingraphslist)
}


getgenes<- function(DE, mygenes){
  # This function subset a set of genes from the original differential expression list
  # correct it for multiple testing
  de = DE[rownames(DE) %in% mygenes,]
  if (length(rownames(de))<1){return(de)}
  de$p_val_adj_FDR = p.adjust(de$p_val, method = "fdr") 
  de$p_val_adj_BH = p.adjust(de$p_val, method = "BH") 
  return(de)  
}


dode <- function(obtmp,celltype,k,ident="severity",nametoadd,namestotest,listtotest){
  # This function writes the excel file of differential expressed genes
  name = celltype
  Idents(obtmp) = ident
  d = FindAllMarkers(obtmp,only.pos=T,min.pct=0.25, logfc.threshold=0,return.thresh=1)
  ptable = getgenes(d,listtotest[[1]])
  start=2
  while(nrow(ptable)<1 && start<length(namestotest)+1){start=start+1;ptable = getgenes(d,listtotest[[start]])}
  if(start<length(namestotest)+1){
    ptable <- ptable %>% drop_na()
    ptable[is.na(ptable)] = 0
    ptable$pathway = namestotest[1]
    ptable$cell_type = celltype
    ptable$dataset = dataset_name
    for (i in start:length(namestotest)){
      pathway = namestotest[i]
      mygenes = listtotest[[i]]
      ptabletmp = getgenes(d,listtotest[[i]])
      if (length(rownames(ptabletmp))>0){
        ptabletmp <- ptabletmp %>% drop_na()
        ptabletmp[is.na(ptabletmp)] = 0
        ptabletmp$pathway = namestotest[i]
        ptabletmp$cell_type = celltype
        ptabletmp$dataset = dataset_name
        ptable = rbind(ptable,ptabletmp)
      }
    }
    
    list_of_datasets <- list("All_DE_genes" = d,"mtDNA"=getgenes(d,mito.genes), "pathways" = ptable)
    write.xlsx(list_of_datasets, file = paste(nametoadd,celltype,ident,k,"zero mtDNA DE genes.xlsx"))
  }
}


findallmarkerssample <- function(obtmp,n,n.cells,findmarkerslist,celltype){
  # This function compares healthy to covid-19 patients by bootstrap.
  # Each iteration applies FindAllMarkers function
  
  Idents(obtmp) = "infection"
  obtmph<- subset(obtmp,idents="healthy")
  obtmpsars2<- subset(obtmp,idents="SARS-CoV-2")
  df<-data.frame(gene=mito.genes)
  dflogfc<- data.frame(gene= mito.genes)
  dfgroup<-data.frame(gene= mito.genes)
  healthy = 0
  covid19 = 0
  for (i in 1:n){
    set.seed(i)
    iter=toString(i)
    cells = c(sample(x = colnames(obtmph), size = n.cells),sample(x = colnames(obtmpsars2), size = n.cells))
    obtmpsamp<-subset(obtmp,cells=cells)
    markers<-FindAllMarkers(obtmpsamp, only.pos=T,min.pct=0,return.thresh=1,logfc.threshold =0,features=mito.genes,cells = colnames(cells))
    markers <- markers[match(intersect(mito.genes,rownames(markers)), rownames(markers)),]
    getg<-getgenes(markers,mito.genes) # returns adjusted p-value
    toadd<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = getg[getg$gene==x,]$p_val_adj_BH}else{x=1}}))
    toaddlogfc<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = getg[getg$gene==x,]$avg_logFC}else{x=0}}))
    toaddgroup<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = as.character(getg[getg$gene==x,]$cluster)}else{x="-"}}))
    df[,paste0("x",i)] = toadd
    dflogfc[,paste0("x",i)] = toaddlogfc
    dfgroup[,paste0("x",i)] = toaddgroup 
  }	
  iterdf = df[,2:ncol(df)]
  nosignum <- rowSums(iterdf > 0.05)
  iterdf$avg_logFC <- rowMeans(dflogfc[,-1])
  iterdf$pval <- as.numeric(nosignum)/n
  iterdf$group_healthy <- rowSums(dfgroup == "healthy")
  iterdf$group_covid19 <- rowSums(dfgroup == "SARS-CoV-2")
  
  print(head(iterdf$pval))
  findmarkerslist[[celltype]] = data.frame(gene=mito.genes,adj_pvalue=iterdf$pval, avg_logFC = iterdf$avg_logFC,
  group_healthy = iterdf$group_healthy, group_covid19 = iterdf$group_covid19, cell_type = celltype, dataset = dataset_name)
  return(findmarkerslist)
}




findallmarkerssample_severity <- function(obtmp,n,n.cells,findmarkerslist,celltype){
  # This function compares severity levels by bootstrap.
  # Each iteration applies FindAllMarkers function
  
  Idents(obtmp) = "severity"
  comparisons = list(c("control","moderate"),c("control","severe"),c("moderate","severe"))
  dffinal = data.frame(gene=mito.genes)
  for (c in comparisons){
    obtmph<- subset(obtmp,idents=c[1])
    obtmpsars2<- subset(obtmp,idents=c[2])
    df<-data.frame(gene=mito.genes)
    dflogfc<- data.frame(gene= mito.genes)
    dfgroup<-data.frame(gene= mito.genes)
    healthy = 0
    covid19 = 0
    for (i in 1:n){
      set.seed(i)
      iter=toString(i)
      cells = c(sample(x = colnames(obtmph), size = n.cells),sample(x = colnames(obtmpsars2), size = n.cells))
      obtmpsamp<-subset(obtmp,cells=cells)
      markers<-FindAllMarkers(obtmpsamp, only.pos=T,min.pct=0,return.thresh=1,logfc.threshold =0,features=mito.genes,cells = colnames(cells))
      markers = markers[match(intersect(mito.genes,rownames(markers)), rownames(markers)),]
      getg<-getgenes(markers,mito.genes) # returns adjusted p-value
      toadd<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = getg[getg$gene==x,]$p_val_adj_BH}else{x=1}}))
      toaddlogfc<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = getg[getg$gene==x,]$avg_logFC}else{x=0}}))
      toaddgroup<- unlist(lapply(mito.genes,function(x){if(x %in% rownames(markers)){x = as.character(getg[getg$gene==x,]$cluster)}else{x="-"}}))
      df[,paste0("x",i)] = toadd
      dflogfc[,paste0("x",i)] = toaddlogfc
      dfgroup[,paste0("x",i)] = toaddgroup
    }	
    nosignum <- rowSums(df[,2:ncol(df)] > 0.05)
    print("rowsum:")
    print(nosignum)
    avg_logFC <- rowMeans(dflogfc[,-1])
    pval <- as.numeric(nosignum)/n
    group_healthy <- rowSums(dfgroup == c[1])
    group_p <- rowSums(dfgroup == c[2])
    print("haha")
    dffinal[,paste0("adj_pvalue ","( ",c[1]," vs ",c[2],")")]=pval
    dffinal[,paste0("avg_logFC ","( ",c[1]," vs ",c[2],")")]=avg_logFC
    dffinal[,paste0("group ",c[1],"( ",c[1]," vs ",c[2],")")]=group_healthy
    dffinal[,paste0("group ",c[2],"( ",c[1]," vs ",c[2],")")]=group_p
  }
  dffinal$cell_type = celltype
  dffinal$dataset= dataset_name
  findmarkerslist[[celltype]] = dffinal
  return(findmarkerslist)
}


celltype_analysis<-function(rds1,mycellstypes,tocountzerogenes,nametoadd){
  # This function iterate per cell type and calls the functions above to create graphs and DE tables.
  
  ob = rds1
  ob = subset_obj(rds1,i,mycellstypes) 
  mt.table<-as.data.frame(ob@assays$RNA@data[mito.genes,])
  table = cbind(ob@meta.data,t(mt.table))
  mymetadata = ob@meta.data
  
  # initialize lists
  geneprop = list() ; cellsprop = list() ; genesnum = list() ; cellsnum = list() ; findmarkerslist = list(); findmarkerslist2 = list()
  qcgraphlist = list() ; maingraphslist = list() ; maingraphslist2 = list() ; maingraphslist_nuc = list()
  pvaluecompare = list() ; pvaluecompare2 = list() ; myepithel = c() ; mydf = list() ; mydf2 = list()
  
  Idents(ob) = "celltype"
  for (celltype in levels(as.factor(ob@meta.data$celltype))){
  
    cells = WhichCells(ob, idents = celltype)
    obtmp = subset(ob,cells=cells)
    metmp = ob@meta.data[cells,]
    tatmp = as.data.frame(ob@assays$RNA@data[,cells])
    if (length(colnames(tatmp))>300){
      if (length(rownames(metmp[metmp$severity=="severe",]))>50 &&
          length(rownames(metmp[metmp$severity=="moderate",]))>50 &&
          length(rownames(metmp[metmp$severity=="control",]))>50){
        Idents(obtmp) = "severity"
        listofmygenes = pathways 
        toadd = pathways_names
        
        todobox = T
        if (todobox){
          table = cbind(metmp[colnames(tatmp),],t(tatmp[mito.genes,]))
          maingraphslist = do_plot(table,metmp,paste(paste(celltype,"healthy vs covid"),"boxplot"),i,toadd,maingraphslist,celltype)
        }
        
        tofindmitomarkers = T
        if (tofindmitomarkers){
          findmarkerslist=findallmarkerssample(obtmp,1000,50,findmarkerslist,celltype)
          findmarkerslist2=findallmarkerssample_severity(obtmp,1000,50,findmarkerslist2,celltype)
        }
        
        totestde = T
        # Differential expression 
        if (totestde){
          dode(obtmp,celltype,i,ident="severity",nametoadd,pathways_names,pathways) # de of genes' sets
          dode(obtmp,celltype,i,ident="infection",nametoadd,pathways_names,pathways) # de of genes' sets
        }
      }
    }
  }
  
  
  if (length(findmarkerslist)>0){
    write.xlsx(findmarkerslist,file = paste(nametoadd,i," markers_perm.xlsx"))
  }
  if (length(findmarkerslist2)>0){
    write.xlsx(findmarkerslist2,file = paste(nametoadd,i," markers_perm severity.xlsx"))
  }
  
  
  if (length(pvaluecompare2)>0){
    write.xlsx(pvaluecompare2, file = paste(nametoadd,i,"mt-zero-cells p-values permutations no-nuc-filter.xlsx"))
  }
  # write and save graphs and tables:
  if(length(geneprop)>0){
    write.xlsx(geneprop, file = paste(nametoadd,i,"zero cells filter proportions division no nuc filter.xlsx"))
    write.xlsx(cellsprop, file = paste(nametoadd,i,"zero genes filter proportions division no nuc filter.xlsx"))
    write.xlsx(genesnum, file = paste(nametoadd,i,"zero genes filter numbers division no nuc filter.xlsx"))
    write.xlsx(cellsnum, file = paste(nametoadd,i,"zero cells filter numbers division no nuc filter.xlsx"))
  }
  
  printplots<-function(mygraphslist,numr,numc){
    print(head(mygraphslist))
    m = 1:length(mygraphslist)
  }
  
  if (length(maingraphslist)>0){
    pdf(paste(nametoadd,"graphs of",i,"zero mtGenes mean expression celltype division MITO no nuc filter.pdf"), height=6, width=20)
    print(marrangeGrob(maingraphslist, nrow=1, ncol=1,top=NULL))
    dev.off()
  }
  saveRDS(maingraphslist,paste(nametoadd,"maingraphslist.rds"))
  
  if (length(qcgraphlist)>0){
    pdf(paste(nametoadd,"graphs of",i,"zero mtGenes mean expression celltype quality control division no nuc filter.pdf"), height=16, width=24)
    print(marrangeGrob(qcgraphlist, nrow=4, ncol=4,top=NULL))   
    dev.off()
  }
  
}


##############  MAIN  ##################

# upload seurat annotated rds file
ob = readRDS("obj_with_anno.rds")
dataset_name = "my_dataset" # change according to the dataset name

# upload genes:
ad<- "genes/"
pathways_names =unlist(lapply(list.files(ad,full.names = F),function(x){x<-strsplit(x,"[.]")[[1]][1]}))
pathways = lapply(list.files(ad,full.names = T),function(x){x<-as.character(read.table(x,header=T)$gene)})
mito.genes =  c("MT-ND1",  "MT-ND2",  "MT-ND3",  "MT-ND4",  "MT-CO1",  "MT-CO2",  "MT-CO3",  "MT-CYB",  "MT-ATP6")

rds1 = ob
i=0
celltype_analysis(rds1,mycellstypes = levels(as.factor(rds1@meta.data$celltype)),tocountzerogenes=F,nametoadd="") 

########################################







