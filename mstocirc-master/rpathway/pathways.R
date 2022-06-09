#!/usr/bin/Rscript



# refer to link,https://www.jianshu.com/p/feaefcbdf986

# which Rscipt
# chmod a+x ???.R

# install clusterProfiler for R
if(require("clusterProfiler")){
  print ('load cluster succ')
} else {
  print ('tempt to install...')
  #source('http://bioconductor.org/biocLite.R')
  #biocLite('clusterProfiler')
  if(!requireNamespace('BiocManager',quietly=TRUE))
     install.packages('BiocManager')
  BiocManager::install('clusterProfiler')
  if(require("clusterProfiler")){
    print('succ install and load')
  } else {
    stop('install failed.')
  }
}


Args <- commandArgs()

spa <- c('hsa','ath')
spb <- c('org.Hs.eg.db','org.At.tair.db')

dd <- which(spa==Args[7])
#dd <- 2
spb[dd]
hx <- read.table(paste(c(Args[6],'query3.tsv'),sep="",collapse="/"),sep='\t',header=T)
#hx <- read.table('/home/qumu/query3.tsv',sep='\t',header=T)
hx
x <- hx$host_gene
x

library(BiocManager)
#BiocManager::install("org.Hs.eg.db")
if(!requireNamespace(spb[dd],quietly=TRUE)){
BiocManager::install(spb[dd])}

library(clusterProfiler)
#library(org.Hs.eg.db)
#library(spb[2])
#if(spb[dd]=='org.At.tair.db'){library(org.At.tair.db)}
#library(org.At.tair.db)
switch(Args[7],'ath'= library(org.At.tair.db),'hsa'=library(org.Hs.eg.db))

#fromtype,ncbi-, TAIR, SYMBOL, 
#egeg <- bitr(x,fromType='SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
if(Args[7]!='ath'){
bitr_gene <- bitr(x,fromType='SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb=spb[dd])}else{
bitr_gene <- bitr(x,fromType='TAIR',toType=c('SYMBOL','ENTREZID','ENZYME'),OrgDb=org.At.tair.db)}
#ARACYC, ARACYCENZYME, ENTREZID, ENZYME, EVIDENCE, EVIDENCEALL, GENENAME, GO, GOALL,
#ONTOLOGY, ONTOLOGYALL, PATH, PMID, REFSEQ, SYMBOL, TAIR.
bitr_gene
#degenes <- read.csv(Args[6],header=T,stringsAsFactors=F,sep='\t')
#genelist <- degenes$Enterez.ID
genelist <- bitr_gene$ENTREZID
genelist[duplicated(genelist)]

genelist
#GO term

#go <- enrichGO(genelist,OrgDb = org.Hs.eg.db, ont = 'ALL', pAdjustMethod = 'BH',pvalueCutoff = 0.5,qvalueCutoff = 0.2,keyType='ENTREZID')
#go <- enrichGO(genelist,OrgDb =org.At.tair.db, ont = 'ALL', pAdjustMethod = 'BH',pvalueCutoff = 0.5,qvalueCutoff = 0.2,keyType='ENTREZID')
go <- enrichGO(genelist,OrgDb =spb[dd], ont = 'ALL', pAdjustMethod = 'BH',pvalueCutoff = 0.5,qvalueCutoff = 0.2,keyType='ENTREZID')
head(go)
go
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])

jpeg(file=paste(c(Args[6],'/1.jpg'),sep="",collapse="/"), width=1200, height=1000, units='px',res=72*2)
#jpeg(file='/home/qumu/1.jpg')
#jpeg(file,width=1200,height=600)

barplot(go, showCateogy=20)
    library(ggplot2)
ggplot(data=go)+geom_bar(aes(x=Description, y=-log10(pvalue),fill=type),stat='identity')+coord_flip()+scale_x_discrete(limits=go$Description)
dev.off()
jpeg(file=paste(c(Args[6],'/2.jpg'),sep="",collapse="/"), width=1200, height=1000, units='px',res=72*2)
dotplot(go,showCategory=20)
dev.off()
jpeg(file=paste(c(Args[6],'/3.jpg'),sep="",collapse="/"), width=2400, height=2800, units='px',res=72*4)
#jpeg(file='/home/qumu/3.jpg')
#go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff = 0.2,keyType='ENTREZID')
go.BP <- enrichGO(genelist,OrgDb=spb[dd],ont='BP',pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff = 0.2,keyType='ENTREZID')
plotGOgraph(go.BP)
dev.off()

#KEGG
#kegg <- enrichKEGG(genelist,organism=Args[7],keyType='kegg',pvalueCutoff=0.05,pAdjustMethod='BH',minGSSize=10,maxGSSize=500,qvalueCutoff=0.2,use_internal_data=FALSE)
#kegg <- enrichKEGG(genelist,organism='ath',keyType='kegg',pvalueCutoff=0.05,pAdjustMethod='BH',minGSSize=10,maxGSSize=500,qvalueCutoff=0.2,use_internal_data=FALSE)
if(Args[7]=='ath'){
kegg_gene <- x[duplicated(x)] }else{kegg_gene <- genelist}
kegg <- enrichKEGG(kegg_gene,organism=Args[7],keyType='kegg',pvalueCutoff=0.05,pAdjustMethod='BH',minGSSize=10,maxGSSize=500,qvalueCutoff=0.2,use_internal_data=FALSE)
head(kegg)
dim(kegg)
jpeg(file=paste(c(Args[6],'/4.jpg'),sep="",collapse="/"), width=1200, height=1000, units='px',res=72*2)
#jpeg(file='/home/qumu/4.jpg')
dotplot(kegg,showCategory=20)
dev.off()   
