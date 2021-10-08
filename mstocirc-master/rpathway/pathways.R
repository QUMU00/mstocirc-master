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

hx <- read.table(paste(c(Args[6],'query3.tsv'),sep="",collapse="/"),sep='\t',header=T)
#x <- read.table('/home/qumu/query3.tsv',sep=',')
hx[,5]
x <- hx$host_gene

BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
egeg <- bitr(x,fromType='SYMBOL',toType=c('ENTREZID','ENSEMBL'),OrgDb='org.Hs.eg.db')
egeg
#degenes <- read.csv(Args[6],header=T,stringsAsFactors=F,sep='\t')
#genelist <- degenes$Enterez.ID
genelist <- egeg$ENTREZID
genelist
genelist[duplicated(genelist)]

genelist
#GO term

go <- enrichGO(genelist,OrgDb = org.Hs.eg.db, ont = 'ALL', pAdjustMethod = 'BH',pvalueCutoff = 0.5,qvalueCutoff = 0.2,keyType='ENTREZID')
head(go)
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])

jpeg(file=paste(c(Args[6],'/1.jpg'),sep="",collapse="/"))
#jpeg(file='/home/qumu/1.jpg')
#jpeg(file,width=1200,height=600)
barplot(go,showCateogy=20,drop=T)
dev.off()
jpeg(file=paste(c(Args[6],'/2.jpg'),sep="",collapse="/"))
dotplot(go,showCategory=50)
dev.off()
jpeg(file=paste(c(Args[6],'/3.jpg'),sep="",collapse="/"))
#jpeg(file='/home/qumu/3.jpg')
go.BP <- enrichGO(genelist,OrgDb=org.Hs.eg.db,ont='BP',pAdjustMethod='BH',pvalueCutoff=0.05,qvalueCutoff = 0.2,keyType='ENTREZID')
plotGOgraph(go.BP)
dev.off()

#KEGG
kegg <- enrichKEGG(genelist,organism='hsa',keyType='kegg',pvalueCutoff=0.05,pAdjustMethod='BH',minGSSize=10,maxGSSize=500,qvalueCutoff=0.2,use_internal_data=FALSE)
head(kegg)
dim(kegg)
jpeg(file=paste(c(Args[6],'/4.jpg'),sep="",collapse="/"))
#jpeg(file='/home/qumu/4.jpg')
dotplot(kegg,showCategory=30)
dev.off()














