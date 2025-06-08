#--------------USA cohort-----------------

library(tidyverse)
library(edgeR)
library(limma)
library(ComplexHeatmap)
#library(annotables)
library(org.Hs.eg.db)
library(DESeq2)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(GOSemSim)


getwd()
#root="/mmfs1/home/david.adeleke/WorkFlow"

#root <- "~"
#main <- "Chapter_Two"
subf="validation"

# Set root and main directories
root <- "/mmfs1/home/david.adeleke/Thesis"
main <- "RNASeq"
proj.scope <- c("PAAD-US", "PACA-CA", "PACA-AU")
# Create directory if it doesn't exist and set working directory
uscontrol= c("DO32751", "DO218973", "DO50272", "DO32769")

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main,  recursive = TRUE)))}

setwd(file.path(root, main))

normseq <- readRDS(file.path(root, "dataset/PAAD-US_ICGC_rnaseq.rds"))

normseq $class=substring(normseq$submitted_sample_id, 14,16)

normseq$icgc_donor_id =ifelse(normseq$class  == "11A", paste0(normseq$icgc_donor_id,"N"), normseq$icgc_donor_id)

rnaseq=normseq[normseq$class !="11A",]

rnaseq=rnaseq[,-23]

normseq=normseq[normseq$class=="11A",]

normseq=normseq[,-23]

rnaseq=rbind(rnaseq, normseq)

rnaseq=rnaseq %>% dplyr::select(c(1,8,10)) %>% filter(!gene_id %in% c("?", "SLC35E2"))

rnaseq= tidyr::pivot_wider(rnaseq, id_cols= "gene_id", names_from = "icgc_donor_id" , values_from = "raw_read_count", values_fn= mean)

rnaseq=as.data.frame(rnaseq)

rownames(rnaseq)= rnaseq$gene_id

rnaseq=rnaseq[,-1]

rnaseq[is.na(rnaseq)] <- 0


#load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/expurity.RData")
#purit=puree[,]
load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

usidtable$grp= ifelse(usidtable$grp == "KnTnHB","KnTn", usidtable$grp)

commonid= intersect(usidtable$icgc_donor_id, colnames(rnaseq))

usidtable= usidtable[usidtable$icgc_donor_id %in%  commonid, ]

rownames(usidtable)= usidtable$icgc_donor_id

rnaseq=rnaseq[,commonid]

dim(rnaseq)
usidtable$grp=as.factor(usidtable$grp)


dds <- DESeqDataSetFromMatrix(countData=rnaseq, colData=usidtable, design = ~ grp )


#srt= 0.05*136
dds75 <- dds[rowSums(counts(dds) >= 5) >= 3,]
nrow(dds75) # 

nrow(dds)
# 

dds <- DESeq(dds75)
vsd <- vst(dds, blind=FALSE)
norm.counts <- assay(vsd)

normalized_counts <- counts(dds, normalized=TRUE)

#compute tumor purity from outside

normid = colnames(norm.counts)[substring(colnames(norm.counts),nchar(colnames(norm.counts)),nchar(colnames(norm.counts)))=="N"]

#rnaseqt=norm.counts[,setdiff(colnames(norm.counts), normid)]

#rnaseqt=t(rnaseqt)
#bkcol= colnames(rnaseqt)
#rnaseqt=as.data.frame(rnaseqt)
#write_csv(rnaseqt, file = "usrnaseq2.csv", col_names=T)
#saveRDS(bkcol, "uspuritycolnms.tsv")


#compute purity using puree

puree=read_tsv(file.path(root,"dataset/us_rnaseq_purity.txt"))

puree=as.data.frame(puree)

mypcoln=readRDS(file.path(root,"dataset/uspuritycolnms.tsv"))

rownames(puree)= mypcoln

puree$id=rownames(puree)

puree=as.data.frame(puree)

puree=puree[,c(3,2)]


normet=norm.counts[,normid]

tum=norm.counts[,setdiff(colnames(norm.counts),normid)]

sub.id=intersect(colnames(tum),rownames(puree) )
tum=tum[,sub.id]
head(tum )
#purity=purie;rm(purit)
#purity$hhh=1


purity=puree[sub.id,]

Tumor.mat = tum
Norm.mat = normet

Norm.matt = rowMeans(Norm.mat)

lamda = as.numeric(purity[colnames(Tumor.mat), 2])
names(lamda)=sub.id

# Initialize the result matrix J
J = matrix(0, nrow(Tumor.mat), ncol(Tumor.mat))
rownames(J) = rownames(Tumor.mat)
colnames(J) = colnames(Tumor.mat)

# Compute J

for (c in 1:ncol(J)) {
  for (r in 1:nrow(J)) {
    normexp.i=Norm.matt[r]
    puri.i=(1-lamda[c])
    nomexp=puri.i*normexp.i
    admx=Tumor.mat[r, c]
    ym= admx-nomexp
    realexp=ym/lamda[c]
    
    J[r, c] = realexp
  }
}


adj.counts=cbind(J,Norm.mat)

sub.id=intersect(colnames(adj.counts), usidtable$icgc_donor_id)

group <-  as.data.frame(usidtable[usidtable$icgc_donor_id %in% sub.id, c("icgc_donor_id","grp")])

rownames(group)=group$icgc_donor_id

group=group[sub.id,2]
group=as.factor( group)

names(group)= sub.id

adj.counts=adj.counts[,sub.id]

#dge=adj.counts
design <- model.matrix(~0 + group)

colnames(design) <- c(levels(group))

rownames(design)=sub.id


# Model fitting and DE calculation 
# DE genes


contr<-  makeContrasts(KnTp-KnTn,
                       KpTn-KnTn,
                       KpTp-KnTn,
                       KpTn - KnTp,
                       KpTp - KnTp,
                       KpTp - KpTn,
                       KnTn -  WT,
                       KnTp - WT,
                       KpTn - WT,
                       KpTp - WT,
                       levels=design)


fit <- lmFit(adj.counts, design)
fit <- contrasts.fit(fit, contr[,1:10])
fit2 <- eBayes(fit)

dt <- decideTests(fit2, lfc=1)
summary(dt)


lap=0

DE =data.frame()
for(i in 1:10){
  s <- topTable(fit2, num=Inf, coef=i, sort.by="p")
  s$gene=rownames(s)
  s$contrast= colnames(contr)[i+lap]
  DE=rbind(DE,s)
}

DMGs=DE[DE$adj.P.Val <=0.05,]
rownames(DMGs)=1:nrow(DMGs)

DMGs$rank= DMGs$adj.P.Val*DMGs$logFC

#save( adj.counts, file="adj.counts022724.RData")

#save( DE, file="DMGsraw022724.RData")

#de.list=list(fit=fit2, DE=DE, voomObj=adj.counts )

#save( de.list, file="de.list115.RData")

DMGs = DMGs[DMGs$contrast %in% colnames(contr)[c(7:10)], ]

genefeatures=unique(DMGs$gene)


#save( adj.counts, file="adj.countsNewIdea.RData")
#save( genefeatures, file="FeatureNewIdea.RData")



dtt= as.data.frame(dt@.Data)

top_gene=colnames(dtt)

mset=list()
dtt$id=rownames(dtt)
snv.set=dtt

for(i in 7:10){
  vv=colnames(snv.set)[i]
  
  lst=list(unique(as.vector(snv.set[snv.set[[i]]==1,length(top_gene)+1])))
  
  lst2=list(unique(as.vector(snv.set[snv.set[[i]]==-1,length(top_gene)+1])))
  
  names(lst)=paste0("UpReg:",vv)
  names(lst2)=paste0("DownReg:",vv)
  mset <- c(mset,  lst)
  mset <- c(mset,  lst2)
}

remove_empty_lists <- function(x) {
  x_filtered <- lapply(x, function(sublist) {
    if (length(sublist) >= 2) {
      sublist
    } else {
      NULL
    }
  })
  x_filtered <- x_filtered[sapply(x_filtered, Negate(is.null))]
  return(x_filtered)
}



# Apply the function to remove empty nested lists
mset <- remove_empty_lists(mset)
reg= ComplexHeatmap::make_comb_mat(mset, mode="distinct")

#,"#00B0F6", "#A3A500"

setordr= c("UpReg:KnTn - WT",  "UpReg:KnTp - WT", "UpReg:KpTn - WT",  "UpReg:KpTp - WT", 
 "DownReg:KnTn - WT", "DownReg:KnTp - WT", "DownReg:KpTn - WT", "DownReg:KpTp - WT")

setordrb = c("UpReg:KpTn - KnTn",  "UpReg:KpTp - KnTn",  "UpReg:KpTn - KnTp",
 "DownReg:KpTn - KnTn","DownReg:KpTp - KnTn","DownReg:KpTn - KnTp")

colcode= c( "#F8766D","#00BF7D",  "#00B0F6", "#E76BF3")
UpSet(reg, set_order = setordr, comb_col = colcode[comb_degree(reg)], bg_col = c("white", "#E6E6E6"), top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
      right_annotation = upset_right_annotation(reg, add_numbers = TRUE, 
                axis_param = list(side = "top"),
                annotation_name_side = "top", 
                gp = gpar(fill = "#A3A500" )))



#colcode[comb_degree(reg)]

#Glimma::glimmaVolcano(tmp, dge = dge)

############################

gset=list()
for(i in 1:length(top_gene)){
  vv=colnames(snv.set)[i]
  lst=list(unique(as.vector(snv.set[snv.set[[i]] !=0,length(top_gene)+1])))
  names(lst)=vv
  gset <- c(gset,  lst)
}

#"KnTn - WT" "KnTp - WT" "KpTn - WT" "KpTp - WT"
# Name of the list you want to extract

target_name <- colnames(snv.set)[7]

# Extract the list with the specified name
target_list <- gset [[target_name]]

# Find elements unique to the target list
unique_elements <- unlist(target_list) # Convert to a vector
other_elements <- unlist(gset[setdiff(names(gset), target_name)]) # Combine other lists and convert to a vector
unique_to_target <- setdiff(unique_elements, other_elements)

unique_to_target=target_list
# Print the elements unique to "UpReg__A"
#print(unique_to_target)
#--------------------GSEA---------------------------------
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

direct = 0
#DMGs$gene %in% unique_to_target &
gsedf= DMGs[ DMGs$contrast == target_name & DMGs$rank > direct, c(1,7,8) ]


gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")


gsedf=na.omit(gsedf)

gsedf=  gsedf[order(gsedf$logFC), ]


de <- gsedf$logFC

names(de)=gsedf$id


nde=as.character(names(de))

#edo <- enrichDGN(nde,minGSSize = 5)

#dotplot(edo, showCategory=10, orderBy = "x") + ggtitle(paste0("dotplot for ORA in ",target_name) )


# x = enrichDO(nde, pAdjustMethod="none", ont="DO",minGSSize = 10)
# y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
# 
# 
# 
# z=ggplot(y, showCategory = 10, 
#          aes(richFactor, fct_reorder(Description, richFactor))) + 
#   geom_segment(aes(xend=0, yend = Description)) +
#   geom_point(aes(color=p.adjust, size = Count)) +
#   scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
#   scale_size_continuous(range=c(2, 10)) +
#   theme_minimal() + 
#   xlab("enrich factor") +
#   ylab(NULL) +  
#   ggtitle(paste0("Top 10 Enriched Disease Ontology in ", target_name))
# 
# z


#--------------------------------------------------------

gsedf= DMGs[DMGs$contrast == target_name , c(1,7,8,9) ]

gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")

gsedf=na.omit(gsedf)


gsedf=  gsedf[order(gsedf$logFC, decreasing = TRUE), ]


#gsedf=gsedf[gsedf$rank < 0, ]

de <- gsedf$logFC

names(de)=gsedf$id


geneList=de

#, scoreType = "pos"

x <- gsePathway(geneList,minGSSize = 10, pvalueCutoff = 0.1)

pathdf=as.data.frame(x)


#viewPathway("Interferon Signaling", 
            #readable = TRUE, foldChange = geneList)



y <- arrange(x, abs(NES)) %>% 
  group_by(sign(NES)) %>% 
  slice(1:5)


ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL) +
  ggtitle(paste0("Enriched Pathways in ", target_name))



###------------------------------------New VIZ---------------------------

pathdf=data.frame()

for(i in c(1:6)){
  
  target_name <- colnames(snv.set)[i]
  
gsedf= DMGs[DMGs$contrast == target_name , c(1,7,8,9) ]


if(nrow(gsedf)==0){
  print(paste0( target_name, " has no differentially expressed gene"))
}else{
  
gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")

gsedf=na.omit(gsedf)

gsedf=  gsedf[order(gsedf$logFC, decreasing = TRUE), ]

de <- gsedf$logFC

names(de)=gsedf$id

geneList=de

x <- gsePathway(geneList,minGSSize = 5, pvalueCutoff = 0.05)
}
if(nrow(x@result)==0){
  print(paste0( target_name, " has no enriched term"))
}else{
  
path=as.data.frame(x)

path$contrst= target_name

pathdf=rbind(pathdf,path )

}
}

pathdf=pathdf[pathdf$setSize >1, ]


#save( pathdf, file="pathdf.RData")


pathdfwt= pathdf[pathdf$contrst %in% colnames(snv.set)[1:6],]

#----------------------------------------------------

mset=list()

for(i in unique(pathdfwt$contrst)){

  lst= pathdfwt %>% filter(contrst==i , enrichmentScore > 0) %>% pull(ID) %>% list()
  
  lst2=pathdfwt %>% filter(contrst==i , enrichmentScore < 0) %>% pull(ID) %>% list()
  
  if(length(lst[[1]])==0){
    print("zero enriched pathways")
  }else{
    names(lst)=paste0("Up:",i)
  }

  
  if(length(lst2[[1]])==0){
    print("zero enriched pathways")
  }else{
    names(lst2)=paste0("Down:",i)
  }
 
  mset <- c(mset,  lst)
  mset <- c(mset,  lst2)
}


# Apply the function to remove empty nested lists

reg= ComplexHeatmap::make_comb_mat(mset, mode="distinct")

UpSet(reg,   comb_col = c( "#F8766D", "#A3A500","#00BF7D","#00B0F6", "#E76BF3")[comb_degree(reg)], bg_col = c("white", "#E6E6E6"), top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
      right_annotation = upset_right_annotation(reg, add_numbers = TRUE))

setordr= names(mset)[c(1,3,5,7,2,4,6,8)]
setordr= names(mset)[c(1,3,5,2,4,6)]



UpSet(reg, set_order =setordr, comb_col = colcode[comb_degree(reg)], bg_col = c("white", "#E6E6E6"), 
      top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
   right_annotation = upset_right_annotation(reg, add_numbers = TRUE, 
  axis_param = list(side = "top"),
  annotation_name_side = "top", 
  gp = gpar(fill = "#A3A500" )))


## common to all mbs
# UP
setupw=intersect(mset$`Up:KpTn - WT`, 
        intersect(
          intersect(mset$`Up:KnTp - WT`, 
                mset$`Up:KnTn - WT`), 
          intersect(mset$`Up:KpTp - WT`,  mset$`Up:KnTn - WT`)))
#Down
setdown=intersect(mset$`Down:KpTn - WT`, 
                 intersect(
                   intersect(mset$`Down:KnTp - WT`, 
                             mset$`Down:KnTn - WT`), 
                   intersect(mset$`Down:KpTp - WT`,  mset$`Down:KnTn - WT`)))


## common to KpTp & KpTn
# UP
setupw=setdiff(intersect(
                   mset$`Up:KpTp - WT`, 
                             mset$`Up:KpTn - WT`), 
                   union(mset$`Up:KnTn - WT`,  mset$`Up:KnTp - WT`))
#Down
setdown=setdiff(intersect(
  mset$`Down:KpTp - WT`, 
  mset$`Down:KpTn - WT`), 
  union(mset$`Down:KnTn - WT`,  mset$`Down:KnTp - WT`))



## common to KpTp only

# UP
setupw=setdiff(mset$`Up:KpTp - WT`,
  union(union(mset$`Up:KnTn - WT`,  mset$`Up:KnTp - WT`),  mset$`Up:KpTn - WT`))

#Down
setdown=setdiff(mset$`Down:KpTp - WT`,
                union(union(mset$`Down:KnTn - WT`,  mset$`Down:KnTp - WT`),
                      mset$`Down:KpTn - WT`))



## common to KpTn only

# UP
setupw=setdiff(mset$`Up:KpTn - WT`,
               union(union(mset$`Up:KnTn - WT`,  mset$`Up:KnTp - WT`),  mset$`Up:KpTp - WT`))

#Down
setdown=setdiff(mset$`Down:KpTn - WT`,
                union(union(mset$`Down:KnTn - WT`,  mset$`Down:KnTp - WT`),
                      mset$`Down:KpTp - WT`))


(unique(pathdfwt[pathdfwt$ID %in%  setupw,2]))

View((unique(pathdfwt[pathdfwt$ID %in%  setdown,])))

setup= setdiff(mset$`Up:KpTn - WT`, 
              union(
                union(mset$`Up:KnTp - WT`, 
                          mset$`Up:KnTn - WT`), 
               union(mset$`Up:KpTp - WT`,  mset$`Up:KnTn - WT`)))

setdown= setdiff(mset$`Down:KpTn - WT`, 
               union(
                 union(mset$`Down:KnTp - WT`, 
                       mset$`Down:KnTn - WT`), 
                 union(mset$`Down:KpTp - WT`,  mset$`Down:KnTn - WT`)))



setd= setdiff(setup, setdown)

View(unique(pathdfwt[pathdfwt$ID %in%  setdown,]))

unique(pathdfwt[pathdfwt$ID %in%  setup,])

xxx=xx[xx$p.adjust < 0.05,]



