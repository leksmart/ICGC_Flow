#--------------USA cohort-----------------
#library(annotables)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(tidyverse)
library(edgeR)
library(limma)
library(ComplexHeatmap)
options(digits=3)


getwd()
root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"
subf="validation"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main,subf)))}

setwd(file.path(root, main, subf))


mirna_query <- GDCquery(project = "TCGA-PAAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        #workflow.type = "BCGSC miRNA Profiling",
                        experimental.strategy = "miRNA-Seq")

GDCdownload(mirna_query, method = "api", files.per.chunk = 100,
            directory = "~/WorkFlow/Chapter_Two/tcga/")


miR_df <- GDCprepare(mirna_query, directory = "~/WorkFlow/Chapter_Two/tcga/")


cross_df<- miR_df %>% dplyr::select(contains("cross"))
cross_df=cbind(miR_df$miRNA_ID,cross_df)
colnames(cross_df) <- gsub(".*_","",colnames(cross_df))
rownames(cross_df)=cross_df$ID
cross_df=cross_df[,-1]

# Function to count "Y" occurrences per row
count_n <- function(row) {
  sum(row == "Y") / length(row)
}

# Filter rows where "N" count is less than 70%
goodmirna <- cross_df %>%
  filter(apply(., 1, count_n) < 0.001) %>% rownames()

rpm_df<- miR_df %>% dplyr::select(contains("_million_"))
rpm_df=cbind(miR_df$miRNA_ID,rpm_df)

colnames(rpm_df) <- gsub(".*_","",colnames(rpm_df))

rpm_df=rpm_df[rpm_df$ID %in% goodmirna, ]

rownames(rpm_df)=rpm_df$ID
rpm_df=rpm_df[,-1]

rm(cross_df,miR_df)



load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/expurity.RData")
purit=puree[,]
load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

usidtable$submitted_sample_id= substring(usidtable$submitted_sample_id,1,16)

colnames(rpm_df)= substring(colnames(rpm_df),1,16)


col_rpm_df=usidtable[ usidtable$submitted_sample_id %in% colnames(rpm_df),1:2]

rpm_df= rpm_df[,col_rpm_df$submitted_sample_id]

colnames(rpm_df)= col_rpm_df$icgc_donor_id


commonid= intersect(usidtable$icgc_donor_id, colnames(rpm_df))

usidtable= usidtable[usidtable$icgc_donor_id %in%  commonid, ]

metas="TCGA-HZ-A9TJ-06A"

usidtable=usidtable[usidtable$submitted_sample_id != metas,]

rpm_df=rpm_df[,commonid]
usidtable=as.data.frame(usidtable)
rownames(usidtable)=usidtable$icgc_donor_id
usidtable=usidtable[commonid,]


library(DESeq2)

usidtable$grp=as.factor(usidtable$grp)

rnaseq=rpm_df

rnaseq=round(rnaseq)

dds <- DESeqDataSetFromMatrix(countData=rnaseq, colData=usidtable, design = ~ grp )


srt= 0.50*172
dds75 <- dds[rowSums(counts(dds) >= 5) >= srt,]
nrow(dds75) # 
# 

dds <- DESeq(dds75)
#vsd <- vst(dds)
norm.counts <- assay(dds)

normalized_counts <- counts(dds, normalized=TRUE)
#compute tumor purity from outside

normid = colnames(norm.counts)[substring(colnames(norm.counts),nchar(colnames(norm.counts)),nchar(colnames(norm.counts)))=="N"]

#rnaseqt=norm.counts[,setdiff(colnames(norm.counts), normid)]

#rnaseqt=t(rnaseqt)
#bkcol= rownames(rnaseqt)
#rnaseqt=as.data.frame(rnaseqt)
#write_csv(rnaseqt, file = "usrnaseq2.csv", col_names=T)
#saveRDS(bkcol, "uspuritycolnms.tsv")


#compute purity using puree

puree=read_tsv("us_rnaseq_purity.txt")

puree=as.data.frame(puree)

mypcoln=readRDS("uspuritycolnms.tsv")

rownames(puree)= mypcoln

puree$id=rownames(puree)

puree=as.data.frame(puree)

puree=puree[,c(3,2)]


normet=norm.counts[,normid]

tum=norm.counts[,setdiff(colnames(norm.counts),normid)]

sub.id=intersect(colnames(tum),rownames(puree) )
tum=tum[,sub.id]

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

usidtable$grp=as.character(usidtable$grp)
usidtable$grp = ifelse(usidtable$grp == "KnTnHB","KnTn", usidtable$grp  )

usidtable$grp=as.factor(usidtable$grp)

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

dt <- decideTests(fit2, lfc=0)
summary(dt)


lap=0

DE =data.frame()
for(i in 1:10){
  s <- topTable(fit2, num=Inf, coef=i, sort.by="p")
  s$gene=rownames(s)
  s$contrast= colnames(contr)[i+lap]
  DE=rbind(DE,s)
}

DMGs=DE[DE$adj.P.Val <=0.05 & abs(DE$logFC) >= 0.5 ,]
rownames(DMGs)=1:nrow(DMGs)
DMGs$rank=DMGs$adj.P.Val*DMGs$logFC

DMmi=DMGs

#save( DMmi, file="DMmi0104.RData")
adj.mi.counts=adj.counts
#save( adj.mi.counts, file="adj.mi.counts0104.RData")


dtt= as.data.frame(dt@.Data)

top_gene=colnames(dtt)

mset=list()
dtt$id=rownames(dtt)
snv.set=dtt

for(i in 1:10){
  vv=colnames(snv.set)[i]
  lst=list(unique(as.vector(snv.set[snv.set[[i]]==1,length(top_gene)+1])))
  lst2=list(unique(as.vector(snv.set[snv.set[[i]]==-1,length(top_gene)+1])))
  
  names(lst)=paste0("UpReg_",vv)
  names(lst2)=paste0("DwnReg_",vv)
  mset <- c(mset,  lst)
  mset <- c(mset,  lst2)
}

remove_empty_lists <- function(x) {
  x_filtered <- lapply(x, function(sublist) {
    if (length(sublist) > 2) {
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

UpSet(reg,   comb_col = c( "#F8766D", "#A3A500","#00BF7D","#00B0F6", "#E76BF3")[comb_degree(reg)], bg_col = c("white", "#E6E6E6"), top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
      right_annotation = upset_right_annotation(reg, add_numbers = TRUE))


mirnalist= intersect(mset$`DwnReg_KnTn - WT`, mset$`DwnReg_KpTp - WT`)
mirnalist=intersect(mirnalist, mset$`DwnReg_KpTn - WT`)

### Identification of Target genes from miRNA.##############

library(multiMiR)

check_mirna_target= get_multimir(url=NULL,org="hsa",mirna=mirnalist,target=NULL,
                                 disease.drug=NULL,
                                 table="validated",
                                 predicted.cutoff=NULL, 
                                 predicted.cutoff.type="p",predicted.site="conserved",
                                 summary=FALSE,add.link=FALSE,use.tibble=FALSE,limit=NULL)

check_mirna_target= as.data.frame(check_mirna_target@data)

ert=load("DMGsraw1221.RData")


val.target= check_mirna_target$target_symbol
DMGs=DE
DMGs=DMGs[DMGs$gene %in% val.target, ]

#DMGsWT= DMGs[DMGs$contrast %in% DMGs$contrast[grep("WT$", DMGs$contrast)],]
library(ReactomePA)

DMGsWT= DMGs[DMGs$contrast == "KpTp - WT" ,]

gsedf= DMGsWT[ , c(1,5,7,8) ]

gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")


gsedf=na.omit(gsedf)
gsedf$rank=gsedf$adj.P.Val * gsedf$logFC
gsedf=  gsedf[order(gsedf$rank, decreasing = TRUE), ]


de <- gsedf$rank

names(de)=gsedf$id


geneList=de

 x <- gsePathway(geneList, pvalueCutoff = 0.05, pAdjustMethod = "BH", 
 verbose = FALSE)

 
 
 y <- arrange(x, abs(NES)) %>% 
   group_by(sign(NES)) %>% 
   slice(1:5)
 
 
 ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
   geom_col(orientation='y') + 
   scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
   theme_minimal() + ylab(NULL)
 
 
 y=as.data.frame(x)


 viewPathway("Immune System", 
             readable = TRUE, foldChange = geneList)

 library("pathview")
 
 hsa04110 <- pathview(gene.data  = geneList,
                      pathway.id = "R-HSA-597592",
                      species    = "hsa",
                      limit      = list(gene=max(abs(geneList)), cpd=1))




edo <- enrichDGN(nde)

dotplot(edo, showCategory=8) + ggtitle(paste0("dotplot for ORA in sig miRNA Target") )


edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=de,  showCategory = 4)
p1



## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=de)
p2


p3 <- cnetplot(edox, foldChange=de, circular = TRUE, colorEdge = TRUE) 
p3
edo <- pairwise_termsim(edo)
q1 <- emapplot(edo, min_edge=0.2)
q2 <- emapplot(edo, cex_category=1.5)
q3 <- emapplot(edo, layout="kk")
q4 <- emapplot(edo, cex_category=1.5,layout="kk") 
q4

edox2 <- pairwise_termsim(edox)
r1 <- treeplot(edox2)


#library(enrichplot)
#library(DOSE)
edo <- enrichDGN(de)


dotplot(edo, showCategory=10) + ggtitle("dotplot for ORA")
#dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList,  showCategory = 4)

## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 

edo <- pairwise_termsim(edo)
q1 <- emapplot(edo, min_edge=0.2)
q2 <- emapplot(edo, cex_category=1.5)
q3 <- emapplot(edo, layout="kk")
q4 <- emapplot(edo, cex_category=1.5,layout="kk") 
q4

edox2 <- pairwise_termsim(edox)
r1 <- treeplot(edox2)

x = enrichDO(de)
y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

z=ggplot(y, showCategory = 20, 
         aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched Disease Ontology")











nde=as.character(names(de))

#edo <- enrichDGN(nde)

#dotplot(edo, showCategory=10, orderBy = "x") + ggtitle(paste0("dotplot for ORA in ",target_name) )



x = enrichDO(nde)
y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

z=ggplot(y, showCategory = 10, 
         aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle(paste0("Top 10 Enriched Disease Ontology in ", substring(target_name,1,nchar(target_name)-6),  " Vs KnTn Significant genes "))
z
#--------------------------------------------------------
library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)
library(ReactomePA)


gsedf= DMGs[DMGs$gene %in% unique_to_target & DMGs$contrast == target_name , c(1,7,8,9) ]

gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")


gsedf=na.omit(gsedf)

gsedf=  gsedf[order(gsedf$rank, decreasing = TRUE), ]


de <- gsedf$rank

names(de)=gsedf$id


geneList=de

#, scoreType = "pos"
x <- gsePathway(geneList)


y <- arrange(x, abs(NES)) %>% 
  group_by(sign(NES)) %>% 
  slice(1:5)


ggplot(y, aes(NES, fct_reorder(Description, NES), fill=qvalues), showCategory=10) + 
  geom_col(orientation='y') + 
  scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
  theme_minimal() + ylab(NULL)
