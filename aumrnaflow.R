library(tidyverse)
library(edgeR)
library(limma)
library(ComplexHeatmap)
#library(annotables)
library(org.Hs.eg.db)
library(DESeq2)
getwd()
root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"
subf="validation"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main,subf)))}

setwd(file.path(root, main, subf))

coho= 3

#metadata=metadata[!is.na(metadata$RNASeq.id),-13]

usmrna <- "dataset/PAAD-US_ICGC_rnaseq.rds"
camrna <-  "dataset/PACA-CA_ICGC_rnaseq.rds"
aumrna<- "dataset/PACA-AU_ICGC_rnaseq.rds"


mrna=c(usmrna,camrna,aumrna)

rnaseq =readRDS(file.path(root, main, subf,mrna[coho]))

rnaseq=rnaseq[rnaseq$submitted_sample_id %in% metadata$RNASeq.id,]

rnaseq=rnaseq[rnaseq$raw_read_count >=1,]

normseq <- readRDS(file.path(root, "dataset/PAAD-US_ICGC_rnaseq.rds"))

normseq $class=substring(normseq$submitted_sample_id, 14,16)

normseq$icgc_donor_id =ifelse(normseq$class  == "11A", paste0(normseq$icgc_donor_id,"N"), normseq$icgc_donor_id)

normseq=normseq[normseq$class=="11A",]

normseq=normseq[,-23]

data = rnaseq[,"gene_id"]

data=unique(data)

library( "biomaRt" ) #example code for mouse gene id mapping
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id",values = data, mart = ensembl )


annots=genemap 

#rnaseq <- cbind(rnaseq, annots$SYMBOL)

annots=na.omit(annots)
annots=unique(annots)

ass=as.data.frame(table(annots$hgnc_symbol))

ass= ass[ass$Freq<2,]

annots=annots[annots$hgnc_symbol %in% ass$Var1,]

colnames(annots)[1]="gene_id"

rnaseq=rnaseq[rnaseq$gene_id %in% annots$gene_id,]


rnaseq=merge(rnaseq, annots, by="gene_id")


rnaseq=rnaseq[,-1]

colnames(rnaseq)[22]="gene_id" 

colnm= colnames(normseq)

rnaseq=rnaseq[,colnm]
rnaseq=rbind(rnaseq, normseq)

rnaseqbk=  rnaseq

rnaseq=rnaseq %>% dplyr::select(c(1,8,10)) %>% filter(!gene_id %in% c("?", "SLC35E2"))

rnaseq= tidyr::pivot_wider(rnaseq, id_cols= "gene_id", names_from = "icgc_donor_id" , values_from = "raw_read_count", values_fn= mean)

rnaseq=as.data.frame(rnaseq)

rownames(rnaseq)= rnaseq$gene_id

rnaseq=rnaseq[,-1]

rnaseq[is.na(rnaseq)] <- 0

#compute tumor purity from outside
#rnaseqt=rnaseq[,1:91]

#rnaseqt=t(rnaseqt)

#rnaseqt=as.data.frame(rnaseqt)
#write_csv(rnaseqt, file = "aurnaseq.csv", col_names=T)
#compute purity using puree

#purronam= rownames(rnaseqt)
#save(purronam, file = "camrnanames.RData")
# puree=read_tsv("aupurity.txt")
# 
# puree$id=purronam
# 
# puree=as.data.frame(puree)
# 
# rownames(puree)=puree$id
#puree=puree[,c(3,2)]

#save(puree, file = "aupuree.RData")
#save(rnaseq, file = "aurrnaseq.RData")
#----------------------------------------------------------------------------------------

load("aurrnaseq.RData")

load("aupuree.RData")

load("dataset/auidtable.RData")

auidpuree= merge(puree, auidtable, by="id")

auidpuree=auidpuree[auidpuree$purity >= 0.4,]

normid = colnames(rnaseq)[substring(colnames(rnaseq),nchar(colnames(rnaseq)),nchar(colnames(rnaseq)))=="N"]

normtable= auidpuree[1:3,]
normtable$chemo=""
normtable$id=normid
normtable$grp="WT"
normtable$time=NA
normtable$k.status=NA
normtable$status=NA

auidpuree= rbind(normtable, auidpuree)


rnaseq=as.data.frame(rnaseq)

#rnaseq =  rnaseq[rnaseq < 0] <- 0


commonid= intersect(auidpuree$id, colnames(rnaseq))

auidpuree = auidpuree[auidpuree$id %in%  commonid, ]

rownames(auidpuree)= auidpuree$id

rnaseq=rnaseq[,commonid]


auidpuree$grp=as.factor(auidpuree$grp)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=round(rnaseq), colData=auidpuree, design = ~ grp )


features =  c("TMOD3","SLC9A3", "BTBD6", "FASN", "LOC647121", "CGB7", "TLX1")

srt= round(1)
dds75 <- dds[rowSums(counts(dds) >= 1) >= srt,]
nrow(dds)
nrow(dds75) 
# 
ddst <- DESeq(dds75)
vsd <- vst(ddst, blind=FALSE)
norm.counts <- assay(vsd)

au.val.rna=norm.counts[features,]
#save( au.val.rna, file="au.val.rna.RData")

normalized_counts <- counts(ddst, normalized=TRUE)
#compute tumor purity from outside



adj.counts=norm.counts



sub.id=intersect(colnames(adj.counts), auidpuree$id)

group <-  as.data.frame(auidpuree[auidpuree$id %in% sub.id, c("id","grp")])

rownames(group)=group$id

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



# Define the list of variables
rownames(auidtable)=auidtable$id
goi <- features[-5]
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
mem.traits.phen= t(au.val.rna)


mem.traits.phen=merge(mem.traits.phen,auidtable, by=0)


modelph <- coxph(formula_str, data = mem.traits.phen)

summary(modelph)


#-------------------------------------------------------------------------------------------------
  load("validation/au.val.rna.RData")
  load("validation/dataset/auidtable.RData")
  
 load("AUDNAmval.RData")
  goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" ,"cg06952671" ,"cg10661002")
  
  
 auidtable=as.data.frame(auidtable)
  
  # Define the list of variables
  rownames(auidtable)=auidtable$id

goi <- goi[-5]
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
mem.traits.phen= t(au.val.rna)


mem.traits.phenau=merge(mem.traits.phen,auidtable, by=0)

rownames(mem.traits.phenau)= mem.traits.phenau$Row.names


mem.traits.phenau=mem.traits.phenau[,-1]

mem.traits.phenau=merge(mem.traits.phenau,mval, by=0)

library(survminer);library(survival)
modelph <- coxph(formula_str, data = mem.traits.phenau[mem.traits.phenau$time >=20,])

summary(modelph)



library(factoextra)

df= mem.traits.phenau[,goi]
# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df, centers = 3)

# Get the cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster


mem.traits.phenau$subtype=cluster_assignments

mem.traits.phenau$subtype = ifelse(mem.traits.phenau$subtype == 1, "Group_A",
                                    ifelse(mem.traits.phenau$subtype == "2", "Group_B", "Group_C"))


sfit <- survfit(Surv(time, status)~subtype, data=mem.traits.phenau[mem.traits.phenau$time >20,])

ggsurvplot(
  sfit,                     # survfit object with calculated statistics.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,         # show confidence intervals for 
  # point estimaes of survival curves.
  #conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 200,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = F,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  #legend.labs = 
  #c("Low Risk", "High Risk"),    # change legend labels.
  palette =  c("#F8766D","#00B0F6",  "#E76BF3", "#E96BF3" ) )

mem.traits.phen$subtype=as.factor(mem.traits.phen$subtype)

mem.traits.phen <- within(mem.traits.phen, subtype <- relevel(subtype, ref = "Group_A"))

