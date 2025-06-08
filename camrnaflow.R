library(tidyverse)
library(edgeR)
library(limma)
library(ComplexHeatmap)
#library(annotables)
library(org.Hs.eg.db)

getwd()
root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"
subf="validation"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main,subf)))}

setwd(file.path(root, main, subf))

coho= 2
metadata=read.csv("dataset/mmc2.csv")[1:319,]
metadata=metadata[metadata$tissue.source=="surgical resection, fresh frozen" & 
                    metadata$cohort=="unpaired",]


metadata=metadata[!is.na(metadata$RNASeq.id),-13]

usmrna <- "dataset/PAAD-US_ICGC_rnaseq.rds"
camrna <-  "dataset/PACA-CA_ICGC_rnaseq.rds"
aumrna<- "dataset/PACA-AU_ICGC_rnaseq.rds"


mrna=c(usmrna,camrna,aumrna)

rnaseq =readRDS(file.path(root, main, subf,mrna[coho]))

rnaseq=rnaseq[rnaseq$submitted_sample_id %in% metadata$RNASeq.id,]

rnaseq=rnaseq[rnaseq$raw_read_count >1,]

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
#rnaseqt=rnaseq[,1:177]

#rnaseqt=t(rnaseqt)

#rnaseqt=as.data.frame(rnaseqt)
#write_csv(rnaseqt, file = "carnaseq.csv", col_names=T)
#compute purity using puree


puree=read_tsv("carnasepurity.txt")

puree$id=colnames(rnaseq)[1:177]

puree=as.data.frame(puree)

rownames(puree)=puree$id
puree=puree[,c(3,2)]

load("dataset/caidtable.RData")

caidpuree= merge(puree, caidtable, by="id")

caidpuree=caidpuree[caidpuree$purity >= 0.7,]

normid = colnames(rnaseq)[substring(colnames(rnaseq),nchar(colnames(rnaseq)),nchar(colnames(rnaseq)))=="N"]

normtable= caidpuree[1:3,]
normtable$chemo=""
normtable$id=normid
normtable$grp="WT"
normtable$time=NA
normtable$k.status=NA
normtable$status=NA

caidpuree= rbind(normtable, caidpuree)


rnaseq=as.data.frame(rnaseq)

#rnaseq =  rnaseq[rnaseq < 0] <- 0
normet=rnaseq[,normid]
tum=rnaseq[,setdiff(colnames(rnaseq),normid)]

sub.id=intersect(colnames(tum),caidpuree$id )


tum=tum[,sub.id]



purity=caidpuree[caidpuree$id %in% sub.id,]
rownames(purity)= purity$id
Tumor.mat = tum
Norm.mat = normet


Norm.matt = rowMeans(Norm.mat)
lamda = as.numeric(purity[colnames(Tumor.mat), 2])
names(lamda)=colnames(Tumor.mat)
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


rnaseq=cbind(Tumor.mat,Norm.mat)


commonid= intersect(caidpuree$id, colnames(rnaseq))

caidpuree = caidpuree[caidpuree$id %in%  commonid, ]

rownames(caidpuree)= caidpuree$id

rnaseq=rnaseq[,commonid]


caidpuree$grp=as.factor(caidpuree$grp)

library(DESeq2)
rnaseq[is.na(rnaseq)] <- 0
rnaseq[rnaseq < 0] <- 0


corfc = abs(min(rnaseq))
rnaseq =rnaseq + corfc
rnaseq=round(rnaseq)
dds <- DESeqDataSetFromMatrix(countData=rnaseq, colData=caidpuree, design = ~ grp )

features =  c("TMOD3","SLC9A3", "BTBD6", "FASN",  "CGB7", "TLX1")

srt= round(3)
dds75 <- dds[rowSums(counts(dds) >= (10+corfc)) >= srt,]
nrow(dds)
nrow(dds75) 


ddst <- DESeq(dds75)
vsd <- vst(ddst, blind=FALSE)
norm.counts <- assay(vsd)

cal.val.rna=norm.counts[features,]

#save( cal.val.rna, file="cal.val.rna.RData")

normalized_counts <- counts(dds, normalized=TRUE)
#compute tumor purity from outside



adj.counts=norm.counts

sub.id=intersect(colnames(adj.counts), caidpuree$id)

group <-  as.data.frame(caidpuree[caidpuree$id %in% sub.id, c("id","grp")])

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

#-------------------------------------------------------------------------------------------------
load("validation/cal.val.rna.RData")
load("validation/dataset/caidtable.RData")
caidtable=as.data.frame(caidtable)
features =  c("TMOD3","SLC9A3", "BTBD6", "FASN", "LOC647121", "CGB7", "TLX1")
# Define the list of variables
rownames(caidtable)=caidtable$id

goi <- features[-5]


library(survminer)
library(survival)
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
mem.traits.phen= t(cal.val.rna)
mem.traits.phenca=merge(mem.traits.phen,caidtable, by=0)
mem.traits.phenca= as.data.frame(mem.traits.phenca)
modelph <- coxph(formula_str, data = mem.traits.phenca)

summary(modelph)


mem.traits.phen=rbind(mem.traits.phenca, mem.traits.phenau)

modelph <- coxph(formula_str, data = mem.traits.phen)

summary(modelph)

