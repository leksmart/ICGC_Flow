
comp.int=c("KpTp", "KnTn")

library(edgeR)
library(limma)
library(ComplexHeatmap)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(tidyverse)

root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main)))}

setwd(file.path(root, main))



load( "validation/adj.counts030324.RData")

load("DMethDE030324.RData")

load("validation/de.list010524.RData")

load("DNAmet122323.RData")   #load("meth030324.RData")

load("validation/DMGsraw1221.RData") #load( "validation/DMGs030324.RData")

de.genes=DMGs;rm(DMGs)
#ordering diffrentially expressed genes

de.genes<-de.genes[with(de.genes, order(abs(logFC), adj.P.Val, decreasing = TRUE)), ]

comp.int=c("KpTp", "KnTn")

contint= paste0(comp.int[1], " - ", comp.int[2])

de.genes=de.genes[de.genes$contrast == contint,]
rownames(de.genes)=de.genes$gene
de.genes=de.genes[de.genes$adj.P.Val < 0.05,]
# voomObj is normalized expression values on the log2 scale
norm.count <- adj.counts; rm(adj.counts)

load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")
clinical=usidtable;rm(usidtable)
clinical$grp=ifelse(clinical$grp =="KnTnHB", "KnTn", clinical$grp)




rnametid= intersect(intersect(colnames(meth), clinical$icgc_donor_id),colnames(norm.count))


low.g_id <- clinical$icgc_donor_id[clinical$icgc_donor_id %in% rnametid &
                                     clinical$grp == comp.int[2]]

high.g_id <-  clinical$icgc_donor_id[clinical$icgc_donor_id %in% rnametid &
                                       clinical$grp == comp.int[1]]

load("rsd.RData")
norm.count=norm.count[unique(rsd[rsd$View==1, 2]),c(low.g_id,high.g_id)]
meth=meth[unique(rsd[rsd$View==2, 2]),c(low.g_id,high.g_id)]

clinical=clinical[clinical$icgc_donor_id %in% c(low.g_id,high.g_id),]

#______________preparing methylation data for cis-regulatory analysis____________#
# get the 450k annotation data

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


bval= t(apply( meth, 1, function(x) (sin(x) + 1) / 2))


# finding probes in promoter of genes
table(data.frame(ann450k)$Regulatory_Feature_Group) ## to find regulatory features of probes

# selecting a subset of probes associated with promoted
promoter.probe <- rownames(data.frame(ann450k))[data.frame(ann450k)$Regulatory_Feature_Group 
                                                %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")]


# find genes probes with significantly different methylation status in 
# KnTn and KpTp PDAC

# low.g_id <- clinical$icgc_donor_id[clinical$grp == comp.int[2]]
# high.g_id <-  clinical$icgc_donor_id[clinical$grp == comp.int[1]]

dbet <- data.frame (low.grade = rowMeans(bval[, low.g_id]),
                    high.grade = rowMeans(bval[, high.g_id]))

dbet$delta <- abs(dbet$low.grade - dbet$high.grade)

db.probe <- rownames(dbet)[dbet$delta > 0.1] # those with deltabeta > 0.2
#db.probe <- db.probe %in% promoter.probe # those resided in promoter

db.probe<- intersect(db.probe, promoter.probe)

# those genes flanked to promote probe
db.genes <- data.frame(ann450k)[rownames(data.frame(ann450k)) %in% db.probe, ]
db.genes <- db.genes[, c("Name","UCSC_RefGene_Name")]
db.genes <- tidyr::separate_rows(db.genes, Name, UCSC_RefGene_Name) # extending collapsed cells
db.genes$comb <- paste(db.genes$Name,db.genes$UCSC_RefGene_Name) # remove duplicates
db.genes <- db.genes[!duplicated(db.genes$comb), ]
db.genes <- db.genes[, -3]

# doing correlation analysis
# polishing matrices to have only high grade samples
rownames(clinical)=clinical$icgc_donor_id
cis.bval.mat <- bval[, high.g_id]
cis.exp.mat <- norm.count[, rownames(clinical)[clinical$grp == comp.int[1]]]

#cis.exp.mat <- cis.exp.mat[, colnames(cis.bval.mat)]

#editing expression matrix rowname
df <- data.frame(name = row.names(cis.exp.mat)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) # this do magic like "text to column" in Excel!
rowName <- df$name
# find duplicates in rowName, if any
table(duplicated(rowName))
#FALSE  TRUE 
#20530     1 
# in order to resolve  duplucation issue
# rowName[duplicated(rowName) == TRUE]
# #[1] "SLC35E2"
# #
# rowName[grep("SLC35E2", rowName)[2]] <- "SLC35E2_2"
# #setting rna row names 
# row.names(cis.exp.mat) <- rowName
# rm(df, rowName) # removing datasets that we do not need anymore

#__________________correlation analysis_________________________________#

cis.reg = data.frame( gene=character(0), cpg=character(0), pval=numeric(0), cor=numeric(0))

for (i in 1:nrow(db.genes)){
  cpg = db.genes[i,][1]
  gene = db.genes[i,][2]
  if (gene %in% rownames(cis.exp.mat)){
    df1 <- data.frame(exp= cis.exp.mat[as.character(gene), ])
    df2 <- (cis.bval.mat[as.character(cpg), ])
    df <- merge(df1,df2, by = 0)
    res <- cor.test(df[,2], df[,3], method = "pearson")
    pval = round(res$p.value, 4)
    cor = round(res$estimate, 4)
    cis.reg[i,] <- c(gene, cpg, pval, cor)
  }
}


cis.reg$adj.P.Val = round(p.adjust(cis.reg$pval, "fdr"),4)
cis.reg <- cis.reg[with(cis.reg, order(cor, adj.P.Val)), ]

# top pair visualization
# inspecting the results, C2orf74 gene has significant correlation with probes:

gen.vis <- merge(data.frame(exp= cis.exp.mat["MGST1", ]), 
                 t(cis.bval.mat[rownames(cis.bval.mat)[1:50], ]),
                 by = 0)

par(mfrow=c(3,2))
sapply(names(gen.vis)[3:8], function(cpg){
  plot(x= gen.vis[ ,cpg], y = gen.vis[,2], xlab = "beta value",
       xlim = c(0,1),
       ylab = "normalized expression" ,
       pch = 19,
       main = paste("MGST1",cpg, sep = "-"),
       frame = FALSE)
  abline(lm(gen.vis[,2] ~ gen.vis[ ,cpg], data = gen.vis), col = "blue")
})




#--------------------- trans-regulation visualization--------------------------
# adding genes to delta beta data 
tran.reg <- data.frame(ann450k)[rownames(data.frame(ann450k)) %in% rownames(dbet), ][, c(4,24)]
tran.reg <- tidyr::separate_rows(tran.reg, Name, UCSC_RefGene_Name) # extending collapsed cells
tran.reg$comb <- paste(tran.reg$Name,tran.reg$UCSC_RefGene_Name) # remove duplicates
tran.reg <- tran.reg[!duplicated(tran.reg$comb), ]
tran.reg <- tran.reg[, -3]
names(tran.reg)[2] <- "gene"

# merging with deltabeta dataframe
dbet$Name <- rownames(dbet)
tran.reg <- merge(tran.reg, dbet, by = "Name")
# joining with differential expression analysis result
#editing expression matrix rowname
df <- data.frame(name = row.names(de.genes)) # keeping rownames as a temporary data frame
df <- data.frame(do.call('rbind', strsplit(as.character(df$name),'|',fixed=TRUE))) # this do magic like "text to column" in Excel!
df$X1[df$X1 == "?"] <- df$X2 # replace "? with entrez gene number
rowName <- df$X1
# find duplicates in rowName, if any
#table(duplicated(rowName))
#FALSE  TRUE 
#16339     1    
# in order to resolve  duplication issue
# rowName[duplicated(rowName) == TRUE]
# grep("SLC35E2", rowName)
# #[1]  9225 15546
# rowName[15546] <- "SLC35E2_2"
# #setting rna row names 
# row.names(de.genes) <- rowName
# rm(df, rowName) # removing datasets that we do not need anymore
de.genes$rownames.exp_mat. <- rownames(de.genes)
#names(de.genes)[8] <- "gene"
# merging
tran.reg <- merge(tran.reg, de.genes, by.x = "gene", by.y="gene")
# inspecting data
hist(tran.reg$logFC)
 # delta was calculated as abs(delta), re-calculate to have original value
tran.reg$delta <- tran.reg$high.grade - tran.reg$low.grade
hist(tran.reg$delta)
# defining a column for coloring
tran.reg$group <- ifelse(tran.reg$delta <= -0.1 & tran.reg$logFC <= -0.5, "hypo-down",
                         ifelse(tran.reg$delta <= -0.1 & tran.reg$logFC >= 0.5, "hypo-up",
                                ifelse(tran.reg$delta >= 0.1 & tran.reg$logFC <= -0.5, "hypr-down",
                                       ifelse(tran.reg$delta >= 0.1 & tran.reg$logFC >= 0.5, "hypr-up", "not-sig"))))

#tran.reg2=tran.reg[tran.reg$gene %in% rsd1$`Variable Name`,]

#intersect(tran.reg2$Name, rsd2$`Variable Name`)


# plotting
cols <- c("hypo-down" = "#B8860B", "hypo-up" = "blue", "not-sig" = "grey", "hypr-down" = "red", "hypr-up" = "springgreen4")

ggplot(tran.reg, aes(x = delta, y = logFC, color = group)) +
  geom_point(size = 2.5, alpha = 1, na.rm = T) +
  scale_colour_manual(values = cols) + 
  theme_bw(base_size = 14) +
  geom_hline(yintercept = 1.0, colour="#990000", linetype="dashed") + 
  geom_hline(yintercept = -1.0, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = 0.1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = -0.1, colour="#990000", linetype="dashed") +
  xlab("mean methylation differences") + 
  ylab("Log2 expression change") +
  ggtitle(paste0(comp.int[2], " Vs ", comp.int[1]))


#################---------------------------------##################################

mval= meth[unique(rsd[rsd$View==2, 2]),c(low.g_id,high.g_id)]
myAnnotation <- cpg.annotate(object = mval, datatype = "array", 
                             what = "M", 
                             analysis.type = "variability", 
                             arraytype = "450K",
                              fdr = 0.05)
# 
DFT=as.data.frame(myAnnotation@ranges)
DFT$probes= rownames(DFT)
str(myAnnotation)
#coef = "paper_Histologic.gradeHigh Grade"
# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)

dmr.tableKpTp <- data.frame(results.ranges)


#############

# Load necessary library
library(curl)

# URL of the file
url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz"

# Download the file
temp <- tempfile()
curl_download(url, temp)
# Unzip the file
data <- (readLines(gzfile(temp)))
# Split each line into separate elements
data <- lapply(strsplit(data, "\t"), as.character)

# Convert to dataframe
cpgdf <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)

colnames(cpgdf)= c("bin",	"chrom",	"chromStart",	"chromEnd"	,"name",
                   "length",	"cpgNum",	"gcNum",	"perCpg"	,"perGc",	"obsExp")

# Clean up
unlink(temp);rm(data)

cpgdf=cpgdf[,-1]
cpgdf=cpgdf[cpgdf$chrom == "chr19", ]
cpgdf$chromStart=as.numeric(cpgdf$chromStart)
cpgdf$chromEnd=as.numeric(cpgdf$chromEnd)


islandData <- GRanges(seqnames=Rle(cpgdf[,1]), 
                      ranges=IRanges(start=cpgdf[,2],
                                     end=cpgdf[,3]),
                      strand=Rle(strand(rep("*",nrow(cpgdf)))))



url= "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeAwgDnaseMasterSites.txt.gz"


# Download the file
temp <- tempfile()
curl_download(url, temp)
# Unzip the file
data <- (readLines(gzfile(temp)))
# Split each line into separate elements
data <- lapply(strsplit(data, "\t"), as.character)

# Convert to dataframe
dnaseidf <- as.data.frame(do.call(rbind, data), stringsAsFactors = FALSE)


colnames(dnaseidf)= c("bin",	"chrom",	"chromStart",	"chromEnd"	,"name",
                   "length",	"cpgNum",	"gcNum",	"perCpg"	,"perGc",	"obsExp")

# Clean up
unlink(temp);rm(data)