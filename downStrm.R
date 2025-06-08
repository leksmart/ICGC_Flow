root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"

setwd(file.path(root, main))

load( "validation/adj.counts1221.RData")

load( "validation/DMGsraw1221.RData")
tt <- DE
target_name= "KpTp - WT"
filnm= "KpTpVsWT"
#get the indices of scored dataset that have FDR < 0.05

tt=tt[tt$contrast == target_name, ]

select_genes = which(tt$logFC >= 0 & tt$adj.P.Val < 0.05 )
tt=tt[select_genes,]

#output how many genes there are in the set that have FDR < 0.05
length(select_genes)

#gene names from the TCGA set contain gene name and entrez gene ids separated by ‘|’
# for all subsequent enrichment analysis we need to have just one id.  Separate the names 
# into their two ids and keep the gene symbols
topgenes_qvalue005 <- tt[select_genes,7]

#output the top 5 entries in the list of top genes
head(topgenes_qvalue005)

#write results out to the file.  This is an example of a set that can be used for
# Protocol 1

working_dir= file.path(getwd(),"validation")
write.table( topgenes_qvalue005, 
            file.path(working_dir, paste0(filnm,"RNAseq_allsignificantgenes.txt")), 
            col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)



##-----------------------------------------------------------------------------------------


#calculate ranks
ranks_RNAseq = sign(tt$logFC) * -log10(tt$adj.P.Val)

#get gene id



tt$geneids= mapIds(org.Hs.eg.db, keys = tt$gene,
                 column = "ENTREZID", keytype = "SYMBOL")

genenames <- tt$gene
geneids <- tt$geneids

#create ranks file
ranks_RNAseq <- cbind(genenames, ranks_RNAseq)
colnames(ranks_RNAseq) <- c("GeneName","rank")

#sort ranks in decreasing order
ranks_RNAseq <- ranks_RNAseq[order(as.numeric(ranks_RNAseq[,2]),decreasing = TRUE),]

write.table(ranks_RNAseq, file.path(working_dir,
                                    paste0(filnm,"Supplementary_Table2_ranks.rnk")), 
            col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)


head(ranks_RNAseq)
tail(ranks_RNAseq)



#fix issue with biomart not working because of url redirection
options(RCurlOptions=list(followlocation=TRUE, postredir=2L))

normalized_expression_RNAseq <- adj.counts


EM_expressionFile_RNAseq <- data.frame(Name = rownames(adj.counts), normalized_expression_RNAseq)
rownames(EM_expressionFile_RNAseq) <- rownames(normalized_expression_RNAseq)

#Add descriptions instead of geneids
tryCatch(expr = { library("biomaRt")}, 
         error = function(e) { 
           source("https://bioconductor.org/biocLite.R")
           biocLite("biomaRt")}, 
         finally = library("biomaRt"))
mart = useMart(biomart="ENSEMBL_MART_ENSEMBL")
mart = useDataset(mart, dataset="hsapiens_gene_ensembl" )

genes = getBM(attributes = c( 'hgnc_symbol', 'description'), filters='hgnc_symbol', 
              values=genenames, mart=mart);
genes$description = gsub("\\[Source.*", "", genes$description);

EM_expressionFile_RNAseq <- merge(genes,EM_expressionFile_RNAseq,  
                                  all.y=TRUE,by.x=1, by.y=1)
colnames(EM_expressionFile_RNAseq)[1] <- "Name"
colnames(EM_expressionFile_RNAseq)[2] <- "Description"

write.table(EM_expressionFile_RNAseq, 
            file.path(working_dir,
                      "Supplementary_Table6_TCGA_PDAC_RNAseq_expression.txt"),
            col.name=TRUE,sep="\t", row.names=FALSE, quote=FALSE)





#write out a GSEA classes file. (optional)
data_classes <- "grp"
Clinical=readRDS("clinical_New.rds")

classDefinitions_RNASeq=Clinical[Clinical$id %in% colnames(adj.counts),]





fileConn <- file(
  file.path(working_dir,"Supplementary_Table7_TCGA_PDAC_RNAseq_classes.cls"))

writeLines(c(paste(length(classDefinitions_RNASeq[,data_classes]), "5 1"), 
             paste("# ", unique(classDefinitions_RNASeq[,data_classes])[1], " ",
                   unique(classDefinitions_RNASeq[,data_classes])[2], " ",
                   unique(classDefinitions_RNASeq[,data_classes])[3], " ",
                   unique(classDefinitions_RNASeq[,data_classes])[4], " ",
                   unique(classDefinitions_RNASeq[,data_classes])[5])), fileConn)


write.table(t(classDefinitions_RNASeq[,data_classes]), 
            file.path(working_dir,"Supplementary_Table7_TCGA_PDAC_RNAseq_classes.cls"), 
            col.name=FALSE, sep="\t",
            row.names=FALSE, quote=FALSE, append=TRUE)
close(fileConn)



#-------------------PEA-------------------------------------------

tryCatch(expr = { library("limma")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("limma")}, 
         finally = library("limma"))

tryCatch(expr = { library("GSA")}, 
         error = function(e) { source("https://bioconductor.org/biocLite.R")
           biocLite("GSA")}, 
         finally = library("GSA"))

tryCatch(expr = { library("RCurl")}, 
         error = function(e) { 
           install.packages("RCurl")}, 
         finally = library("RCurl"))

# This protocol can use RNA-seq expression data or microarray expression data. 
#Specify which you would like to use
dataType_rnaseq <- TRUE

#The field in the class definition file that defines the classes of the data.
data_classes <- "grp"

#string to name the analysis
analysis_name <-  "KpTpVsKnTn"

#from Supplementary protocol 1
expression_file <- "Supplementary_Table6_TCGA_PDAC_RNAseq_expression.txt"


gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"

#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include terms 
#inferred from electronic annotations(IEA)
#start with gmt file that has pathways only
rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(working_dir,paste("Supplementary_Table3_",gmt_file,sep="") )

download.file(
  paste(gmt_url,gmt_file,sep=""),
  destfile=dest_gmt_file
)




#if you haven't automatically downloaded the gmt file set 
#the path to the gmt file below.
gmt_file <- dest_gmt_file

capture.output(
  genesets <- GSA.read.gmt(gmt_file),
  file="./gsea_load_output.txt"
)
names(genesets$genesets) <- genesets$geneset.names




classDefinitions <-data.frame(classDefinitions_RNASeq[,c(1,5)])
rownames(classDefinitions)= classDefinitions$id

colnames(classDefinitions)[1]= "patient"

colnames(classDefinitions)[2]= "data_classes"

cid = intersect(colnames(normalized_expression_RNAseq),classDefinitions$patient)

identical(colnames(normalized_expression_RNAseq), classDefinitions$patient)

minimalSet <- ExpressionSet(assayData=normalized_expression_RNAseq[,cid])
classDefinitions=classDefinitions[cid,]



#create minimal set

classes <- factor(classDefinitions[,"data_classes"])

#create model
modelDesign <- model.matrix(~ 0 + classes)

#assign data set for the analysis
data_for_gs_analysis <- minimalSet

genesets_filtered <- ids2indices(genesets$genesets, rownames(data_for_gs_analysis), 
                                 remove.empty=TRUE)


geneset_sizes <- unlist(lapply(genesets_filtered, length))
geneset_indices <- which(geneset_sizes>=10 & geneset_sizes<500)




design <- model.matrix(~ 0 + classes)

contrast_kptpvskntn <- makeContrasts(
  kptpvskntn ="classesKpTp-classesKnTn",levels=design)




#Run ROAST
mroast_results <- mroast(data_for_gs_analysis, genesets_filtered[geneset_indices],
                         design,contrast=contrast_kptpvskntn, nrot=10000)



mroast_descr <- unlist(lapply(rownames(mroast_results), 
                              function(x){unlist(strsplit(x,"\\%"))[1]}))




mroast_results_file <- "mroast_results_generic_em.txt"

Phenotype <- unlist(lapply(mroast_results[,"Direction"],function(x)
{if(x=="Up"){1}else{(-1)}}))
genes <- c()
for(i in 1:length(rownames(mroast_results))){
  current_geneset <- unlist(genesets_filtered
                            [ which( names(genesets_filtered) %in% rownames(mroast_results)[i])])
  current_genes <- c()
  for(j in 1:length(current_geneset)){
    if(j==length(current_geneset)){
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]], 
                             sep="")
    } else {
      current_genes <- paste(current_genes, 
                             rownames(data_for_gs_analysis)[current_geneset[j]],
                             ",", sep="")
    }
  }
  genes <- rbind(genes, current_genes)
}
rownames(genes) <- rownames(mroast_results)

mroast_results_generic_em <- data.frame( rownames(mroast_results), mroast_descr, 
                                         PValue=mroast_results[,"PValue"], FDR=mroast_results[,"FDR"], Phenotype, genes)
write.table(mroast_results_generic_em, file.path(working_dir,mroast_results_file), 
            col.name=TRUE, sep="\t", row.names=FALSE, quote=FALSE) 



##--------------------CAMARA-----------------------------------------------------------------#

camera_results_file <- "camera_results_generic_em.txt"



camera_results <- camera(data_for_gs_analysis, 
                         genesets_filtered[geneset_indices], design, contrast=contrast_kptpvskntn)
camera_descr <- unlist(lapply(rownames(camera_results), 
                              function(x){unlist(strsplit(x,"\\%"))[1]}))
camera_Phenotype <- unlist(lapply(camera_results[,"Direction"], 
                                  function(x){if(x=="Up"){1}else{(-1)}}))


camera_genes <- c()
for(i in 1:length(rownames(camera_results))){
  current_geneset <- unlist( 
    genesets_filtered[ which( names( genesets_filtered ) %in% 
                                rownames(camera_results)[i])])
  current_genes <- c()
  for(j in 1:length(current_geneset)){
    if(j==length(current_geneset)){
      current_genes <- paste( current_genes, 
                              rownames(data_for_gs_analysis) [current_geneset[j]],
                              sep="")
    } else {
      current_genes <- paste( current_genes, 
                              rownames(data_for_gs_analysis)[ current_geneset[j]], ",", 
                              sep="")
    }
  }
  camera_genes <- rbind(camera_genes, current_genes)
}
rownames(camera_genes) <- rownames(camera_results)

camera_results_generic_em <- data.frame(rownames(camera_results), camera_descr, 
                                        PValue = camera_results[,"PValue"], FDR=camera_results[,"FDR"], Phenotype, genes )
write.table(camera_results_generic_em, file.path(working_dir,camera_results_file), 
            col.name=TRUE, sep="\t", row.names=FALSE, quote=FALSE)





#library(org.Hs.eg.db)


#library(enrichplot)
#library(DOSE)

#library(clusterProfiler)
#library(GOSemSim)

library(forcats)
library(ggplot2)
library(ggstance)
library(enrichplot)
library(ReactomePA)

load("validation/rsd.RData")
target_name= "KpTp - KnTn"
gsedf=rsd[rsd$grp== target_name & rsd$View ==1,]

colnames(gsedf)[2]="gene"
gsedf$id= mapIds(org.Hs.eg.db, keys = gsedf$gene,
                 column = "ENTREZID", keytype = "SYMBOL")


gsedf=na.omit(gsedf)

gsedf=  gsedf[order(gsedf$`Normalized Relative Importance`), ]


de <- gsedf$`Normalized Relative Importance`

names(de)=gsedf$id


nde=as.character(names(de))

edo <- enrichDGN(nde)

dotplot(edo, showCategory=10, orderBy = "x") + ggtitle(paste0("dotplot for ORA in ", target_name ))



x = enrichDO(nde)
y <- mutate(x, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

library(ggplot2)
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
