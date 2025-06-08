library(edgeR)
library(knitr)
library(limma)
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

getwd()
root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main)))}

setwd(file.path(root, main))

message((paste0("Working directory is set successfully at ", getwd())))

junck.list =  c()

broom = function(){
  
  for(file_name in junck.list){
    if (file.exists(file_name)) {
      unlink(file_name)
      print(paste0(file_name," deleted.."))
    } else{
      print(paste0(file_name," not exists.."))
    }
  }
  
}


DNAMeth<- readRDS(file.path(root, "dataset/PAAD-US_ICGC_DNAmeth.rds"))
#DNAMethb=DNAMeth
DNAMeth$class=substring(DNAMeth$submitted_sample_id, 14,16)

DNAMeth$icgc_donor_id =ifelse(DNAMeth$class  == "11A", paste0(DNAMeth$icgc_donor_id,"N"), DNAMeth$icgc_donor_id)

#normDNAMeth=  DNAMeth %>% filter(class=="11A") %>% pull(icgc_donor_id) %>% unique()

DNAMeth <- DNAMeth[,c("icgc_donor_id", "probe_id","methylation_value")]

dnaMet= tidyr::pivot_wider(DNAMeth, id_cols= "probe_id", names_from = "icgc_donor_id" , 
                           values_from = "methylation_value", values_fn=max)
dnaMet=as.data.frame(dnaMet)
rownames(dnaMet)=dnaMet$probe_id
dnaMet=dnaMet[,-1]

## remove probes with NA
probe.na <- rowSums(is.na(dnaMet))

probe <- probe.na[probe.na == 0]
met <- dnaMet[row.names(dnaMet) %in% names(probe), ]

dim(met)
dim(dnaMet)


# get the 450k annotation data
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)


## remove probes that match to chromosome  X and Y 
keep <- !(row.names(met) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
table(keep)
met <- met[keep, ]
rm(keep) # remove no further needed probes.

## remove SNPs overlapped probe
table(is.na(ann450k$Probe_rs))
# probes without snp
no.snp.probe <- ann450k$Name[is.na(ann450k$Probe_rs)]

snp.probe <- ann450k[!is.na(ann450k$Probe_rs), ]
#snps with maf <= 0.05
snp5.probe <- snp.probe$Name[snp.probe$Probe_maf <= 0.05]

# filtre met
met <- met[row.names(met) %in% c(no.snp.probe, snp5.probe), ]

crs.reac <- read.csv("cross_reactive_probe.chen2013.csv")
crs.reac <- crs.reac$TargetID[-1]

# filtre met
met <- met[ -which(row.names(met) %in% crs.reac), ]


####---------------START HERE-------################################
#save(met, file = "usmetinput.RData")

load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/usmetinput.RData")

load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/usidtable.RData")

normid = colnames(met)[substring(colnames(met),nchar(colnames(met)),nchar(colnames(met)))=="N"]
normet=met[,normid]

met=met[,setdiff(colnames(met),normid)]

#save(normet, file="normet.RData")
 library(InfiniumPurify)
puree<-InfiniumPurify::getPurity(tumor.data=met, tumor.type="PAAD")
puree=as.data.frame(puree)

 puree$id=rownames(puree)
# 
 puree=as.data.frame(puree)
# 
puree=puree[,c(2,1)]

colnames(puree)[2]= "purity"
#save(purit,file="/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/methpurity.RData")
#load("expurity.RData")
#purit=puree
# 
# 
# puree=read_tsv("validation/us_rnaseq_purity.txt")
# 
# puree=as.data.frame(puree)
# 
# mypcoln=readRDS("validation/uspuritycolnms.tsv")
# 
# rownames(puree)= mypcoln
# 
# puree$id=rownames(puree)
# 
# puree=as.data.frame(puree)
# 
# puree=puree[,c(3,2)]


purit=puree; 

Tumor.mat = met
Norm.mat = normet
purit=purit[,2]
names(purit)=puree$id
DMC=InfiniumDMC(tumor.data=Tumor.mat,normal.data=Norm.mat,purity=purit)

Norm.mat1 = Norm.mat
colnames(Norm.mat1) = paste0(colnames(Norm.mat1),1)
Norm.mat2 = Norm.mat
colnames(Norm.mat2) = paste0(colnames(Norm.mat2),2)
Norm.mat3 = Norm.mat
colnames(Norm.mat3) = paste0(colnames(Norm.mat3),3)

Normdf= cbind(Norm.mat2, cbind(Norm.mat1,Norm.mat3))

purifiedmval= InfiniumPurify(tumor.data=Tumor.mat,normal.data=Normdf,purity=purit)

bval= purifiedmval

mval = t(apply(bval, 1, function(x) asin(x)))


#save(mval, file="pureMval.RData")

#save(bval, file="pureBval.RData")

#save(Normdf, file="norBval.RData")

#Mvalnormdf= t(apply(Normdf, 1, function(x) asin(x)))

#Mvalnormdf= t(apply(Normdf, 1,  function(x)asin(2*x - 1))

#save(Mvalnormdf, file="norMval.RData")


###############Start Here#################################################

#  Norm.mat=as.matrix(Norm.mat)
#  Norm.mat= t(apply( Norm.mat, 1, function(x) asin(x)))
#  
#  Norm.matt = rowMeans(Norm.mat)
#  
# 
#  
#  
#  ids= intersect(names(purit), colnames(Tumor.mat))
#  lamda = purit[ids]
# lamda=na.omit(lamda)
#  Tumor.mat=Tumor.mat[,ids]
# 
# # # Initialize the result matrix J
# 
# J = matrix(0, nrow(Tumor.mat), ncol(Tumor.mat))
# rownames(J) = rownames(Tumor.mat)
# colnames(J) = colnames(Tumor.mat)
# 
# # Compute J
# 
# 
# for (c in 1:ncol(J)) {
# 
#   for (r in 1:nrow(J)) {
# 
#     normexp.i=Norm.matt[r]
#     puri.i=1-lamda[c]
#     nomexp=puri.i*normexp.i
#     admx=Tumor.mat[r, c]
#     ym= admx-nomexp
#     realexp=ym/lamda[c]
# 
#     J[r, c] = realexp
#   }
# }


meth= cbind(mval,Mvalnormdf)

tumtable=usidtable[!(usidtable$icgc_donor_id %in% normid),]

normtable=usidtable[1:ncol(Mvalnormdf),]
normtable$icgc_donor_id=colnames(Mvalnormdf)

normtable$grp="WT"


Clinical=rbind(tumtable, normtable)

Clinical= Clinical[,c(1,6)]
Clinical= unique(Clinical)
rownames(Clinical)= Clinical$icgc_donor_id
Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn", Clinical$grp)

#purit=rbind(purit, data.frame(id=normid, purity= rep(0,length(normid))))
#rownames(purit)= purit$id

#Clinical= merge(Clinical, purit, by=0)
sub.id=colnames(meth)
sub.id=intersect(Clinical$icgc_donor_id, sub.id)

#sub.id=setdiff(sub.id,excl)

meth=meth[,sub.id]
Clinical=Clinical[Clinical$icgc_donor_id %in% sub.id,]
setequal(colnames(meth), Clinical$icgc_donor_id)


#mval=t(apply(meth, 1, function(x) log2(x/(1-x))))
mval=meth
#mval[is.nan(mval)] <- 0
#save(mval, file="DNAval021624.RData")
#save(mval, file="DNAval021624.RData")

#save(mval, file="DNAval131.RData")
#save(meth, file="DNAmet131.RData")

#save(mval, file="DNAval2224.RData")
#save(meth, file="DNAmet2224.RData")
#save(meth, file="DNAmet20624.RData")
#save(mval, file="DNAval020624.RData")

#Making grouping variable
Clinical$grp <- as.factor(Clinical$grp)
#levels(clinical$paper_Histologic.grade)
Clinical$grp<- relevel(Clinical$grp, ref = "WT")

#_____________ DMC analysis________________#

colnames(Clinical)[1]= "icgc_donor_id"

group <-  as.data.frame(Clinical[Clinical$icgc_donor_id %in% sub.id, c("icgc_donor_id","grp")])

rownames(group)=group$icgc_donor_id
group=group[sub.id,2]
group=as.factor( group)


design <- model.matrix(~0 + group)


colnames(design) <- c(levels(group))

contr<-  makeContrasts(KnTp - KnTn,
                       KpTn - KnTn,
                       KpTp - KnTn,
                       KpTn - KnTp,
                       KpTp - KnTp,
                       KpTp - KpTn,
                       KnTn - WT,
                       KnTp - WT,
                       KpTn - WT,
                       KpTp - WT,
                       levels=design)

fit <- lmFit(mval, design)

fit <- contrasts.fit(fit, contr[,1:10])
fit2 <- eBayes(fit)

#, lfc=log2(1.5)

dt=decideTests(fit2  , lfc=0)

summary(dt)

lap=0

DP =data.frame()
for(i in seq(1+lap,10)){
  s <- topTable(fit2, num=Inf, coef=i, sort.by="p")
  s$probe=rownames(s)
  s$contrast= colnames(contr)[i]
  DP=rbind(DP,s)
}


DP=DP[DP$adj.P.Val <=0.05 ,]
rownames(DP)=1:nrow(DP)

save(DP, file="DMethDE022724.RData")

#save(DP, file="DMethDE011124.RData")
#save(DP, file="DMethDE122323.RData")
#save(DP, file="DMethDE131.RData")
#save(DP, file="DMethDE129.RData")
#save(DP, file="DMethDE2224.RData")

#save(DP, file="DMethDE020624.RData")


#save( probefeatures, file="DFeatureNewIdea.RData")


# extracting significantly methylated probes
deff.meth = topTable(fit2, coef=1, sort.by="p",number = nrow(mval), adjust.method = "BY")


# Visualization
# plot the top 10 most significantly deferentially methylated CpGs 
par(mfrow=c(2,5))
sapply(rownames(deff.meth)[1:10], function(cpg){
  plotCpg(meth, cpg=cpg, pheno=Clinical$grp, ylab = "Beta values")
})



# making a volcano plot
#making dataset

10#KnTn      "#F8766D"   green  
11#KnTnHB  "#A3A500" brown
12#KnTp     "#00BF7D"   red
13#KpTn    "#00B0F6"  blue
14#KpTp  "#E76BF3"  pink

contrscol=7

dat <- data.frame(foldchange = fit2[["coefficients"]][,contrscol], 
               logPvalue =  -log10(fit2[["p.value"]][,contrscol]),
               adj.Pvalue = p.adjust(fit2[["p.value"]][,contrscol],method = "BH"))


dat$threshold <- as.factor(abs(dat$foldchange) < 1)
dat$direction = sign(dat$foldchange)
hypo = nrow(dat[dat$direction == -1 & dat$adj.Pvalue <= 0.05,])
hyper = nrow(dat[dat$direction == 1 & dat$adj.Pvalue <= 0.05,])

file <- tempfile()
ggsave(file, device = "png")
unlink(file)
#Visualization
#library(ggplot2)
cols <- c("TRUE" = "grey", "FALSE" =  "#A3A500" )
pllt= ggplot(data=dat, aes(x=foldchange, y = logPvalue, color=threshold)) +
  geom_point(alpha=.6, size=1.2) +
  scale_colour_manual(values = cols) +
  geom_vline(xintercept =1, colour="#990000", linetype="dashed") + 
  geom_vline(xintercept = - 1, colour="#990000", linetype="dashed") +
  geom_hline(yintercept = 1.5, colour="#990000", linetype="dashed") +
  geom_text(x = -4, y = 20, label = paste0(" hypomethylated: ", hypo), color = "#3D3D3D", size = 4) +
  geom_text(x = 4, y = 20, label = paste0(" hypermethylated: ", hyper), color = "#3D3D3D", size = 4) +
  theme(legend.position="none") +
  xlab(paste0(colnames(contr)[contrscol],"  Fold Change")) +
  ylab("-log10 p value") +
  theme_bw() +
  theme(legend.position = "none") + ggsave("methplot7.png")


#KnTn - WT KnTp - WT KpTn - WT KpTp - WT
# setting some annotation

myAnnotation <- cpg.annotate(object = mval, datatype = "array", 
                             what = "M", 
                             analysis.type = "differential", 
                             design = design, 
                             contrasts = TRUE,
                             cont.matrix = contr,
                             arraytype = "450K",
                             coef="KpTp - WT",
                             fdr = 0.05)
# 
str(myAnnotation)
#coef = "paper_Histologic.gradeHigh Grade"
# DMR analysis
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
#results.ranges
#----------------------------------------------------------------------

# visualization
dmr.tableKpTp <- data.frame(results.ranges)

# setting up variable for grouping and color

# set up the grouping variables and colours
pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(Clinical$grp))]
names(groups) <- levels(factor(Clinical$grp))



#setting up the genomic region 
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 6363
# coordinates are stored under results.ranges[dmrIndex]

chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))

# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))



# defining CpG islands track
# download cpgislands for chromosome number 6 from ucsc
chr6.cpg <- read.csv("chr6-cpg.csv")

islandData <- GRanges(seqnames=Rle(chr6.cpg[,1]), 
                      ranges=IRanges(start=chr6.cpg[,2],
                                     end=chr6.cpg[,3]),
                      strand=Rle(strand(rep("*",nrow(chr6.cpg)))))

# DNAseI hypersensitive sites track
#downloaded from ucsc
chr6.dnase <- read.csv("chr6-dnase.csv")

dnaseData <- GRanges(seqnames=chr6.dnase[,1],
                     ranges=IRanges(start=chr6.dnase[,2], end=chr6.dnase[,3]),
                     strand=Rle(rep("*",nrow(chr6.dnase))),
                     data=chr6.dnase[,5])

#Setting up the ideogram, genome and RefSeq tracks 

iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name=paste0(chrom))
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

#Ensure that the methylation data is ordered by chromosome and base position.

ann450kOrd <- ann450k[order(ann450k$chr,ann450k$pos),]
bvalOrd <- meth[match(ann450kOrd$Name,rownames(meth)),]

#Create the data tracks:
# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(ann450kOrd$chr),
                   ranges=IRanges(start=ann450kOrd$pos, end=ann450kOrd$pos),
                   strand=Rle(rep("*",nrow(ann450kOrd))),
                   betas=bvalOrd)

# methylation data track
methTrack <- DataTrack(range=cpgData, 
                       groups=Clinical$grp, # change this if your groups are diffrent
                       genome = gen,
                       chromosome=chrom,
                       ylim=c(-0.05,1.05),
                       col=pal,
                       type=c("a","p"), 
                       name="DNA Meth.\n(beta value)",
                       background.panel="white", 
                       legend=TRUE, 
                       cex.title=0.8,
                       cex.axis=0.8, 
                       cex.legend=0.8)

# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")


# Set up the track list and indicate the relative sizes of the different tracks. 
# Finally, draw the plot using the plotTracks function
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack,
               rTrack)
sizes <- c(2,2,5,2,2,2,3) # set up the relative sizes of the tracks

tiff( filename = "dmrtkptp.tiff", width = 15, height = 10, units = "in", res = 400)
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
dev.off()

#----------------------------------


#__________________________Differential variability_________________#
fitvar <- varFit(mval, design = design)

ff=topVar(fitvar,coef=6)

fit <- contrasts.fit(fitvar , contr[,1:10])
fit2 <- eBayes(fit)

# Summary of differential variability
dtv=(decideTests(fit2))

summary(dtv)


lap=0
DVP =data.frame()
for(i in seq(1+lap,10)){
  s <- topTable(fit2, num=Inf, coef=i, sort.by="p")
  s$probe=rownames(s)
  s$contrast= colnames(contr)[i]
  DVP=rbind(DVP,s)
}

DVP=DVP[DVP$adj.P.Val <=0.05 ,]
rownames(DVP)=1:nrow(DVP)


#save(DVP, file="DVP.Rdata")


dtt= as.data.frame(dtv@.Data)

top_gene=colnames(dtt)

mset=list()
dtt$id=rownames(dtt)
snv.set=dtt

for(i in 1:4){
  vv=colnames(snv.set)[i]
  lst=list(unique(as.vector(snv.set[snv.set[[i]]==1,length(top_gene)+1])))
  lst2=list(unique(as.vector(snv.set[snv.set[[i]]==-1,length(top_gene)+1])))
  
  names(lst)=paste0("Hyper.",vv)
  names(lst2)=paste0("Hypo.",vv)
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




setordrb = c("Hyper.KnTp - KnTn", "Hypo.KnTp - KnTn",  "Hyper.KpTn - KnTn",
"Hypo.KpTn - KnTn",  "Hyper.KpTp - KnTn", "Hypo.KpTp - KnTn",
"Hyper.KpTn - KnTp", "Hypo.KpTn - KnTp" )

colcode= c( "#F8766D","#00BF7D",  "#00B0F6", "#E76BF3")

UpSet(reg, set_order = setordrb[c(1,3,5,7,2,4,6,8)], comb_col = colcode[comb_degree(reg)], bg_col = c("white", "#E6E6E6"), top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
      right_annotation = upset_right_annotation(reg, add_numbers = TRUE, 
                                                axis_param = list(side = "top"),
                                                annotation_name_side = "top", 
                                                gp = gpar(fill = "#A3A500" )))



#------------------Integrated analysis-----------------------------------

DMPs=DP[DP$contrast =="KnTnHB - KnTn",]
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$probe[DMPs$adj.P.Val<0.05]
all <- rownames(mval)

# Run enrichment - Can take a bit of time...
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all)
#gstk <- gometh(sig.cpg=sigCpGs, all.cpg=all, collection = "KEGG")
# Top 10 GO categories
term4=topGSA(gst, number=100)




DMPs=DVP[DVP$contrast =="KpTp - KnTn",]
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$probe[DMPs$adj.P.Val<0.05]
all <- rownames(mval)

# Run enrichment - Can take a bit of time...
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all)
#gstk <- gometh(sig.cpg=sigCpGs, all.cpg=all, collection = "KEGG")
# Top 10 GO categories
topGSA(gst, number=10)

