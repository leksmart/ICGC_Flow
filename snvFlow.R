#/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/validation
getwd()
root="/mmfs1/home/david.adeleke/WorkFlow"
#root <- "~"
main <- "Chapter_Two"

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main)))}

setwd(file.path(root, main))

library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(RTCGA.clinical)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(ggplot2)
library(DT)
library(tidyr)
library(reshape2)
library(plyr)
library(plotly)
library(tibble)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(GGally)
library(kableExtra)
library(survival)
library(survminer)
library(ggplot2)
library(grid)
library(gridExtra)
library(rmdformats)


setwd(  file.path(root, "ChapterTwo"))

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

##Start

initia= function(snv,i){
  
  excl_var= c("synonymous_variant", "Silent")
  proj.scope = c("PAAD-US", "PACA-CA", "PACA-AU")
  
  clinical=clinical %>% filter(project_code %in% proj.scope[c(i)])
  snv<- snv %>% filter(project_code %in% proj.scope[i]) %>%
    filter(!(Variant_Classification %in% excl_var))
  snv= as.data.frame(snv)
  
  cn.id = intersect(clinical$icgc_donor_id,snv$icgc_donor_id)
  
  clinical=clinical[clinical$icgc_donor_id %in% cn.id,]
  snv=snv[snv$icgc_donor_id %in% cn.id,]
  
  snv$mutsta=ifelse(snv$Tumor_Seq_Allele2 != snv$Reference_Allele,1,0)
  snv.cols= snv %>% dplyr::select(c(1,22, 50)) %>% unique()
  
  names(snv.cols)= c("gene", "id", "mut")
  
  snv.mat= snv.cols %>% pivot_wider(names_from = id, values_from = mut, values_fill = 0)
  
  mat= as.matrix(snv.mat[,-1])
  rn= snv.mat %>% pull(gene)
  rownames(mat)= rn
  rm(rn)
  mat <<- mat
}


exPlt= function(file_name="plt.png"){
  # Specify the folder path where you want to save the plot
  folder_path <- file.path(getwd(), "images")
  # Specify the fil_ename for the plot
  file_name <- file_name
  #deparse(quote(pp))
  # Create the complete file path
  file_path <- file.path(folder_path, file_name)
  # Save the plot as a PNG file
  ggsave(file_path, pp, width = 6, height = 4, dpi = 300)
  cat("Plot saved at", file_path, "\n")
}

wd= getwd()

imgwd <- file.path(wd, "images")

datawd <- file.path(getwd(),"output")

# usidtable=tumoridtable
# save(usidtable, file="usidtable.RData")
# 
# asd=load("usidtable.RData")
# normid=c("DO32769N", "DO32751N", "DO50272N")
# 
# dd=data.frame(icgc_donor_id=normid, submitted_sample_id=c("TCGA-H6-A45N-11A-11D-A26I-08",
#                                                           "TCGA-H6-8124-11A-11D-2396-08",
#                                                           "TCGA-YB-A89D-11A-12D-A36O-08") , submitted_matched_sample_id=NA, grp="WT")
# 
# tumoridtable=rbind(tumoridtable,dd)

# "TCGA-YB-A89D-11A-11R-A36B-13"
# "TCGA-H6-A45N-11A-12R-A26Y-13" 
# "TCGA-H6-8124-11A-01R-2401-13" 
# "TCGA-HV-A5A3-11A-11R-A26Y-13"

coho=3

snv1 <- "dataset/PAAD-US_ICGC_snv.rds"
snv2 <-  "dataset/PACA-CA_ICGC_snv.rds"
snv3 <-  "dataset/PACA-AU_ICGC_snv.rds"


snv=c(snv1,snv2, snv3)
snv=readRDS(file.path(root,snv[coho]))
snv= snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")


#load("usidtable.RData")
tumoridtable=unique(snv[,c(22,26,27)])

# Return here after runnining Kntn, Kntp vectors


#snv=readRDS(file.path(root,snvall[coho]))
#snv= snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")

clinical1 <- "dataset/PAAD-US_ICGC_clin.rds"
clinical2 <-  "dataset/PACA-CA_ICGC_clin.rds"
clinical3 <- "dataset/PACA-AU_ICGC_clin.rds"


clinical=c(clinical1,clinical2,clinical3)
clinical=readRDS(file.path(root,clinical[coho]))

snv=initia(snv=snv, i=coho)

snv=as.data.frame(snv)

KnTn = colnames(snv[, which(snv["KRAS", ] == 0 & snv["TP53", ] == 0 )]);length(KnTn)
KnTp = colnames(snv[, which(snv["KRAS", ] == 0 & snv["TP53", ] == 1 )]);length(KnTp)
KpTp = colnames(snv[, which(snv["KRAS", ] == 1 & snv["TP53", ] == 1 )]);length(KpTp)
KpTn = colnames(snv[, which(snv["KRAS", ] == 1 & snv["TP53", ] == 0 )]);length(KpTn)

tumoridtable$grp=NA
tumoridtable$grp= ifelse(tumoridtable$icgc_donor_id %in% KnTn,"KnTn",tumoridtable$grp)
tumoridtable$grp= ifelse(tumoridtable$icgc_donor_id %in% KnTp,"KnTp",tumoridtable$grp)
tumoridtable$grp= ifelse(tumoridtable$icgc_donor_id %in% KpTp,"KpTp",tumoridtable$grp)
tumoridtable$grp=ifelse(tumoridtable$icgc_donor_id %in% KpTn,"KpTn",tumoridtable$grp)
tumoridtable=na.omit(tumoridtable)

#clinicals

clin=clinical[,c(1,7,22,23)]
colnames(clin)= c("id","chemo","status", "time")  
clin$grp=NA
clin$grp= ifelse(clin$id %in% KnTn,"KnTn",clin$grp)
clin$grp= ifelse(clin$id %in% KnTp,"KnTp",clin$grp)
clin$grp= ifelse(clin$id %in% KpTp,"KpTp",clin$grp)
clin$grp= ifelse(clin$id %in% KpTn,"KpTn",clin$grp)

table(clin$grp)
clin=clin[!is.na(clin$grp),]
ref="KnTn"
clin$k.status= ifelse(substring(clin$grp,1,2) == "Kp","KRAS", clin$grp)
clin$k.status=as.factor(clin$k.status)
clin$grp=as.factor(clin$grp)
clin <- within(clin, grp <- relevel(grp, ref = ref))
clin3 <- within(clin, k.status <- relevel(k.status, ref = ref))

clin1$cohort= "USA"
clin2$cohort= "CAN"
clin3$cohort= "AUS"
clin= rbind(clin1,clin2)
clin= rbind(clin, clin3)
sfit <- survfit(Surv(time, status)~grp, data=clin[clin$time >20,])

#[clin$grp %in% c("Non_KRAS_Lo_TMB", "Non_KRAS_Hi_TMB")

pp=ggsurvplot(
  sfit,   
  # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,  # show confidence intervals for 
  palette= c( "#47BC82","#FAA928","#B95D4B", "#A18CBA"),
  # point estimaes of survival curves.
  xlim = c(0,2500),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 250,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),
  # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE,
  risk.table.height=0.4,
  legend.title="US+AU+CA") # show bars instead of names in text annotations

pp


#-------------------------------------------------------------------kNtNHB inclusion------

load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

clinical=readRDS(file.path(root,"dataset/PAAD-US_ICGC_clin.rds"))

clinical=merge(clinical, usidtable, by="icgc_donor_id")

clinical=clinical[,c(1,22,23,28)]


clin=clinical
colnames(clin)[2]= "status"
table(clin$grp)
clin=clin[!is.na(clin$grp),]
ref="KnTn"
clin$k.status= ifelse(substring(clin$grp,1,2) == "Kp","KRAS", clin$grp)
clin$k.status=as.factor(clin$k.status)
clin$grp=as.character(clin$grp)
clin$grp= ifelse(clin$grp == "KnTnHB", "KnTn", clin$grp)
clin$grp=as.factor(clin$grp)
clin <- within(clin, grp <- relevel(grp, ref = ref))
clin <- within(clin, k.status <- relevel(k.status, ref = ref))
sfit <- survfit(Surv(time, status)~grp, data=clin[clin$time >5,])

#[clin$grp %in% c("Non_KRAS_Lo_TMB", "Non_KRAS_Hi_TMB")
#"#00B0F6",
mycol= c(  "#F8766D",  "#00BF7D" , "#00B0F6" , "#E76BF3" )

#c( "#F8766D", "#A3A500","#00BF7D", "#E76BF3")
pp=suppressMessages(ggsurvplot(
  sfit,   
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = F,  # show confidence intervals for 
  palette= mycol,
  # point estimaes of survival curves.
  xlim = c(0,2500),        # present narrower X axis, but not affect
  # survival estimates.
  break.time.by = 250,    # break X axis in time intervals by 500.
  ggtheme = theme_RTCGA(),
  # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE,
  risk.table.height=0.35,
  legend.title="MBS")) # show bars instead of names in text annotations

pp


#--------------------------------------------------------------------------------------------#

clin <- within(clin, grp <- relevel(grp, ref = ref))
tt=coxph(Surv(time, status) ~ grp, data = clin[clin$time >5,]) 
ttt = tt %>%  gtsummary::tbl_regression(exp = TRUE) 

ttt

suppressWarnings(ggforest(tt))

#saveRDS(clin, file="allclinical.rds")



################################################################

top_gene = c("KRAS", "TP53", "TTN", "MUC16","FLG", "OBSCN", "RYR1", "SYNE1", "ARID1A",
             "RYR2", "LRP1B")

snv.set=mat[top_gene,]

mset=list()

snv.set=as.data.frame(snv.set)
snv.set$id=rownames(snv.set)

for(i in 1:(length(top_gene)-1)){
  vv=colnames(snv.set)[i]
  lst=list(unique(as.vector(snv.set[snv.set[[i]]==1,length(top_gene)+1])))
  names(lst)=vv
  mset <- c(mset,  lst)
}

reg=mset[c(1:length(top_gene)-1)]
reg=mset

reg= ComplexHeatmap::make_comb_mat(reg, mode="distinct")



UpSet(reg, top_annotation = upset_top_annotation(reg, add_numbers=TRUE),
      right_annotation = upset_right_annotation(reg, add_number=TRUE),
      comb_col = c("red", "blue", "orange", "black", "green")[comb_degree(reg)])





##############################
snv=c(snv1,snv2, snv3)
snv=readRDS(file.path(root,snv[coho]))
snv= snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")
clinical=readRDS(file.path(root,"dataset/PAAD-US_ICGC_clin.rds"))

snv=initia(snv=snv, i=coho)

snv=as.data.frame(snv)
top_gene=c("KRAS", "TP53")

mat=snv[top_gene,]

mat=t(mat)
mat=as.data.frame(mat)
mat$None=ifelse((mat$KRAS + mat$TP53)==0,1,0)


top_gene=c("KRAS", "TP53", "None")

mset=list()
mat$id=rownames(mat)
snv.set=mat

for(i in 1:length(top_gene)){
  vv=colnames(snv.set)[i]
  lst=list(unique(as.vector(snv.set[snv.set[[i]]==1,length(top_gene)+1])))
  names(lst)=vv
  mset <- c(mset,  lst)
}

reg=mset[c(1:length(top_gene))]
reg=mset
reg= ComplexHeatmap::make_comb_mat(reg, mode="distinct")


m=reg 

ppt= UpSet(m, comb_col = c("#E76BF3", "#00B0F6","#00BF7D","#F8766D"), top_annotation = HeatmapAnnotation(show_annotation_name = T,
                                            annotation_label = c("Subtype",""),
                                            col = list(Subtype = c("KpTp" = "#E76BF3","KpTn" = "#00B0F6", 
                                                                   "KnTp" = "#00BF7D", "KnTn"="#F8766D")),
                                            Subtype = c("KpTp", "KpTn","KnTp", "KnTn"),
  "distinct\nsize" = anno_barplot(comb_size(m), add_numbers = T,
                                      border = FALSE, 
                                      gp = gpar(fill = c("white")), 
                                      height = unit(2, "cm")),  annotation_name_side = "right",
  annotation_name_rot = 0),right_annotation = upset_right_annotation(reg, add_numbers = TRUE))





UpSet(m,   comb_col = c("red", "blue", "orange", "black", "green"),
      bg_col = c("white", "#E6E6E6"), top_annotation = upset_top_annotation(reg, add_numbers = TRUE),
      right_annotation = upset_right_annotation(reg, add_numbers = FALSE))


###------------------------------Mutation Burden Heat Map-------------------------------


load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

clinical=readRDS(file.path(root,"dataset/PAAD-US_ICGC_clin.rds"))

clinical=merge(clinical, usidtable, by="icgc_donor_id")

clin=clinical[,c(1,22,23,28)]

root="/mmfs1/home/david.adeleke/WorkFlow"
snv1 <- "dataset/PAAD-US_ICGC_snv.rds"
snv2 <-  "dataset/PACA-CA_ICGC_snv.rds"
snv3 <-  "dataset/PACA-AU_ICGC_snv.rds"
coho=1
snv=c(snv1)
snv=readRDS(file.path(root,snv[coho]))
snv= snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")

snv=initia(snv=snv, i=coho)

snv=as.data.frame(snv)

snv$grp= rowSums(mat, na.rm = TRUE)

snv= snv[snv$grp >=8,-177]

mat= as.matrix(snv)

# Calculate row sums
rowSumsVector <- rowSums(mat)

# Sort indices by decreasing row sum
sortedIndices <- order(-rowSumsVector)

# Use the sorted indices to reorder the rows of the matrix
mat <- mat[sortedIndices, , drop = FALSE]



load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

usidtable= usidtable[usidtable$icgc_donor_id %in% colnames(mat),]

usidtable$grp=ifelse(usidtable$grp == "KnTnHB", "KnTn",  usidtable$grp)

KnTn=usidtable[usidtable$grp=="KnTn",1]

KnTnHB=usidtable[usidtable$grp=="KnTnHB",1]

KnTp=usidtable[usidtable$grp=="KnTp",1]

KpTn= usidtable[usidtable$grp=="KpTn",1]
KpTp = usidtable[usidtable$grp=="KpTp",1]

mat= mat[,c( KnTn, KnTp,  KpTn, KpTp)]

split = c( rep(1, each = length(KnTn)),
                    rep(2, each = length(KnTp)),
                       rep(3, each = length(KpTn)),
                           #rep(4, each = length(KnTnHB)),
                               rep(5, each = length(KpTp)))


#"#A3A500",

mycol= c(  "#F8766D",  "#00BF7D" , "#00B0F6" , "#E76BF3" )


#KnTn "#47BC82"
#KnTnHB "#936FF9"
#KnTp   #FAA928 
#KpTn  #B95D4B
#KpTp  #A18CBA

ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE),
  foo = anno_block(gp = gpar(fill = mycol), 
  labels = c("KnTn",  "KnTp", "KpTn",  "KpTp")))

colnames(mat)= NULL

Heatmap(mat, name = "Mutation", cluster_rows = F, cluster_columns = F,
        col=c("#EEEEEE", "#117C81") , column_split = split, top_annotation = ha, 
        column_title = NULL)

#-------------------------------------------------------------------------------

load("/mmfs1/home/david.adeleke/WorkFlow/ChapterTwo/usidtable.RData")

clinical=readRDS(file.path(root,"dataset/PAAD-US_ICGC_clin.rds"))

clinical=merge(clinical, usidtable, by="icgc_donor_id")

clin=clinical[,c(1,22,23,28)]

root="/mmfs1/home/david.adeleke/WorkFlow"
snv1 <- "dataset/PAAD-US_ICGC_snv.rds"
snv2 <-  "dataset/PACA-CA_ICGC_snv.rds"
snv3 <-  "dataset/PACA-AU_ICGC_snv.rds"
coho=1
snv=c(snv1)
snv=readRDS(file.path(root,snv[coho]))
snv= snv %>% dplyr::filter(Hugo_Symbol != "UnknownGene")

snv=initia(snv=snv, i=coho)

snv=as.data.frame(snv)
snv=t(snv)

mutb= data.frame(Total=rowSums(snv), id= rownames(snv))

mutb$MBS= ifelse(mutb$id %in% KnTn, "KnTn",
                 ifelse(mutb$id %in% KnTp, "KnTp",
                        ifelse(mutb$id %in% KpTn, "KpTn",
                               ifelse(mutb$id %in% KpTp, "KpTp",
                                      ifelse(mutb$id %in% KnTnHB, "KnTnHB","UN")))))
                        
outli="DO46657"

mutb=mutb[mutb$id != outli,]

# 
# # Visualize: Specify the comparisons you want
# my_comparisons <- list( c("KnTn", "KnTnHB"), c("KnTn", "KnTp"), c("KnTn", "KpTn"), c("KnTn", "KpTp") )
# ggboxplot(mutb, x = "MBS", y = "Total",
#           color = "MBS", palette = "jco")+
#   stat_compare_means(method = "anova", label.y = 120)+      # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test",
#                      ref.group = "KnTn")  



ggplot(mutb, aes(x = MBS, y =Total, fill = MBS)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  scale_fill_manual(values = mycol) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Total Mutation Burden  "))+
  stat_compare_means(method = "anova", label.y = 180)+      # Add global p-value
  stat_compare_means(label = "p.signif", method = "t.test",
                     ref.group = "KnTn")




