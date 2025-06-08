library(mvlearnR)

setwd("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/")

#load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/validation/auvalidationset.RData")

contrst= c("KnTn", "KnTnHB", "KnTp", "KpTn", "KpTp" )
Clinical=readRDS("clinical_New.rds")
rownames(Clinical)=Clinical$idd
j=contrst[1]
k=contrst[3]
z=contrst[4]
w=contrst[5]

dustin=k
jwd=z


#contrs= paste0(dustin, " - ", jwd)
#contrs= paste0(jwd, " - ", dustin)

contrs= "KpTp - KnTn"

contrx= paste0(nat, " - ", dustin)

contrx=""
nat=""
#"KnTp - KnTn"
#k cc dustin,
Clinical$grp=as.character(Clinical$grp)
Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn",  Clinical$grp )

Clinical$class= ifelse(Clinical$grp %in% c(dustin),0,
                       ifelse(Clinical$grp %in% c(jwd),1,
                              ifelse(Clinical$grp %in% c(""),3,3)))



Clinical= Clinical[ Clinical$class < 3,]


load( "validation/adj.counts115.RData")

load( "validation/DMGsraw115.RData")

adj.counts<- edgeR::cpm(adj.counts, log=TRUE)

RNA=t(adj.counts)
rm(adj.counts)



#arcsin value from bvalue
load("DNAmet011124.RData")

#DE from mval from arcsin of bval
load("DMethDE011124.RData")
#meth=mval

Meth=t(mval)
rm(mval)

de.genes=DE
rm(DE)
de.genes=de.genes[de.genes$contrast %in% c(contrs,contrx) & 
                    abs(de.genes$logFC) >= 1 & de.genes$adj.P.Val < 0.05 ,7]

#DP=DVP
DP=DP[DP$contrast %in% c(contrs,contrx)  & 
        abs(DP$logFC) >= 1.5 & DP$adj.P.Val < 0.05 ,7]


#save(promoter.probe, file="promoter.probe.RData")


load("promoter.probe.RData")


DP=intersect(DP, promoter.probe)

ids=intersect(Clinical$id, rownames(Meth))
ids=intersect(ids, rownames(RNA))
Clinical=Clinical[Clinical$id %in% ids,]
Meth=Meth[Clinical$id,DP]

RNA=RNA[Clinical$id,intersect(colnames(RNA),de.genes )]

X=list( RNA, Meth)


Y=Clinical$class
Y=Y+1

omics= filter.unsupervised(X, method = "variance",
                           pct.keep = 85,
                           center = FALSE,
                           scale = FALSE,
                           standardize = FALSE,
                           log2TransForm = FALSE,
                           Xtest = NULL)

#rm(Meth, RNA)

Xdata1 <- omics$X[[1]]
Xdata2 <- omics$X[[2]] 

#Xdata3 <- t(adj.mi.counts[,rownames(Xdata2)])

dim(Xdata1)

dim(Xdata2)

#dim(Xdata3)

Xdata=list(Xdata1, Xdata2)

########################################################################################

fit.cvsida <- cvSIDA(Xdata, Y,
                     withCov = F,
                     Xtestdata = Xdata,
                     Ytest = Y, plotIt = F,
                     isParallel = TRUE,
                     ncores = 32,
                     AssignClassMethod="Joint")

CorrelationPlots(Xdata, Ytest=Y, fit.cvsida$hatalpha)

needt=networkPlot(fit.cvsida,cutoff=0.6)

mymat = networkplotinner(fit.cvsida)

myComb = mymat$ViewCombinations
nComb = dim(mymat$ViewCombinations)[2]


networkPlot(fit.cvsida,cutoff=0.5)

DiscriminantPlots(Xdata, Y, fit.cvsida$hatalpha)

#rsd=data.frame()
rsd1= VarImportancePlot(fit.cvsida)

rsd1=rsd1[[1]][[1]]

rsd1=rsd1[rsd1$Loading !=0,]

rsd1$grp= paste0(dustin," - ", jwd)

rsd2= VarImportancePlot(fit.cvsida)
rsd2=rsd2[[1]][[2]]

rsd2=rsd2[rsd2$Loading != 0,]

rsd2$grp=  paste0(dustin," - ", jwd)



rsd=rbind(rsd,rsd1)

rsd=rbind(rsd,rsd2)


#save(rsd, file="validation/rsd116.RData")
write.csv(rsd, file="validation/rsd.csv")

impvar= data.frame()
impvar=rbind(impvar,rsd1)


"joint"
tdfs=data.frame(fit.cvsida$PredictedClass, Y)
table(tdfs$fit.cvsida.PredictedClass, tdfs$Y)


#-------------------------------PLOT-----------------------------------

Clinical=readRDS("clinical_New.rds")


nomclin= Clinical[1:3,]
normid= c("DO32751N", "DO32769N", "DO50272N")

nomclin$id= normid
rownames(nomclin)= normid
nomclin$grp= "WType"
Clinical=rbind(Clinical, nomclin)


Clinical$grp=as.character(Clinical$grp)
Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn",  Clinical$grp )



load( "validation/adj.counts115.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=RNAlogCPM[,unique(impvar$`Variable Name`)]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,unique(rsd2$`Variable Name`)]
Dmem.traits.phen=merge(meth, Clinical, by=0)


library(ggplot2)
library("ggpubr")


#if(!require(devtools)) install.packages("devtools")

#devtools::install_github("kassambara/ggpubr")


mem.traits.phen$grp <- factor(mem.traits.phen$grp, levels=c("WType","KnTn", "KnTnHB",   "KnTp",   "KpTn",   "KpTp" ))

Dmem.traits.phen$grp <- factor(Dmem.traits.phen$grp, levels=c("WType","KnTn", "KnTnHB",   "KnTp",   "KpTn",   "KpTp" ))


#compare_means(DHRS3 ~ grp, data = mem.traits.phen)

iii=1


ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Expression of  ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)



#-------------------------SURVIVAR-------------------------------------------------------------------#
library(survminer)
library(survival)
library(utils)

covariates= unique(impvar$`Variable Name`)
Dcovariates= unique(rsd2$`Variable Name`)
#covariates<- paste0("ME", covariates)

mem.traits.phenbk=mem.traits.phen

mem.traits.phen=mem.traits.phen[mem.traits.phen$grp %in% c(dustin, jwd),]

Dmem.traits.phen=Dmem.traits.phen[Dmem.traits.phen$grp %in% c(dustin, jwd),]

# ---------------------------LASSO REGRESSION------------------------------------------#
Dmem.traits.phen=Dmem.traits.phen[Dmem.traits.phen$time > 5,]

endl= ncol(Dmem.traits.phen)-5
library(glmnet)
x <- as.matrix(Dmem.traits.phen[, 2:endl])
y <- Surv(Dmem.traits.phen$time, Dmem.traits.phen$status)
lasso_model <- cv.glmnet(x, y, family = "cox")


# Plot the cross-validated mean squared error (optional)
plot(lasso_model)

# Identify selected genes based on the optimal lambda value
optimal_lambda <- lasso_model$lambda.min
selected_genes <- coef(lasso_model, s = optimal_lambda)
selected_genes <- selected_genes[selected_genes[, 1] != 0, ]

# Print the selected genes
print(selected_genes)

goi="cg27583010"



# Define the list of variables
goi <- names(selected_genes)
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = Dmem.traits.phen)

summary(modelph)


#-------------------

library(glmnet)
mem.traits.phen=mem.traits.phen[mem.traits.phen$time > 5,]

endl= ncol(mem.traits.phen)-5

x <- as.matrix(mem.traits.phen[, 2:endl])
y <- Surv(mem.traits.phen$time, mem.traits.phen$status)
lasso_model <- cv.glmnet(x, y, family = "cox")


# Plot the cross-validated mean squared error (optional)
plot(lasso_model)

# Identify selected genes based on the optimal lambda value
optimal_lambda <- lasso_model$lambda.min
selected_genes <- coef(lasso_model, s = optimal_lambda)
selected_genes <- selected_genes[selected_genes[, 1] != 0, ]

# Print the selected genes
print(selected_genes)



# Define the list of variables
goi <- names(selected_genes)
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
mem.traits.phendk=  mem.traits.phenbk[ mem.traits.phenbk$grp %in% c(dustin, jwd),]

modelph <- coxph(formula_str, data = mem.traits.phendk)

summary(modelph)



#-------------------complex heat-----------------
matdf= mem.traits.phen[,2:endl]
matdf=as.matrix(matdf)
matdf=t(matdf)
base_mean = rowMeans(matdf)
mat_scaled = t(apply(matdf, 1, scale))

type =mem.traits.phen$grp

ha = HeatmapAnnotation(type = type, annotation_name_side = "left")


ht_list = Heatmap(mat_scaled, name = "expression", row_km = 5, 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) +
  Heatmap(base_mean, name = "base mean", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), 
                                                                    height = unit(2, "cm"))),
          width = unit(15, "mm"))


library(ComplexHeatmap)
library(circlize)

expr = readRDS(system.file(package = "ComplexHeatmap", "extdata", "gene_expression.rds"))
mat = as.matrix(expr[, grep("cell", colnames(expr))])
base_mean = rowMeans(mat)
mat_scaled = t(apply(mat, 1, scale))

type = gsub("s\\d+_", "", colnames(mat))
ha = HeatmapAnnotation(type = type, annotation_name_side = "left")

ht_list = Heatmap(mat_scaled, name = "expression", row_km = 5, 
                  col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
                  top_annotation = ha, 
                  show_column_names = FALSE, row_title = NULL, show_row_dend = FALSE) +
  Heatmap(base_mean, name = "base mean", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6), 
                                                                    height = unit(2, "cm"))),
          width = unit(15, "mm")) +
  rowAnnotation(length = anno_points(expr$length, pch = 16, size = unit(1, "mm"), 
                                     axis_param = list(at = c(0, 2e5, 4e5, 6e5), 
                                                       labels = c("0kb", "200kb", "400kb", "600kb")),
                                     width = unit(2, "cm"))) +
  Heatmap(expr$type, name = "gene type", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
          width = unit(15, "mm"))

ht_list = rowAnnotation(block = anno_block(gp = gpar(fill = 2:6, col = NA)), 
                        width = unit(2, "mm")) + ht_list

draw(ht_list)




KntnKpTp.DNA= c(cg11571304, cg18610930)
KntnKpTp.RNA= c(SLC9A3,BTBD6,  FHAD1, TLX1, SLC39A10, CDK6)

KntnKpTn.RNA=c(SLC9A3, GPR126, HIST1H1E, NEK10)

KntnKpTn.DNA =  "cg27583010"

KntnKnTp.RNA= FOXL1

#table(rsd$grp)

# KnTp - KnTn KnTp - KpTn KpTn - KnTn KpTp - KnTn

KpTp - KnTn
rsd1= rsd[rsd$View == 1, ]
rsd2= rsd[rsd$View == 2, ]

length(unique(rsd2$`Variable Name`))
length((rsd2$`Variable Name`))

rsd1= rsd[rsd$View == 1 & rsd$grp == "KpTp - KnTn", 2]

rsd3= rsd[rsd$View == 1 & rsd$grp == "KnTp - KnTn", 2]

rsd2= rsd[rsd$View == 1 & rsd$grp == "KpTn - KnTn", 2]

rsd2= rsd[rsd$View == 2, ]

#-----------------------------------------------------------------
NOTCH1 regulation of endothelial cell calcification

ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Expression of  ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)

