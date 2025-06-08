library(mvlearnR)
library(ggplot2)
library("ggpubr")
library(survminer)
library(survival)
library(utils)

setwd("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/")

#goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" ,"cg06952671" ,"cg10661002")
#"DO218962" "DO32763"  "DO46622"  "DO51677" 
contrst= c("KnTn", "KnTnHB", "KnTp", "KpTn", "KpTp" )
Clinical=readRDS("clinical_New.rds")
rownames(Clinical)=Clinical$idd
j=contrst[1]
cc=contrst[2]
k=contrst[3]
z=contrst[4]
w=contrst[5]

dustin=w
jwd=j
nat= "KpTn"

contrs= paste0(dustin, " - ", jwd)
#contrs= paste0(jwd, " - ", dustin)
contrx= paste0(dustin, " - ", jwd)

contrx=""

#"KnTp - KnTn"
#k cc dustin,
#Clinical$grp=as.character(Clinical$grp)
#Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn",  Clinical$grp )

Clinical$class= ifelse(Clinical$grp %in% c(jwd),0,
                       ifelse(Clinical$grp %in% c(dustin),1,
                              ifelse(Clinical$grp %in% c(nat),3,3)))

Clinical= Clinical[ Clinical$class < 3,]


#load( "validation/adj.counts1221.RData")
#load( "validation/adj.counts022724.RData")
#load( "validation/DMGsraw022724.RData")


load(file=file.path("/mmfs1/home/david.adeleke/Thesis/dataset","PAAD-US_adj.counts.RData" ))

load(file=file.path("/mmfs1/home/david.adeleke/Thesis/dataset","PAAD-US_DMGs.RData"))



RNA=t(adj.counts)
head(RNA[1:5,1:5])
rm(adj.counts)

load("validation/pureMval.RData")

load("DMethDE022724.RData")

meth=mval
Meth=t(meth)
head(Meth[1:5,1:5])
rm(meth)

#de.genes=DE
de.genes=DMGs
# de.genes=unique(de.genes[de.genes$contrast %in% c(contrs) & 
#                     abs(de.genes$logFC) >= 1 & de.genes$adj.P.Val < 0.05 ,7])

#load( "validation/FeatureNewIdea.RData")


de.genes=unique(de.genes[de.genes$contrast %in% c("KpTp - KnTn", "KpTn - KnTn") & 
        abs(de.genes$logFC) >= 0.5 & de.genes$adj.P.Val < 0.05 ,7])


#de.genes=unique(de.genes[de.genes$gene %in% genefeatures & 
 #                          abs(de.genes$logFC) >= 0.5 & de.genes$adj.P.Val < 0.05 ,7])

#DP=DVP

load("DFeatureNewIdea.RData")
DP=unique((DP[DP$contrast  %in% c("KpTp - KnTn", "KpTn - KnTn") &    abs(DP$logFC) >= 0.1 & DP$adj.P.Val < 0.05 ,7]))

#DP=unique((DP[DP$probe %in% probefeatures  &
            #    abs(DP$logFC) >= 0.1 & DP$adj.P.Val < 0.05 ,7]))

#save(promoter.probe, file="promoter.probe.RData")

load("promoter.probe.RData")


DP=intersect(DP, promoter.probe)

ids=intersect(Clinical$id, rownames(Meth))
ids=intersect(ids, rownames(RNA))
Clinical=Clinical[Clinical$id %in% ids,]
Meth=Meth[Clinical$id,DP]

RNA=RNA[Clinical$id,intersect(colnames(RNA),de.genes )]

Meth <- replace(Meth, is.infinite(Meth), 0)
Meth <- replace(Meth, is.nan(Meth), 0)
Meth <- replace(Meth, is.na(Meth), 0)

RNA[is.na(RNA)] <- 0

X=list( RNA, Meth)

Y=Clinical$class
Y=Y+1

omics= filter.unsupervised(X, method = "variance",
                           pct.keep = 60,
                           center = FALSE,
                           scale = FALSE,
                           standardize = T,
                           log2TransForm = FALSE,
                           Xtest = NULL)


#rm(Meth, RNA)

Xdata1 <- omics$X[[1]]
Xdata2 <- omics$X[[2]] 
hist(Xdata1[,"TMOD3"])
hist(Xdata2[,"cg06952671"])
dim(Xdata1)
dim(Xdata2)

snv = readRDS(file.path("/mmfs1/home/david.adeleke/Thesis","dataset/SnvPathPAAD-US.rds"))
snv[snv == 9] <- 0
snv= snv[which(rowSums(snv) >=8),]
snv=snv[,ids]
snv=t(snv)
snv=snv[,c("KRAS", "TP53", "C9orf53")]
dim(snv)

Xdata=list( Xdata1,  Xdata2)

fit.cvsida <- cvSIDA(Xdata, Y,
                     withCov = F,
                     Xtestdata = Xdata,
                     Ytest = Y, plotIt = F,
                     isParallel = TRUE,
                     ncores = 32,
                     AssignClassMethod="Joint")

rsd=data.frame()

rsd0= VarImportancePlot(fit.cvsida)
rsd1=rsd0[[1]][[1]]
rsd1=rsd1[rsd1$`Absolute Loading` !=0,]
rsd1$grp= "KpTp - KnTn"

rsd2=rsd0[[1]][[2]]
rsd2=rsd2[rsd2$`Absolute Loading` != 0,]
rsd2$grp= "KpTp - KnTn"
  
rsd=rbind(rsd,rsd1)
rsd=rbind(rsd,rsd2)
  
Xdata1b=  Xdata1[, rsd[rsd$View==1,2]]
Xdata2b= Xdata2[, rsd[rsd$View==2,2]]
Xdatab= list( Xdata1b,  Xdata2b)



fit.cvsidab <- cvSIDA(Xdatab, Y,
                     withCov = F,
                     Xtestdata = Xdatab,
                     Ytest = Y, plotIt = F,
                     isParallel = TRUE,
                     ncores = 32,
                     AssignClassMethod="Joint")

  
colpal= c( "#F8766D",  "#E76BF3")
  CorrelationPlots(Xdata, Ytest=Y, fit.cvsida$hatalpha, color.palette =  colpal)
  
 netobj=networkPlot(fit.cvsidab,cutoff=0.5) 
  
  DiscriminantPlots(Xdatab, Y, fit.cvsidab$hatalpha)
  
  rsd=data.frame()
  
  rsd0= VarImportancePlot(fit.cvsidab)
  rsd1=rsd0[[1]][[1]]
  rsd1=rsd1[rsd1$`Absolute Loading` !=0,]
  rsd1$grp= "KpTp - KnTn"
  
  rsd2=rsd0[[1]][[2]]
  rsd2=rsd2[rsd2$`Absolute Loading` != 0,]
  rsd2$grp= "KpTp - KnTn"
  
  rsd=rbind(rsd,rsd1)
  rsd=rbind(rsd,rsd2)
########################################################################################

save(rsd,fit.cvsidab,Y,Xdatab, Xdata,Clinical, file="/mmfs1/home/david.adeleke/Thesis/dataset/chapt2out.RData")
  
Clinical=readRDS("clinical_New.rds")

load( "validation/adj.counts022724.RData")
load("pureMval.RData")
RNAlogCPM=t(adj.counts)

meth= t(mval); rm(mval)

cid= intersect(intersect(Clinical$id, 
    rownames(RNAlogCPM)), rownames(meth))

rownames(Clinical)=Clinical$id
Clinical=Clinical[Clinical$id %in% cid,]

mem1=RNAlogCPM[cid,unique(rsd[rsd$View==1,2])]
mem2=meth[cid,unique(rsd[rsd$View==2,2])]

mem= merge(mem1, mem2, by=0); rm(mem1, mem2)
rownames(mem)=mem$Row.names
mem=mem[,-1]

mem.clin =  merge(mem, Clinical, by=0)

mem.clinb=mem.clin 


#save(mem.clin, file="mem.clin.RData")



# ---------------------------LASSO REGRESSION------------------------------------------#
rsd0= VarImportancePlot(fit.cvsida)
rsd1=rsd0[[1]][[1]]

rsd2=rsd0[[1]][[2]]

goi = rsd1[,2][1:20]
poi = rsd2[,2][1:20]

foi = c(goi, poi)

library(glmnet)
load("mem.clin.RData")

mem.clin5=mem.clin[mem.clin$time > 5 & mem.clin$grp %in% c("KpTp" , "KnTn"),]

endl= ncol(mem.clin5)-5

allcol= colnames(mem.clin5[, 2:endl])
pint <- colnames(mem.clin5[, grepl("^cg", names(mem.clin5))])

gint =setdiff(allcol,pint)
foi= c(pint, gint)

x <- as.matrix(mem.clin5[, pint])
y <- Surv(mem.clin5$time, mem.clin5$status)
lasso_model <- cv.glmnet(x, y, family = "cox")

# Plot the cross-validated mean squared error (optional)
plot(lasso_model)

# Identify selected genes based on the optimal lambda value
optimal_lambda <- lasso_model$lambda.min
selected_probes <- coef(lasso_model, s = optimal_lambda)
selected_probes <- selected_probes[selected_probes[, 1] != 0, ]

# Print the selected genes
print(selected_probes)

topp = names(selected_probes)

# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(topp, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = mem.clin5[mem.clin5$grp %in% c("KpTp" , "KnTn"),])

summary(modelph)
suppressWarnings(ggforest(modelph))


x <- as.matrix(mem.clin5[, gint])
y <- Surv(mem.clin5$time, mem.clin5$status)
lasso_model <- cv.glmnet(x, y, family = "cox")

# Plot the cross-validated mean squared error (optional)
plot(lasso_model)

# Identify selected genes based on the optimal lambda value
optimal_lambda <- lasso_model$lambda.min
selected_genes <- coef(lasso_model, s = optimal_lambda)
selected_genes <- selected_genes[selected_genes[, 1] != 0, ]

# Print the selected genes
print(selected_genes)
topg=  names(selected_genes)

# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(topg, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = mem.clin5[mem.clin5$grp %in% c("KpTp" , "KnTn"),])

summary(modelph)
suppressWarnings(ggforest(modelph))

# Define the list of variables
goi <- c(topp, topg) 
#,

#"cg24994127",
goi = c( "SLC9A3", "BTBD6", "SLC39A10", "ANKRD57", "PADI1" )
mem.clindf= mem.clin5[, c(goi, "id", "chemo", "status", "time", "grp")]

#mem.clindf[, 2:6] <- mem.clindf[, 2:6] * mem.clindf[, 1]


#goi = c("GABRA2", "SLC9A3", "BTBD6", "SLC39A10", "ANKRD57" )
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))

modelph <- coxph(formula_str, data = mem.clindf)


# Fit Cox proportional hazards model
#modelph <- coxph(formula_str, data = mem.clin5[mem.clin5$grp %in% c("KpTp" , "KnTn"),])

summary(modelph)
suppressWarnings(ggforest(modelph))


###---------------------####################---------------------
goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" ,"cg06952671" ,"cg10661002")



Clinical=readRDS("clinical_New.rds")

load( "validation/adj.counts022724.RData")
load("pureMval.RData")
RNAlogCPM=t(adj.counts)

#load("DNAmet2224.RData")

meth= t(mval); rm(mval)

cid= intersect(intersect(Clinical$id, 
                         rownames(RNAlogCPM)), rownames(meth))

rownames(Clinical)=Clinical$id
Clinical=Clinical[Clinical$id %in% cid,]

mem1=RNAlogCPM[cid, goi]
mem2=meth[cid,poi]

mem= merge(mem1, mem2, by=0); rm(mem1, mem2)
rownames(mem)=mem$Row.names
mem=mem[,-1]

mem.clin =  merge(mem, Clinical, by=0)


library(factoextra)
#
#goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" )
#"cg24994127",
#goi = c( "cg24994127", "SLC9A3", "BTBD6", "SLC39A10", "ANKRD57", "PADI1" )

#goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" ,"cg06952671" ,"cg10661002")

df=mem.clin[,foi]

# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)

gap_stat <- cluster::clusGap(df, FUN = hcut, K.max = 10, B = 10)

fviz_gap_stat(gap_stat)

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df, centers = 2)

# Get the cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster


mem.clin$subtype=cluster_assignments

mem.clin$subtype= ifelse(mem.clin$subtype == 1, "Group_A",
                         ifelse(mem.clin$subtype == "2", "Group_B", "Group_C"))

table(mem.clin$grp, mem.clin$subtype)


sfit <- survfit(Surv(time, status)~subtype, data=mem.clin)

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


####---------------------------------------------------------------------------###






















filepath = c("KpTp130.RData","KpTn130.RData" ,"KpTnKptp130.RData", "none")

colpal= c(  "#F8766D",  "#00BF7D" , "#00B0F6" , "#E76BF3" )

nampal = c("red", "green", "blue", "pink")

mbs = c("KnTn","KnTp","KpTn","KpTp")
contrs= c("KpTp-KnTn", "KpTn-KnTn", "KpTp-KpTn", "none")

metdf= data.frame(colpal, nampal, mbs, filepath, contrs);rm(filepath, colpal, nampal, mbs, contrs)

ii = 1

#colpal= metdf[,"colpal"][c(1,4)] 
rsd=data.frame()

for(ii in 1:3){
#cvsida = load(metdf[,"filepath"][ii])
rsd0= VarImportancePlot(fit.cvsida)
rsd1=rsd0[[1]][[1]]
rsd1=rsd1[rsd1$Loading !=0,]
rsd1$grp= metdf$contrs[ii]

rsd2=rsd0[[1]][[2]]
rsd2=rsd2[rsd2$Loading != 0,]
rsd2$grp= rsd1$grp= metdf$contrs[ii]

rsd=rbind(rsd,rsd1)

rsd=rbind(rsd,rsd2)

}



#save(rsd, file="validation/rsd.RData")


library(ggplot2)
library("ggpubr")
library(survminer)
library(survival)
library(utils)

#load("validation/rsd.RData")
Clinical=readRDS("clinical_New.rds")


load( "validation/adj.counts1221.RData")
RNAlogCPM=t(adj.counts)


load("DNAval20624.RData")

#load("DNAmet2224.RData")

meth= t(mval); rm(mval)

cid= intersect(intersect(Clinical$id, 
                rownames(RNAlogCPM)),
               rownames(meth))

rownames(Clinical)=Clinical$id
Clinical=Clinical[Clinical$id %in% cid,]

#RNAlogCPM=t(adj.counts)
mem1=RNAlogCPM[cid,unique(rsd[rsd$View==1,2])]
mem2=meth[cid,unique(rsd[rsd$View==2,2])]

mem= merge(mem1, mem2, by=0); rm(mem1, mem2)
rownames(mem)=mem$Row.names
mem=mem[,-1]

mem.clin =  merge(mem, Clinical, by=0)

save(mem.clin, file="mem.clin.RData")



# ---------------------------LASSO REGRESSION------------------------------------------#
load("mem.clin.RData")
mem.clin5=mem.clin[mem.clin$time > 5,]

endl= ncol(mem.clin5)-5
#342
allcol= colnames(mem.clin5[, 2:endl])
pint <- colnames(mem.clin5[, grepl("^cg", names(mem.clin5))])

gint =setdiff(allcol,pint)
foi= c(pint, gint)

library(glmnet)
x <- as.matrix(mem.clin5[, pint])
y <- Surv(mem.clin5$time, mem.clin5$status)
lasso_model <- cv.glmnet(x, y, family = "cox")

# Plot the cross-validated mean squared error (optional)
plot(lasso_model)

# Identify selected genes based on the optimal lambda value
optimal_lambda <- lasso_model$lambda.min
selected_genes <- coef(lasso_model, s = optimal_lambda)
selected_genes <- selected_genes[selected_genes[, 1] != 0, ]

# Print the selected genes
print(selected_genes)

topp= names(selected_genes)
topg=  names(selected_genes)
tgoi= c("FGD6", "DPP4",  "SLC4A11" )
# Define the list of variables
goi <- c(topp, topg) 

goi=rsd[rsd$View ==1,]

goi = rsd[rsd$View == 2,]
goi = head(goi[order(-goi$`Absolute Loading`),2],20)
goi= tail(goi,10)
goi= c("SLC9A3", "BTBD6" ,  "LPCAT2", "cg03758150")
#goi<- c("SLC9A3", "PCSK6", "cg11462533" , "cg21511203", "cg17701373")

# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = mem.clin5)

summary(modelph)
suppressWarnings(ggforest(modelph))

#---------------------------------------------------------------------

Clinical=readRDS("clinical_New.rds")

nomclin= Clinical[1:3,]
normid= c("DO32751N", "DO32769N", "DO50272N")

nomclin$id= normid
rownames(nomclin)= normid
nomclin$grp= "WType"
Clinical=rbind(Clinical, nomclin)


load( "validation/adj.counts1221.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)

poi= c("cg06952671" ,"cg10661002")

goi= c("TMOD3", "SLC9A3","BTBD6","CGB7" , "TLX1" ,  "MAPK15" )

#"LDB3", "HDC", "LPL", "PDE3B" ,
#, "cg12741420", "cg09453116", "cg27499361"

mem.traits.phen=RNAlogCPM[,goi]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,poi]
Dmem.traits.phen=merge(meth, Clinical, by=0)
allmem.traits.phen = merge(mem.traits.phen[,1:(ncol(mem.traits.phen)-5)], Dmem.traits.phen, by="Row.names")
allmem.traits.phen$grp= as.character(allmem.traits.phen$grp)

allmem.traits.phen$grp = ifelse(allmem.traits.phen$grp  == "KnTnHB", 
                                "KnTn", allmem.traits.phen$grp)

allmem.traits.phen = allmem.traits.phen[allmem.traits.phen$grp != "WType",]


library(factoextra)
#
#goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" )

load("mem.clin.RData")

df=mem.clin[,c(2:269)]

df=mem.clin[,goi]

# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 3, linetype = 2)


gap_stat <- cluster::clusGap(df, FUN = hcut, K.max = 10, B = 10)
fviz_gap_stat(gap_stat)

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df, centers = 3)

# Get the cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster


mem.clin$subtype=cluster_assignments

mem.clin$subtype= ifelse(mem.clin$subtype == 1, "Group_A",
                                ifelse(mem.clin$subtype == "2", "Group_B", "Group_C"))

table(mem.clin$grp, mem.clin$subtype)


sfit <- survfit(Surv(time, status)~subtype, data=mem.clin)

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


selected_genes = colnames((allmem.traits.phen))[2:7]



goi <- selected_genes 

goi = c("TMOD3", "SLC9A3", "BTBD6", "CGB7" , "MAPK15")
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = allmem.traits.phen)

summary(modelph)
suppressWarnings(ggforest(modelph))

###----------------------------------------------

#load("validation/rsd.RData")
Clinical=readRDS("clinical_New.rds")

load( "validation/adj.counts1221.RData")


load("DNAmet011124.RData")

meth= t(mval); rm(mval)

cid= intersect(intersect(Clinical$id, 
                         rownames(RNAlogCPM)),
               rownames(meth))

rownames(Clinical)=Clinical$id
Clinical=Clinical[Clinical$id %in% cid,]

RNAlogCPM=t(adj.counts)
mem1=RNAlogCPM[cid,goi[1:7]]
mem2=meth[cid,goi[8:9]]

mem= merge(mem1, mem2, by=0); rm(mem1, mem2)
rownames(mem)=mem$Row.names
mem=mem[,-1]

mem.clin =  merge(mem, Clinical, by=0)




























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


















covariates= unique(impvar$`Variable Name`)
Dcovariates= unique(rsd2$`Variable Name`)
#covariates<- paste0("ME", covariates)








#-------------------------------PLOT-----------------------------------

Clinical=readRDS("clinical_New.rds")

nomclin= Clinical[1:3,]
normid= c("DO32751N", "DO32769N", "DO50272N")

nomclin$id= normid
rownames(nomclin)= normid
nomclin$grp= "WType"
Clinical=rbind(Clinical, nomclin)


load( "validation/adj.counts1221.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=RNAlogCPM[,unique(impvar$`Variable Name`)]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,unique(rsd2$`Variable Name`)]
Dmem.traits.phen=merge(meth, Clinical, by=0)
allmem.traits.phen = merge(mem.traits.phen[,1:(ncol(mem.traits.phen)-5)], Dmem.traits.phen, by="Row.names")


library(ggplot2)
library("ggpubr")
library(survminer)
library(survival)
library(utils)

#if(!require(devtools)) install.packages("devtools")

#devtools::install_github("kassambara/ggpubr")

mem.traits.phen$grp =as.character(mem.traits.phen$grp )
Dmem.traits.phen$grp =as.character(Dmem.traits.phen$grp )

mem.traits.phen$grp=ifelse(mem.traits.phen$grp =="KnTnHB", "KnTn",  mem.traits.phen$grp)

Dmem.traits.phen$grp=ifelse(Dmem.traits.phen$grp =="KnTnHB", "KnTn",  Dmem.traits.phen$grp)


mem.traits.phen$grp <- factor(mem.traits.phen$grp, levels=c("WType","KnTn",  "KnTp",   "KpTn",   "KpTp" ))

Dmem.traits.phen$grp <- factor(Dmem.traits.phen$grp, levels=c("WType","KnTn",  "KnTp",   "KpTn",   "KpTp" ))


#compare_means(DHRS3 ~ grp, data = mem.traits.phen)

iii=1


ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Gene Expression: ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)


ggplot(Dmem.traits.phen, aes(x = grp, y =get(rsd2$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Methylation:  ",  rsd2$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 1)

#-------------------------SURVIVAR-------------------------------------------------------------------#

covariates= unique(impvar$`Variable Name`)
Dcovariates= unique(rsd2$`Variable Name`)
#covariates<- paste0("ME", covariates)

mem.traits.phenbk=mem.traits.phen

mem.traits.phen=mem.traits.phen[mem.traits.phen$grp %in% c(dustin, jwd),]

Dmem.traits.phen=Dmem.traits.phen[Dmem.traits.phen$grp %in% c(dustin, jwd),]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mem.traits.phen)})

{
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  
  p.df <- t(as.data.frame(univ_results, check.names = FALSE))
  p.df=as.data.frame(p.df)
  
  p.df= p.df[p.df$p.value <=0.05,]
  
  
  library(kableExtra)
  p.df %>%
    kbl() %>%
    kable_styling()
  
  sig.modu.univ= rownames(p.df[p.df$p.value <=0.05,])
  
  
  ## generate all comginations of module
  sig.modu.univ.combo<- combn(sig.modu.univ, 2, FUN = NULL, simplify = FALSE)
  
  
  c.df <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(c.df) <- c( "Mod.A","Mod.B", "pvalue.A","pvalue.B",  "HR.A","HR.B", "CI.A", "CI.B")
  
  for(gen in sig.modu.univ.combo){
    d=unlist((gen))
    g1=(d[1:1])
    g2=(d[2:2])
    
    ddh=coxph(Surv(time, status) ~ get(g1) + get(g2), data = mem.traits.phen)
    
    Mod.A=g1
    Mod.B=g2
    pvalue.A= summary(ddh)$coefficients[1, 5]
    pvalue.B=summary(ddh)$coefficients[2, 5]
    HR.A=round(exp(coef(ddh))[1],2)
    HR.B=round(exp(coef(ddh))[2],2)
    CI.A=paste0(round(exp(confint(ddh))[1,1],2), " , " ,round(exp(confint(ddh))[1,2],2))
    CI.B=paste0(round(exp(confint(ddh))[2,1],2), " , " ,round(exp(confint(ddh))[2,2],2))
    
    df<- data.frame(matrix(ncol = 8, nrow = 0))
    
    colnames(df) <- c( "Mod.A","Mod.B", "pvalue.A","pvalue.B",  "HR.A","HR.B", "CI.A", "CI.B" )
    
    df[1,] <- c(Mod.A, Mod.B, round(pvalue.A,2), round(pvalue.B,2), HR.A,  HR.B, CI.A, CI.B)
    
    assign("c.df",rbind(c.df,df), envir = .GlobalEnv)
    
  }
  
  c.df=c.df[c.df$pvalue.A <0.05 & c.df$pvalue.B <0.05, ]
  
  c.df %>%
    kbl() %>%
    kable_styling()
  
}
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


# Define the list of variables
goi <- names(selected_genes)  #"cg11462533" "cg23560159" "cg21200382"
goi= goi[c(1,4,5)]
# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = Dmem.traits.phen)

summary(modelph)
suppressWarnings(ggforest(modelph))

#-------------------

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
goi <- names(selected_genes)[3:length(selected_genes)]
goi=goi[-c(2,6,8,9)]  #"SLC9A3" "BTBD6"  "FASN" , "DHRS3",  "CDK6"   "PCSK6"  "TLX1"  Top six 
# Create a formula dynamically

formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model


mem.traits.phendk=  mem.traits.phenbk[ mem.traits.phenbk$grp %in% c("KpTp", "KnTn"),]
mem.traits.phendk$grp = ifelse(mem.traits.phendk$grp == "KnTnHB", 
                               "KnTn",  mem.traits.phendk$grp )
modelph <- coxph(formula_str, data = mem.traits.phendk)

summary(modelph)


suppressWarnings(ggforest(modelph))





DRS= -0.41545*SLC9A3  -0.79314*BTBD6 -0.33056*FASN +0.62441*CDK6 + 0.44733*PCSK6 + 0.75832*TLX1 


#-------------------------------------------------Risk Score______________________________________
goi= c("SLC9A3" ,"BTBD6" , "FASN" , "DHRS3",  "CDK6" ,  "PCSK6" , "TLX1")
poi= c( "cg11462533", "cg23560159", "cg21200382")
features =  c(goi, poi)




Clinical=readRDS("clinical_New.rds")

nomclin= Clinical[1:3,]
normid= c("DO32751N", "DO32769N", "DO50272N")

nomclin$id= normid
rownames(nomclin)= normid
nomclin$grp= "WType"
Clinical=rbind(Clinical, nomclin)


load( "validation/adj.counts1221.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=RNAlogCPM[,unique(impvar$`Variable Name`)]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,unique(rsd2$`Variable Name`)]
Dmem.traits.phen=merge(meth, Clinical, by=0)
allmem.traits.phen = merge(mem.traits.phen[,1:(ncol(mem.traits.phen)-5)], Dmem.traits.phen, by="Row.names")




mem.traits.phen=allmem.traits.phen[,c(colnames(allmem.traits.phen)[137:141],features)]

mem.traits.phen$DRS= with(mem.traits.phen, -0.45322*(SLC9A3) -0.74002*(BTBD6) -0.54414*(FASN)  
                          + 0.89434*(TLX1) + 0.79992*(CDK6) + 0.63090*(PCSK6) - 0.22428*(cg11462533) 
                          -0.38313*(cg23560159)+ 0.30200*(cg21200382))


library(fabricatr)
qsubtype<-split_quantile(mem.traits.phen$DRS, type = 3)
mem.traits.phen<- cbind(mem.traits.phen,qsubtype)


library(ggplot2)
library("ggpubr")
library(survminer)
library(survival)
library(utils)

ggplot(mem.traits.phen, aes(x = grp, y =DRS, fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Risk Score: "))+
  stat_compare_means(method = "anova", label.y = 20)

mem.traits.phen$DRSs=mem.traits.phen$DRS^2
mem.traits.phen=mem.traits.phen[mem.traits.phen$grp != "WType",]
res.cut <- surv_cutpoint(mem.traits.phen, time = "time", event = "status",  variables = c("DRS"))
plot(res.cut, "DRS", palette = "npg")

mem.traits.phen$DRScat = ifelse(mem.traits.phen$DRS<= 2.9, "High", "Low")

mem.traits.phen$grp= as.character(mem.traits.phen$grp)

mem.traits.phen$grp= ifelse(mem.traits.phen$grp=="KnTnHB", "KnTn", mem.traits.phen$grp)



res.cat <- surv_categorize(res.cut)
sfit <- survfit(Surv(time, status)~DRS, data=res.cat[res.cat$time >5,])

#sfit <- survfit(Surv(time, status)~DRScat, data=mem.traits.phen)


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
  surv.median.line = "hv"  # add the median survival pointer.
  #legend.labs = 
  #c("Low Risk", "High Risk"),    # change legend labels.
  #palette =  c(colpal) # custom color palettes.
)














Clinical=readRDS("clinical_New.rds")

nomclin= Clinical[1:3,]
normid= c("DO32751N", "DO32769N", "DO50272N")

nomclin$id= normid
rownames(nomclin)= normid
nomclin$grp= "WType"
Clinical=rbind(Clinical, nomclin)


load( "validation/adj.counts1221.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=[,features]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)
rownames(mem.traits.phen)=mem.traits.phen$Row.names
mem.traits.phen=mem.traits.phen[,-1]
meth=meth[, features[8:9]]
mem.traits.phen=merge(meth,mem.traits.phen, by=0)

mem.traits.phen$DRS= with(mem.traits.phen, 0.45222*(TMOD3)  -0.32910*(SLC9A3) -0.91721*(BTBD6) -0.83312*(FASN)  
                          -1.01812*LOC647121 + 0.37013*CGB7 + 0.88256*TLX1 + 0.3902*(cg06952671) -0.2715*(cg10661002))


library(fabricatr)
qsubtype<-split_quantile(mem.traits.phen$DRS, type = 3)
mem.traits.phen<- cbind(mem.traits.phen,qsubtype)




library(ggplot2)
library("ggpubr")
library(survminer)
library(survival)
library(utils)

ggplot(mem.traits.phen, aes(x = grp, y =DRS, fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Risk Score: "))+
  stat_compare_means(method = "anova", label.y = 20)

mem.traits.phen$DRSs=mem.traits.phen$DRS^2
mem.traits.phen=mem.traits.phen[mem.traits.phen$grp != "WType",]
res.cut <- surv_cutpoint(mem.traits.phen, time = "time", event = "status",  variables = c("DRSs"))
plot(res.cut, "DRSs", palette = "npg")

mem.traits.phen$DRScat = ifelse(mem.traits.phen$DRSs<= 394.42, "High", "Low")

mem.traits.phen$grp= as.character(mem.traits.phen$grp)

mem.traits.phen$grp= ifelse(mem.traits.phen$grp=="KnTnHB", "KnTn", mem.traits.phen$grp)


ggplot(mem.traits.phen, aes(x = grp, y =1/sqrt(DRSs), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#f4f4f4"))+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Risk Score: "))+
  stat_compare_means(method = "anova", label.y = 0.1, label.x = 2)


res.cat <- surv_categorize(res.cut)
sfit <- survfit(Surv(time, status)~DRSs, data=res.cat)

#sfit <- survfit(Surv(time, status)~qsubtype, data=mem.traits.phen)


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
  surv.median.line = "hv"  # add the median survival pointer.
  #legend.labs = 
  #c("Low Risk", "High Risk"),    # change legend labels.
  #palette =  c(colpal) # custom color palettes.
)


#-------------------------------------

# Extracting the features for clustering


goi1= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1")
goi2= c("cg06952671" ,"cg10661002")


load( "validation/adj.counts1221.RData")

#load("DNAmet122323.RData")

load("DNAmet011124.RData")

meth= t(mval)


rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=RNAlogCPM[,unique(goi1)]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,unique(goi2)]
Dmem.traits.phen=merge(meth, Clinical, by=0)
allmem.traits.phen = merge(mem.traits.phen[,1:(ncol(mem.traits.phen)-5)], Dmem.traits.phen, by="Row.names")




df <- allmem.traits.phen[, c("id",goi)]
rownames(df)=df$id

df= df[, -1]


library(factoextra)

# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df, centers = 3)

# Get the cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster


allmem.traits.phen$subtype=cluster_assignments

allmem.traits.phen$subtype= ifelse(allmem.traits.phen$subtype == 1, "Group_A",
                                ifelse(allmem.traits.phen$subtype == "2", "Group_B", "Group_C"))

library(survminer)
library(survival)
sfit <- survfit(Surv(time, status)~subtype, data=allmem.traits.phen)

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

allmem.traits.phen$grp=as.character(allmem.traits.phen$grp)
allmem.traits.phen$grp= ifelse(allmem.traits.phen$grp=="KnTnHB", "KnTn", allmem.traits.phen$grp)

tb.df=allmem.traits.phen[allmem.traits.phen$grp != "WType",c(15,16)]


tb.df = tb.df %>%
  group_by(grp, subtype) %>%
  summarize(frequency = n())
tb.df$percentage = (tb.df$frequency/110)*100

ggplot(tb.df, aes(x=subtype, y=frequency, fill=grp)) +
  geom_col(position="dodge")


ggplot(tb.df, aes(x=subtype, y=percentage, fill=grp)) + geom_col()


ggplot(tb.df, aes(x = subtype, y = percentage, fill = grp)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(title = "Histogram of multi-omics based subtypes Frequencies by SNV types",
       x = "Cluster group",
       y = "Percentage",
       fill = " SNV subtype") +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  theme_minimal()










require("survival")
model <- coxph( Surv(time, status) ~ TMOD3 + SLC9A3 + BTBD6 + FASN + 
                  LOC647121 + CGB7 + TLX1, data = mem.traits.phen[mem.traits.phen$grp %in% c("KnTn", "KpTp"),] )

suppressWarnings(ggforest(model))



#if(!require(devtools)) install.packages("devtools")

#devtools::install_github("kassambara/ggpubr")

#--------------------Kmeans Clusttering -------------------------

# Extracting the features for clustering
features =  c("TMOD3","SLC9A3", "BTBD6", "FASN", "LOC647121", "CGB7", "TLX1", "cg06952671" ,"cg10661002")

df <- mem.traits.phen[, c("Row.names",features)]
rownames(df)=df$Row.names

df= df[, -1]


library(factoextra)

# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df[,1:9], centers = 3)

# Get the cluster assignments for each sample
cluster_assignments <- kmeans_result$cluster


mem.traits.phen$subtype=cluster_assignments

mem.traits.phen$subtype= ifelse(mem.traits.phen$subtype == 1, "Group_A",
                                ifelse(mem.traits.phen$subtype == "2", "Group_B", "Group_C"))


ggplot(mem.traits.phen, aes(x = subtype, y =1/sqrt(DRS^2), fill = subtype)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#f4f4f4"))+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Risk Score: "))+
  stat_compare_means(method = "anova", label.y = 0.1, label.x = 2)




sfit <- survfit(Surv(time, status)~subtype, data=mem.traits.phen)

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

model <- coxph( Surv(time, status) ~ qsubtype, data = mem.traits.phen[mem.traits.phen$time >10,] )

suppressWarnings(ggforest(model))

mem.traits.phen$grp= as.character(mem.traits.phen$grp)
mem.traits.phen$grp= ifelse(mem.traits.phen$grp=="KnT", "KnTn", mem.traits.phen$grp )

mem.traits.phen=mem.traits.phen[mem.traits.phen$grp != "WType",]
table(mem.traits.phen$grp, mem.traits.phen$subtype)


ref="Group_A"
mem.traits.phen$subtype= as.factor(mem.traits.phen$subtype)
mem.traits.phen <- within(mem.traits.phen, subtype <- relevel(subtype, ref = ref))
ttt=coxph(Surv(time, status) ~ subtype, data = mem.traits.phen) %>% 
  gtsummary::tbl_regression(exp = TRUE) 


#------------------------------------------------------------------------------------------
#Try weighing using features using LASSO coefficient



#Try weighing using Features using Sandra coefficients

######################################################################################


mem.traits.phen$grp =as.character(mem.traits.phen$grp )

Dmem.traits.phen$grp =as.character(Dmem.traits.phen$grp)

mem.traits.phen$grp=ifelse(mem.traits.phen$grp =="KnTnHB", "KnTn",  mem.traits.phen$grp)

Dmem.traits.phen$grp=ifelse(Dmem.traits.phen$grp =="KnTnHB", "KnTn",  Dmem.traits.phen$grp)


mem.traits.phen$grp <- factor(mem.traits.phen$grp, levels=c("WType","KnTn",  "KnTp",   "KpTn",   "KpTp" ))

Dmem.traits.phen$grp <- factor(Dmem.traits.phen$grp, levels=c("WType","KnTn",  "KnTp",   "KpTn",   "KpTp" ))


#compare_means(DHRS3 ~ grp, data = mem.traits.phen)

iii=1


ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Gene Expression: ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)



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
