library(mvlearnR)
excontrs= c("KpTn - WT", "KpTp - KnTn" ,  "KpTp - KnTp" ,  "KpTp - KpTn", "KpTp - WT")

setwd("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/")
goi= c("TMOD3", "SLC9A3","BTBD6","FASN" ,"LOC647121" ,"CGB7" , "TLX1" ,"cg06952671" ,"cg10661002")

#load("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/validation/auvalidationset.RData")

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
nat= cc
nat="KnTp"
contrs= paste0(dustin, " - ", jwd)
#contrs= paste0(jwd, " - ", dustin)
contrx= paste0(nat, " - ", dustin)


#contrx=""

#"KnTp - KnTn"
#k cc dustin,
Clinical$grp=as.character(Clinical$grp)
Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn",  Clinical$grp )

Clinical$class= ifelse(Clinical$grp %in% c(dustin),0,
                       ifelse(Clinical$grp %in% c(jwd),1,
                              ifelse(Clinical$grp %in% c(nat),2,3)))


Clin= Clinical[ Clinical$class < 3,]

#KnTn KnTp KpTp 

load( "validation/adj.counts1221.RData")

load( "validation/DMGsraw1221.RData")

#arcsin value from bvalue
#load("DNAmet122323.RData")
load("DNAmet011124.RData")

#DE from mval from arcsin of bval
load("DMethDE011124.RData")


#adj.counts<- edgeR::cpm(adj.counts, log=TRUE)

RNA=t(adj.counts)


meth=mval

Meth=t(meth)
rm(meth)



de.genes=DE[DE$contrast %in% c(contrs, contrx, "KpTn - KnTn") & 
                    abs(DE$logFC) >= 1 & DE$adj.P.Val < 0.05 ,7]

de.genes=unique(de.genes)
#DP=DVP
de.probes=DP[DP$contrast %in% c(contrs, contrx, "KpTn - KnTn")  & 
        abs(DP$logFC) >= 1.5 & DP$adj.P.Val < 0.05 ,7]

de.probes=unique(de.probes)
#save(promoter.probe, file="promoter.probe.RData")

load("promoter.probe.RData")


de.probes=intersect(de.probes, promoter.probe)

#clin.id = Clinical[Clinical$grp %in% c(dustin, jwd),"id"]
clin.id =Clin$id
ids=intersect(clin.id, rownames(Meth))
ids=intersect(ids, rownames(RNA))
Clin=Clinical[Clinical$id %in% ids,]
Meth=Meth[Clin$id,de.probes]

RNA=RNA[Clin$id,de.genes]

X=list( RNA, Meth)


dim(RNA); dim(Meth)

Y=(Clin$class)+1

omics= filter.unsupervised(X, method = "variance",
                           pct.keep = 80,
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


# 
# myresult=selpscca.pred( Xdata1, Xdata2, Y=Clinical$time, fitselpCCA=NULL, family="survival",
#                        event= Clinical$status,model.separately=FALSE, ncancorr=1,
#                        CovStructure="Iden", isParallel=TRUE, ncores=30,
#                        nfolds=5, ngrid=10, standardize=TRUE,thresh=0.0001,
#                        maxiteration=20, showProgress=T)
# 
# 
# #check output
# train.correlation=myresult$selp.fit$maxcorr
# optTau=myresult$selp.fit$optTau
# hatalpha=myresult$selp.fit$hatalpha
# hatbeta=myresult$selp.fit$hatbeta
# predictionModel=summary(myresult$mod.fit)

#, "#FF0000"
colpal= c( "#F8766D",  "#E76BF3" , "#FF0000")

CorrelationPlots(Xdata, Ytest=Y, fit.cvsida$hatalpha)

netobj=networkPlot(fit.cvsida,cutoff=0.5) 

DiscriminantPlots(Xdata, Y, fit.cvsida$hatalpha)

#"edges", "vertices",
KnTnKpTpnet=netobj[[1]]$NetworkGraph
netdf=igraph::as_data_frame(KnTnKpTpnet, what = c("both"))

netdfvt=as.data.frame(netdf$vertices)

netdfed=as.data.frame(netdf$edges)


load("DMethDE011124.RData")

load( "validation/DMGsraw1221.RData")


write.table(netdfvt, "netvtKnTnKpTp.txt",   col.names=TRUE,
            sep="\t", row.names=FALSE, quote=FALSE)


write.table(netdfed, "netedKnTnKpTp.txt",   col.names=TRUE,
            sep="\t", row.names=FALSE, quote=FALSE)


DP= DP[DP$probe %in% unique(c(netdfed$from, netdfed$to)) & DP$contrast ==contrs ,]
DE= DE[DE$gene %in% unique(c(netdfed$from, netdfed$to)) & DE$contrast ==contrs ,]

DP=DP[,c(7,1)]
DE=DE[,c(7,1)]

colnames(DE)= c("name", "FoldChange")
colnames(DP)= c("name", "FoldChange")

DD=rbind(DP, DE)

write.table(DD, "DEGPKnTnKpTp.txt",   col.names=TRUE,
            sep="\t", row.names=FALSE, quote=FALSE)



rsd=data.frame()
rsd1= VarImportancePlot(fit.cvsida)

rsd1=rsd1[[1]][[1]]

rsd1=rsd1[rsd1$`Mean Loading` !=0,]

rsd1$grp= ""

rsd2= VarImportancePlot(fit.cvsida)
rsd2=rsd2[[1]][[2]]

rsd2=rsd2[rsd2$`Mean Loading` >=0.1,]

rsd2$grp=""

rsd=rbind(rsd,rsd1)
rsd=rbind(rsd,rsd2)
save(rsd, file="validation/rsdn.RData")


impvar= rsd1

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


load( "validation/adj.counts1221.RData")

load("DNAmet011124.RData")

meth= t(mval)

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)

mem.traits.phen=RNAlogCPM[,rsd[rsd$View ==1,2]]

mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)

meth=meth[,rsd[rsd$View ==2 ,2]]
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
  stat_compare_means(method = "anova", label.y = 10)


ggplot(Dmem.traits.phen, aes(x = grp, y =get(rsd2$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Methylation:  ",  rsd2$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 1)

#-------------------------SURVIVAR-------------------------------------------------------------------#

covariates= rsd[rsd$View ==1, 2]
Dcovariates=  rsd[rsd$View ==2, 2]
#covariates<- paste0("ME", covariates)

mem.traits.phenbk=mem.traits.phen

#mem.traits.phen=mem.traits.phen[mem.traits.phen$grp %in% c(dustin, jwd),]

#Dmem.traits.phen=Dmem.traits.phen[Dmem.traits.phen$grp %in% c(dustin, jwd),]

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
Dmem.traits.phenb=Dmem.traits.phen[Dmem.traits.phen$time > 5,]

endl= ncol(Dmem.traits.phenb)-5
library(glmnet)
x <- as.matrix(Dmem.traits.phenb[, 2:endl])
y <- Surv(Dmem.traits.phenb$time, Dmem.traits.phenb$status)
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
#goi= goi[c(1,4,5)]
# Create a formula dynamically
goi=c("cg27467383", "cg11462533", "cg01891818", "cg18504358", "cg17202781" )
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = Dmem.traits.phen)

summary(modelph)
suppressWarnings(ggforest(modelph))

#---------------------------------------------------------------------------------

mem.traits.phenb=mem.traits.phen[mem.traits.phen$time > 5,]

endl= ncol(mem.traits.phenb)-5

x <- as.matrix(mem.traits.phenb[, 2:endl])
y <- Surv(mem.traits.phenb$time, mem.traits.phenb$status)
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
goi <- names(selected_genes)[1:length(selected_genes)]
#goi=goi[-c(2,6,8,9)]  #"SLC9A3" "BTBD6"  "FASN" , "DHRS3",  "CDK6"   "PCSK6"  "TLX1"  Top six 
# Create a formula dynamically
goi = c("SLC26A6", "AMIGO2", "PSMB9", "BTBD6", "SLC9A3", "RARRES3")
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model


#mem.traits.phendk=  mem.traits.phenbk[ mem.traits.phenbk$grp %in% c("KpTp", "KnTn"),]
mem.traits.phendk$grp = ifelse(mem.traits.phendk$grp == "KnTnHB", 
                               "KnTn",  mem.traits.phendk$grp )
modelph <- coxph(formula_str, data = mem.traits.phenb)

summary(modelph)


suppressWarnings(ggforest(modelph))



DRS= -0.41545*SLC9A3  -0.79314*BTBD6 -0.33056*FASN +0.62441*CDK6 + 0.44733*PCSK6 + 0.75832*TLX1 


#-------------------------------------------------Risk Score______________________________________
#goi= c("SLC9A3" ,"BTBD6" , "FASN" , "DHRS3",  "CDK6" ,  "PCSK6" , "TLX1")

goi= c("SLC9A3","SP6","BTBD6","CGB7" , "BAI1" , "TMOD3" ,"TLX1","PCSK6", "TACR1", "NFE2L3")

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




mem.traits.phen=allmem.traits.phen[,c("chemo","status", "time", "grp",features)]

mem.traits.phen$DRS= with(mem.traits.phen, -0.45322*(SLC9A3) -0.74002*(BTBD6) -0.54414*(FASN)  
                 + 0.3238*(TLX1) + 0.79992*(CDK6) + 0.63090*(PCSK6) -0.22046*(cg11462533) 
                 -0.31708*(cg23560159)+ 0.19876*(cg21200382) + 0.2682*(CGB7) -0.5535*(BAI1))


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

df <- mem.traits.phenbk[, c(features)]
rownames(df)=df$id

df= df[, -1]


library(factoextra)

# Extract the optimal number of clusters
fviz_nbclust(df, kmeans, method = "silhouette")

# Run K-means clustering with the optimal number of clusters
kmeans_result <- kmeans(df, centers = 3)

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
