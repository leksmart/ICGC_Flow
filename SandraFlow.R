library(mvlearnR)

library(survminer)
library(survival)
library(utils)

library(ggplot2)
library("ggpubr")

setwd("/mmfs1/home/david.adeleke/WorkFlow/Chapter_Two/")

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

contrs= paste0(dustin, " - ", jwd)
#contrs= paste0(jwd, " - ", dustin)
contrx= paste0(nat, " - ", dustin)

contrx=""
nat=""
#"KnTp - KnTn"
#k cc dustin,
Clinical$grp=as.character(Clinical$grp)
Clinical$grp = ifelse(Clinical$grp == "KnTnHB", "KnTn",  Clinical$grp )

Clinical$class= ifelse(Clinical$grp %in% c(dustin),0,
                       ifelse(Clinical$grp %in% c(jwd),1,
                              ifelse(Clinical$grp %in% c(nat),3,3)))


Clinical= Clinical[ Clinical$class < 3,]


load( "validation/adj.counts1221.RData")

load( "validation/DMGsraw1221.RData")

adj.counts<- edgeR::cpm(adj.counts, log=TRUE)

RNA=t(adj.counts)
rm(adj.counts)


load("DNAmet122323.RData")

load("DMethDE122323.RData")



Meth=t(meth)
rm(meth)



de.genes=DE

de.genes=de.genes[de.genes$contrast %in% c(contrs,contrx) & 
                    abs(de.genes$logFC) >= 0.5 & de.genes$adj.P.Val < 0.05 ,7]

DP=DP[DP$contrast %in% c(contrs,contrx)  & 
        abs(DP$logFC) >= 0.2 & DP$adj.P.Val < 0.05 ,7]


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
                           pct.keep = 75,
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
# myresult=selpscca.pred(Xdata1, Xdata2, Y=Clinical$time, fitselpCCA=NULL, family="survival",
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


CorrelationPlots(Xdata, Ytest=Y, fit.cvsida$hatalpha)

needt=networkPlot(fit.cvsida,cutoff=0.67)

mymat = networkplotinner(fit.cvsida)

myComb = mymat$ViewCombinations
nComb = dim(mymat$ViewCombinations)[2]


networkPlot(fit.cvsida,cutoff=0.55)

DiscriminantPlots(Xdata, Y, fit.cvsida$hatalpha)

rsd=data.frame()
rsd1= VarImportancePlot(fit.cvsida)
rsd1=rsd1[[1]][[1]]

rsd1=rsd1[rsd1$Loading !=0,]

#rsd1=rsd1[rsd1$Loading > 0.01,]

rsd1$grp= paste0(w," - ", j)

rsd2= VarImportancePlot(fit.cvsida)
rsd2=rsd2[[1]][[2]]

rsd2=rsd2[rsd2$Loading > 0,]

rsd2$grp= paste0(w," - ", j)

rsd=rbind(rsd,rsd1)


save(rsd, file="validation/rsd.RData")


impvar= data.frame()
impvar=rbind(impvar,rsd1)


"Separate"
tdf=data.frame(fit.cvsida$PredictedClass)
tdf=cbind(tdf,Y)

table(tdf$X1, tdf$Y)
table(tdf$X2, tdf$Y)


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

rownames(Clinical)=Clinical$id
RNAlogCPM=t(adj.counts)
mem.traits.phen=RNAlogCPM[,unique(impvar$`Variable Name`)]
mem.traits.phen=merge(mem.traits.phen, Clinical, by=0)




#if(!require(devtools)) install.packages("devtools")

#devtools::install_github("kassambara/ggpubr")


mem.traits.phen$grp <- factor(mem.traits.phen$grp, levels=c("WType","KnTn", "KnTnHB",   "KnTp",   "KpTn",   "KpTp" ))


#compare_means(DHRS3 ~ grp, data = mem.traits.phen)

iii=1


ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Expression of  ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)


ggplot(mem.traits.phen, aes(x = grp, y =get(impvar$`Variable Name`[iii]), fill = grp)) +
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.07) + 
  theme(legend.position = "none")+
  xlab("Mutation-based Subtypes") + 
  ylab(paste0("Mean Expression of  ",  impvar$`Variable Name`[iii], " gene"))+
  stat_compare_means(method = "anova", label.y = 20)

#-------------------------SURVIVAR-------------------------------------------------------------------#


covariates= unique(impvar$`Variable Name`)
#covariates<- paste0("ME", covariates)
mem.traits.phen=mem.traits.phen[mem.traits.phen$grp %in% c(dustin, jwd),]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = mem.traits.phen)})

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



ddh=coxph(Surv(time, status) ~ SLC9A3 + NFE2L3  + DHRS3 +  TLX1 +  LPCAT2 + CDKN2B , data = mem.traits.phen)

ddh=coxph(Surv(time, status) ~ JSRP1 + REM1  + KCNN3 +  CSTA +  SP6 + LINGO4 + RNF157 + PCSK6 + CGB7 + NFE2L3 + IRX2, TACR1 + TDRD9 + FERMT1 + FLJ40330 , data = mem.traits.phen)

summary(ddh)

ddh=coxph(Surv(time, status) ~ SLC9A3 + NFE2L3  + TLX1 + CDKN2B , data = mem.traits.phen)


ddh=coxph(Surv(time, status) ~ SLC9A3   + SCARNA6 +GL   , data = mem.traits.phen)

# ---------------------------LASSO REGRESSION------------------------------------------#
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
#goi <-c(goi ,"grp")

# Create a formula dynamically
formula_str <- as.formula(paste("Surv(time, status) ~", paste(goi, collapse = " + ")))
# Fit Cox proportional hazards model
modelph <- coxph(formula_str, data = mem.traits.phen)

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
