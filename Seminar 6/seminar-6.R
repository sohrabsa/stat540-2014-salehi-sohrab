## Seminar 6
# ref: http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar06_highVolumeLinearModelling.html

source("http://www.bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("statmod")

library(limma)
library(lattice)

dataDir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Second Semester/Stat 540/Lab/stat540_2014/examples/photoRec/data/"
prDat <- read.table(paste(dataDir, "GSE4051_data.tsv", sep = ""))
str(prDat, max.level = 0)
prDesc <- readRDS(paste(dataDir, "GSE4051_design.rds", sep=""))
str(prDesc)

prepareData <- function(gene.list = NULL) {
  if (is.null(gene.list))
    return(NULL);
  
  filteredPrDat <- prDat[rownames(prDat) %in% gene.list, ]
  (newDat <- with(filteredPrDat, 
                  data.frame(prDesc,
                             gExp = as.vector(t(filteredPrDat)),
                             gene = factor(rep(gene.list, each = ncol(filteredPrDat)), 
                                           levels = gene.list))))
  
  
}

makeStripplot <- function(a.data.frame, pch = 19, cex = 1) {
  # check that the input has the minimum requirements
  reqVariables <- c("gExp", "devStage", "gType", "gene")
  if (!all(sapply(reqVariables, function(x) !is.null(a.data.frame[[x]])))) {
    print("Error! The input data.frame doesn't contain required fields!")
    return(NULL)
  }
  
  (f <- stripplot(gExp ~ devStage | gene, a.data.frame,
                  group = gType, jitter.data = TRUE,
#                   auto.key = TRUE, type = c('p', 'a'), grid = TRUE, pch = pch, cex = cex))
                  auto.key = TRUE, type = c('p', 'a'), grid = TRUE))
}

m <- 1000
n <- 3
x <- matrix(rnorm(m * n), nrow = m)


obsVars <- apply(x, 1, var)
summary(obsVars)
mean(obsVars < 1/3)
densityplot(~obsVars, n = 200)

# check the number of observations that have a variance less than 1/3 of the true variance
length(which(obsVars < 1/3))

# optional take home exercise 
# add two more groups, cond1, cond2, and cond3
numOfSamples <- 3
numOfGenes <- 1000
cond1 <- matrix(rnorm(numOfSamples * numOfGenes), byrow=T, ncol=3)
cond2 <- matrix(rnorm(numOfSamples * numOfGenes, mean=10, sd=3), byrow=T, ncol=3)
cond3 <- matrix(rgamma(numOfSamples * numOfGenes, shape=1, scale=3), byrow=T, ncol=3)

# check out the different models
densityplot(as.vector(cond1))
densityplot(as.vector(cond2))
densityplot(as.vector(cond3))

# stack the models together
conds <- cbind(cond1, cond2, cond3)
head(conds)

# look at the variance across all conditions
condVar <- apply(conds, 1, var)
summary(condVar)

# look at variance within each group
vars <- lapply(seq(3), function(x) { apply(conds[, (1+(x-1)*3):(x*3)], 1, var) })
lapply(vars, summary)


# Fit a linear model: explain gene expression in 
# the wild type mice as a function of developmental stage (one-way ANOVA)
wtDes <- subset(prDesc, gType == "wt")
str(wtDes)

# this is something new! Its elegant! You just specify from prDesc and the subset works like a join in associative DBs!
wtDat <- subset(prDat, select = prDesc$gType == "wt")
str(wtDat, max.level = 0)

wtDesMat <- model.matrix(~devStage, wtDes)
str(wtDesMat)

wtFit <- lmFit(wtDat, wtDesMat)
wtEbFit <- eBayes(wtFit)
topTable(wtEbFit)


# this call will work but I don't like it
topTable(wtEbFit, coef = 2:5)  # cryptic! error-prone!

colnames(coef(wtEbFit))
# extract the 
(dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit)))))

# plot hits 3, 6, and 9
gene.list <- dsHits$ID[c(3, 6, 9)]
gene.list

## optional exercise
# extract the data.frame, filter by wildtype, and then plot
t3Dat <- prepareData(gene.list)
t3Dat <- t3Dat[t3Dat$gType == "wt", ]
makeStripplot(t3Dat)

# use lm and compare F-statistic with p-values
lm.models <- lapply(gene.list, function(x) lm(formula=gExp~devStage, t3Dat, subset = gene == x))

calcPvalueForLM <- function(some.lm.fit) {
  pf(summary(some.lm.fit)$fstatistic[1], summary(some.lm.fit)$fstatistic[2], summary(some.lm.fit)$fstatistic[3], lower=FALSE)
}


(p.value <- dsHits[dsHits$ID %in% gene.list, c('F', 'P.Value')])

lmSummary <- sapply( seq(length(lm.models)), function(x) { c(summary(lm.models[[x]])$fstatistic['value'], calcPvalueForLM(lm.models[[x]])) } )
rownames(lmSummary) <- c("F", "P.Value")
colnames(lmSummary) <- gene.list
t(lmSummary)

# TODO: ask this question.
print("They don't seem similar, not even sorted. Is it expected?")

## Be the boss of topTable()
# how many q-values(?) less than 1e-5?
colnames(dsHits)
sum(dsHits$adj.P.Val < 1e-5)

# dsHits is too small, only has 10 rows
nrow(dsHits)

# tell topTable to return everyone
dsHits <- topTable(wtEbFit, coef = grep("devStage", colnames(coef(wtEbFit))), number=Inf)
nrow(dsHits)

# perform the test again, which yields 350 as requested
sum(dsHits$adj.P.Val < 1e-5)

# find out 63rd hit (assuming that it means sorted by adjusted.p.value)
dsHits[63, c("F", "adj.P.Val", "devStageP6")]

## Let's consider the effects alone (the first case-scenario according to the lectures)
niceScatter <- function(x.var, y.var, x.label, y.label) {
  lims <- c(-20, 15)
  xyplot(y.var~x.var, asp = 1, 
         panel = function(x,y, ...){
           panel.smoothScatter(x, y, ...)
           panel.abline(0, 1, col="orange")
         } , nbin = 150, ylim=lims, xlim=lims, xlab=x.label, ylab=y.label)
}

# make sure they're not sorted in a wierd way, i.e., corresponding rows belong to same IDs
p10Hits <- topTable(wtEbFit, coef = "devStageP10", number=Inf, sort.by="none")
p02Hits <- topTable(wtEbFit, coef = "devStageP2", number=Inf, sort.by="none")

niceScatter(p02Hits$t, p10Hits$t, "t-statistics for P2 effect", "t-statistic for p10 effect")

# now look at them together, its a very similar plot, but looks scaled
niceScatter(dsHits$devStageP2, dsHits$devStageP10, "t-statistics for P2 effect", "t-statistic for p10 effect")

# looks like a scaled version, eh? 
summary(p10Hits$t)
summary(dsHits$devStageP10)d

# the difference in their max.value
#16.75000/5.2790 = 3.172949

## continue with the associated p-values
tail(pss)
pss <- data.frame(vals = c(p10Hits$adj.P.Val, p02Hits$adj.P.Val), source=rep(c("p10", "p02"), each=nrow(p02Hits)) )
densityplot(~vals, pss, groups=source,  auto.key=T)

# p10 seems to have more mass under .2 thus more distant from the ref
temp.table <- data.frame("P02"= p02Hits$adj.P.Val < 1e-03, "P10"=p10Hits$adj.P.Val < 1e-03)
addmargins(table((temp.table)))

head(p10Hits)
p10Hits.BY <- topTable(wtEbFit, coef = "devStageP10", number=Inf, sort.by="none", adjust.method="BY")
data <- cbind(p10Hits[, c("P.Value", "adj.P.Val")], "BY"=p10Hits.BY$adj.P.Val)

splom(data[1:100000, ], panel=function(x,y, ...) {
  panel.xyplot(x, y, ...)
  panel.abline(0,1, col="orange")
}) 

# BH adjusted P-values are larger than raw P-values, but not that much which implies their true importance
# which matches are previous observations

## Perform inference for some contrasts

(cont.matrix <- makeContrasts(P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - 
                                devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
res <- topTable(wtEbFitCont)

# plot top 4
top.genes <- res[1:4, ]
head(top.genes)
makeStripplot(subset(prepareData(top.genes$ID), gType == "wt"))

# adjust P-values globally & use a cut-of
cutoff <- 1e-04
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

# those that decline from P10 to Fourweeks, thus have a negative ttest
declining.genes <- rownames(wtResCont)[which(wtResCont[,"fourweeksVsP10"] < 0)]
makeStripplot(subset(prepareData(declining.genes[1:4]), gType == "wt"))

# no intersection 
intersect(top.genes, declining.genes[1:4])

# Rise from P10 to 4-weeks?
rising.genes <- rownames(wtResCont)[which(wtResCont[,"fourweeksVsP10"] > 0)]
makeStripplot(subset(prepareData(rising.genes[1:4]), gType == "wt"))

# interaction? no...
intersect(top.genes, rising.genes)
intersect(declining.genes, rising.genes)

# redo the analysis with a less treshold
cutoff <- 0.01
nHits <- 8
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

hits1 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] < 0)]
makeStripplot(subset(prepareData(hits1[1:nHits]), gType=="wt" ))

hits2 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] < 0)]
makeStripplot(subset(prepareData(hits2[1:nHits]), gType == "wt"))

hits3 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0)]
makeStripplot(subset(prepareData(hits3[1:nHits]), gType == "wt"))

hits4 <- rownames(prDat)[which(wtResCont[, "fourweeksVsP10"] > 0)]
makeStripplot(subset(prepareData(hits4[1:nHits]), gType == "wt"))

# indeed very useful
vennDiagram(wtResCont)

hits5 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] != 0 & wtResCont[, "fourweeksVsP10"] != 0)]
makeStripplot(subset(prepareData(hits5), gType == "wt"))

hits6 <- rownames(prDat)[which(wtResCont[, "P10VsP6"] > 0 & wtResCont[, "fourweeksVsP10"] < 0)]
makeStripplot(subset(prepareData(hits6), gType == "wt"))

## Take-home exercise
## make a new contrast matrix to enforce some change for previous stages
cutoff <- 0.01
nHits <- 8
colnames(wtFit)
(cont.matrix <- makeContrasts(P2VsE16 = devStageP2 - (Intercept), P6VsP2 = devStageP6 - devStageP2, P10VsP6 = devStageP10 - devStageP6, fourweeksVsP10 = devStage4_weeks - 
                                devStageP10, levels = wtDesMat))
wtFitCont <- contrasts.fit(wtFit, cont.matrix)
wtEbFitCont <- eBayes(wtFitCont)
wtResCont <- decideTests(wtEbFitCont, p.value = cutoff, method = "global")
summary(wtResCont)

# force some change both the way
hits7 <- rownames(prDat)[which( (wtResCont[, "P2VsE16"] != 0 & wtResCont[, "P6VsP2"] != 0) &  
                                  wtResCont[, "P10VsP6"] == 0 & wtResCont[, "fourweeksVsP10"] == 0)]
makeStripplot(subset(prepareData(hits7[4:6]), gType == "wt"))

head(wtResCont)
