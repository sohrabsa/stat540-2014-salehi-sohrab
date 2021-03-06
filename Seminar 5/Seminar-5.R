# Seminar 5
# Fitting and interpretting linear models (low volume)
# ref: http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar05_lowVolumeLinearModelling.html

library(lattice)
library(ggplot2)
library(car)

# load the data
dataDir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Second Semester/Stat 540/Lab/stat540_2014/examples/photoRec/data/"
prDat <- read.table(paste(dataDir, "GSE4051_data.tsv", sep = ""))
str(prDat, max.level = 0)
prDesc <- readRDS(paste(dataDir, "GSE4051_design.rds", sep=""))
str(prDesc)

(luckyGenes <- c("1419655_at","1438815_at"))


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

jDat <- prepareData(luckyGenes)

str(jDat)
head(jDat)
tail(jDat)


# the graphs match...
stripplot(gExp ~ devStage | gene, jDat,
          group = gType, jitter.data = TRUE,
          auto.key = TRUE, type = c('p', 'a'), grid = TRUE)

# the plotting function

makeStripplot <- function(a.data.frame, pch = 1, cex = 3) {
  # check that the input has the minimum requirements
  reqVariables <- c("gExp", "devStage", "gType", "gene")
  if (!all(sapply(reqVariables, function(x) !is.null(a.data.frame[[x]])))) {
    print("Error! The input data.frame doesn't contain required fields!")
    return(NULL)
  }

  (f <- stripplot(gExp ~ devStage | gene, a.data.frame,
            group = gType, jitter.data = TRUE,
            auto.key = TRUE, type = c('p', 'a'), grid = TRUE, pch = pch, cex = cex))
}

# fugly version...
makeStripplot(jDat, pch = 17, cex = 3)

makeStripplot(jDat)
str(newDat)
head(newDat)


## Do a two-sample t-test
targetProbesetName <- "1456341_a_at"
targetDevStage1 <- "P2"
targetDevStage2 <- "4_weeks"

levels(prDesc$devStage)
# assuming a common variance
ttDat <- prepareData(list(targetProbesetName))
t.test(gExp~devStage, subset(ttDat, devStage %in% c(targetDevStage1, targetDevStage2)), var.equal=T )
  
## Fit a linear model with a categorical covariate
targetProbesetName <- "1438786_a_at"
ttDat <- prepareData(list(targetProbesetName))
makeStripplot(ttDat)

# use lm to fit a linear model
the.fit <-lm(formula=gExp~devStage, subset(ttDat, gType == "wt"))
summary(the.fit)

# the Estimate for Intercept make sence, as its the predicted value of gExp when the independent variable
# (in this case devStage) is zero (as categorical variables, they're encoded in such a way where the first 
# one is encoded as zero and so on) and given the graph for the wt, which is above 8.5 for the wt at E16, its 
# sounds a reasonable.

## Perform inference for a contrast
# P2 and P10
coef(the.fit)
contMat <- matrix(c(0,1,0,-1,0), nr=1, nc=5)
(obsDiff <- contMat %*% (coef(the.fit)))

(sampMeans <- aggregate(gExp ~ devStage, ttDat, FUN = mean,
                        subset = gType == "wt"))
with(sampMeans, gExp[devStage == "P2"] - gExp[devStage == "P10"])

# get the variance covariance matrix
vcov(the.fit)

# check 
summary(the.fit)$coefficients[ , "Std. Error"]
sqrt(diag(vcov(the.fit)))

# calcualte the standard error for the statistic
(estSe <- contMat %*% vcov(the.fit) %*% t(contMat))
(testStat <- obsDiff/estSe)

# calculate a two-way p-value, given that the statistic  has a Student's distribution and is symetric
2 * pt(abs(testStat), df = df.residual(the.fit), lower.tail = FALSE)

## Fit a linear model with two categorical covariates
makeStripplot(oDat <- prepareData("1448690_at"))
str(oDat)


oFitBig <- lm(gExp~devStage*gType, dat=oDat)
summary(oFitBig)
summary(oFitBig)$coef

# now omit the interaction term
oFitSmall <- lm(gExp~devStage+gType, dat=oDat)
summary(oFitSmall)$coef

# compare the models
(a <- anova(oFitBig, oFitSmall))


# test another probeset, 1429225_at
makeStripplot(oDat <- prepareData("1429225_at"))
str(oDat)


oFitBig <- lm(gExp~devStage*gType, dat=oDat)
summary(oFitBig)
summary(oFitBig)$coef

# now omit the interaction term
oFitSmall <- lm(gExp~devStage+gType, dat=oDat)
summary(oFitSmall)$coef

# compare the models
(a <- anova(oFitBig, oFitSmall))


## Ideas for further work
# inputs a list of genes, conducts lm over them all and returns the results in a list
aggregateTesting <- function(gene.list) {
  theDat <- prepareData(gene.list)
  lapply(gene.list, function(x) {lm(gExp~devStage, data=theDat, subset= gene == x)})
}

k = aggregateTesting(list("1429225_at", "1448690_at"))

summary(k[[1]])

## corresponding to the lecture method, 
install.packages("car")
library(car)
# copy the data
agDat <- oDat

# make it quantitative
agDat$devStage <-recode(agDat$devStage, 
                   "'E16'=-2; 'P2'=2; 'P6'=6; 'P10'=10; '4_weeks'=28",
                   as.factor.result = FALSE)

agDat
# make a function to plot devStage
plotLinearModel <- function(the.data, formula.string) {
  ggplot(the.data, aes(x=devStage, y=gExp)) +
    geom_point() +
    stat_smooth(aes(x=devStage, y=gExp), formula=reformulate(formula.string, 'y'), method="lm")
    #stat_smooth(formula=y~x, method="lm")
}

plotLinearModel(agDat, "y~poly(x,3)")

## one method is to introduce a completely different covariate
# generate random values from a normal distribution and use them as age
agDat <- oDat
agDat$age <- rnorm(nrow(agDat), m=20, sd=5)

# fit a linear model

# fit a quadratic model
ggplot(agDat, aes(age, gExp)) +
  geom_point() +
  stat_smooth(formula=y~poly(x, 2), method="lm")


# dropping the 4-week data and analyze using quadratic model
agDat <- subset(agDat, devStage != "28")
plotLinearModel(agDat, "y~poly(x,2)")



##########################################
?pt
t = matrix(2, nrow = 20, ncol = 1)
t[10] <- 2

mu1 = 5
mu2 = 10
mu = matrix(c(mu1, mu2), nrow = 2, ncol = 1)

createDesignMatrix <- function(paramNum, len) {
  designMatrix <- matrix(0, nrow  = len, ncol = paramNum)
  k = lapply(seq(paramNum), function(x) { 
    startPoint <- (x-1) * len/paramNum + 1
    stopPoint <- startPoint - 1 + len/paramNum
    print(paste(startPoint, stopPoint, sep="_"))
    designMatrix[startPoint:stopPoint, x] <<- 1
    })
  designMatrix
}

createDesignMatrix(3, 30)
