# Two group comparison
# ref: http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar04_compileNotebook-dataAggregation-twoGroupComparison.html

library(lattice)

COMMON.SEED = 987

dataDir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Second Semester/Stat 540/Lab/stat540_2014/examples/photoRec/data/"

prDat <- read.table(paste(dataDir, "GSE4051_data.tsv", sep=""))
prDes <- readRDS(paste(dataDir, "GSE4051_design.rds", sep=""))

str(prDat)
str(prDat, max.level = 0)
str(prDes)

set.seed(COMMON.SEED)

theGene <- sample(1:nrow(prDat), 1)
theGene
pDat <- data.frame(prDes, gExp = unlist(prDat[theGene, ]))
str(pDat)

head(pDat)
aggregate(gExp~gType, pDat, FUN = mean)
stripplot(gType~gExp, pDat)


theHtest <- t.test(gExp~gType, pDat)
str(theHtest)

theHtest$data.name
head(pDat)

# doesn't like the formula interface (has no s3 ...)
ks.test(gExp~gType, pDat)

# use the hard way
ks.test(pDat$gExp[pDat$gType == "wt"], pDat$gExp[pDat$gType !="wt"], pDat)



# Aggregation -------------------------------------------------------------

kDat <- readRDS(paste(dataDir, "GSE4051_MINI.rds", sep=""))
kMat <- as.matrix(kDat[c('crabHammer', 'eggBomb', 'poisonFang')])

head(kDat)
str(kMat)
median(kMat[, 'crabHammer'])
median(kMat[, 'eggBomb'])
median(kMat[, 'poisonFang'])

apply(kMat, 2, median)
apply(kMat, 2, quantile, probs = .99)

head(kMat)
# find min Sample for each gene
apply(kMat, 2, min)
which(kMat[, 'crabHammer'] == 8.214)

# fine min gene Expression for each sample
apply(kMat, 1, min)
colnames(kMat)
colnames(kMat)[apply(kMat, 1, which.min)]

colSums(kMat)
rowSums(kMat)

# sanity check?
all.equal(rowSums(kMat), apply(kMat, 1, sum))


# Aggregation on groups of observations -----------------------------------

aggregate(eggBomb~devStage, kDat, FUN=mean)
aggregate(eggBomb~devStage * gType, kDat, FUN=mean)

aggregate(eggBomb~devStage * gType, kDat, FUN=range)


# Two sample test, on a handful of genes -----------------------------------
keepGenes <- c("1431708_a_at", "1424336_at", "1454696_at", "1416119_at", "1432141_x_at", "1429226_at" )
# rownames(prDat) %in% keepGenes # Don't run this line, rowNames of course tries to print out every row
miniDat <- subset(prDat, rownames(prDat) %in% keepGenes)

head(miniDat)
la.transpose = t(miniDat) # transpose returns a matrix

la.transpose

# levels are the possible values that the factor could be set to, and factors are like categories
miniDat <- data.frame(gExp = as.vector(la.transpose),
                      gene = factor(rep(rownames(miniDat), each = ncol(miniDat)),
                                    levels = keepGenes))

# add the genotypes etc
# this is actually an amazing default tbehaviour by R
# One thing to note, if you supply a list, instead of a 
# vector, it will append each element of the list as a separate column
# while vectors would be repeated (recycled)
miniDat <- data.frame(prDes, miniDat)

# The help entry on stripplot is quite instructive, e.g., 
# To interpret y1 + y2 as a sum, one can either set 
#     allow.multiple=FALSE or use I(y1+y2).
stripplot(gType ~ gExp | gene, miniDat,
          scales = list(x = list(relation = "same")),
          group = gType, auto.key = F)



# Now multiple aggregation, into the plyr ------------------------------------------------

someData <- droplevels(subset(miniDat, gene == keepGenes[1]))
t.test(gExp ~ gType, someData)

library(plyr)
# data.frame implies functions starting with d in plyr
head(miniDat)

# output is discarded, dumped onto the screen
d_ply(miniDat, ~gene, function(x) { t.test (gExp~gType, x)}, .print=T)

results <- dlply(miniDat, ~gene, function(x) { t.test (gExp~gType, x)})
names(results)

# dl -> data.frame to list
# dd -> data.frame to data.frame

# for each gene, apply t.test
ttRes <- ddply( miniDat,  ~ gene, function(y) {
  yy <- t.test(gExp ~ gType, y)
  round(c(tStat = yy$statistic, est = yy$estimate), 3)
} )
ttRes


# take home ideas

# report p.values for t.test, wil.cox test and K.S test
ttRes <- ddply( miniDat,  ~ gene, function(y) {
  g1 = y$gExp[y$gType == "wt"] 
  g2 = y$gExp[y$gType != "wt"]
  yy1 <- t.test(gExp ~ gType, y)
  yy2 <- ks.test(g1, g2, y)
  yy3 <- wilcox.test(gExp ~ gType, y)
  round(c(t.pvalue = yy1$p.value, ks.pvalue = yy2$p.value, wilcox.pvalue = yy3$p.value), 3)
} )
ttRes

## more genes, lets sample a 100
# reload the data
prDat <- read.table(paste(dataDir, "GSE4051_data.tsv", sep=""))
prDes <- readRDS(paste(dataDir, "GSE4051_design.rds", sep=""))


sample.indexes <- sample(nrow(prDat), 100)
newData <- prDat[sample.indexes, ]

# reshape the data (remember to name the columns!)
newData <- data.frame(gExp = as.vector(t(newData)),
                      gene = factor(rep(rownames(newData), each=ncol(newData)), 
                             levels = rownames(newData))
                      )

newData <- data.frame(prDes, newData)

tests <- ddply(newData, ~gene, function(y) {
  g1 = y$gExp[y$gType == "wt"] 
  g2 = y$gExp[y$gType != "wt"]
  yy1 <- t.test(gExp ~ gType, y)
  yy2 <- ks.test(g1, g2, y)
  yy3 <- wilcox.test(gExp ~ gType, y)
  round(c(t.pvalue = yy1$p.value, ks.pvalue = yy2$p.value, wilcox.pvalue = yy3$p.value), 3)
})

tt <- data.matrix(tests[, -1])
rownames(tt) = tests$gene

# plot the pvalues etc
library(ggplot2)
require(gridExtra)
# ref: http://www.cookbook-r.com/Graphs/Scatterplots_(ggplot2)/
p1 <- ggplot(tests, aes(x=ks.pvalue, y=t.pvalue)) +
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line 
                          #  (by default includes 95% confidence region)

p2 <- ggplot(tests, aes(x=ks.pvalue, y=wilcox.pvalue)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm)

p3 <- ggplot(tests, aes(x=wilcox.pvalue, y=t.pvalue)) +
  geom_point(shape=1) +    
  geom_smooth(method=lm)

grid.arrange(p1, p2, p3, ncol=3)

# tranform into true/false matrix
TRESHOLD = .05
ttt <- aaply(tt, 1, function(y) {
  y < TRESHOLD
})
head(ttt)

# for each gene, how many tests indicate its significant?
rowSums(ttt)

# how many genes are hits by only two  test?
sum( aaply(ttt, 1, function(y) { ifelse(sum(y) == 2 ,1,0)} ))

# how many genes are hits by no test?
sum( aaply(ttt, 1, function(y) { ifelse(sum(y) == 0 ,1,0)} ))

## timing study

# put it all into a function
DEAnalyze <- function(sample.size = 10) {
  
  # get the samples
  sample.indexes <- sample(nrow(prDat), 100)
  newData <- prDat[sample.indexes, ]
  
  # reshape the data (remember to name the columns!)
  newData <- data.frame(gExp = as.vector(t(newData)),
                        gene = factor(rep(rownames(newData), each=ncol(newData)), 
                                      levels = rownames(newData))
  )
  
  newData <- data.frame(prDes, newData)
  
  ddply(newData, ~gene, function(y) {
    g1 = y$gExp[y$gType == "wt"] 
    g2 = y$gExp[y$gType != "wt"]
    yy1 <- t.test(gExp ~ gType, y)
    yy2 <- ks.test(g1, g2, y)
    yy3 <- wilcox.test(gExp ~ gType, y)
    round(c(t.pvalue = yy1$p.value, ks.pvalue = yy2$p.value, wilcox.pvalue = yy3$p.value), 3)
  })
}

# time it for different sample sizes
sampleSizes = c(10, 50, 100, 500)
de.times <- adply(sampleSizes, 1, function(y) {system.time(DEAnalyze(y)) })
de.times$sampleSize = sampleSizes
de.times
?adply
