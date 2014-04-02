Seminar 9
========================================================

Ref1: http://www.ugrad.stat.ubc.ca/~stat540/seminars/seminar09_clustering-pca.html

*Excerpt from ref*
Clustering on the photoreceptor data set.


```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(plyr)
library(gplots)
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(RColorBrewer)
library(cluster)
library(pvclust)
library(xtable)
library(limma)
```



```r
source("style.R")
```


# load the data

```r
dataDir <- "/Users/sohrab/Me/Apply/Canada Apply/Courses/Second Semester/Stat 540/Lab/stat540_2014/examples/photoRec/data/"
prDat <- read.table(paste(dataDir, "GSE4051_data.tsv", sep = ""), header = T, 
    row.names = 1)
str(prDat, max.level = 0)
```

```
## 'data.frame':	29949 obs. of  39 variables:
```

```r
prDesc <- readRDS(paste(dataDir, "GSE4051_design.rds", sep = ""))
str(prDesc)
```

```
## 'data.frame':	39 obs. of  4 variables:
##  $ sidChar : chr  "Sample_20" "Sample_21" "Sample_22" "Sample_23" ...
##  $ sidNum  : num  20 21 22 23 16 17 6 24 25 26 ...
##  $ devStage: Factor w/ 5 levels "E16","P2","P6",..: 1 1 1 1 1 1 1 2 2 2 ...
##  $ gType   : Factor w/ 2 levels "wt","NrlKO": 1 1 1 1 2 2 2 1 1 1 ...
```


Rescale the data (the rows, i.e., genes across samples) to aid visualization

```r
sprDat <- t(scale(t(prDat)))
str(sprDat, max.level = 0, give.attr = FALSE)
```

```
##  num [1:29949, 1:39] 0.0838 0.1758 0.7797 -0.3196 0.8358 ...
```

```r

# check the first 10...
round(data.frame(avgBefore = rowMeans(head(prDat)), avgAfter = rowMeans(head(sprDat)), 
    varBefore = apply(head(prDat), 1, var), varAfter = apply(head(sprDat), 1, 
        var)), 2)
```

```
##              avgBefore avgAfter varBefore varAfter
## 1415670_at        7.22        0      0.02        1
## 1415671_at        9.37        0      0.35        1
## 1415672_at        9.70        0      0.15        1
## 1415673_at        8.42        0      0.03        1
## 1415674_a_at      8.47        0      0.02        1
## 1415675_at        9.67        0      0.03        1
```


It worked, the rows have mean = 0 and variance 1.

# Sample Clustering

## Hierarchical Clustering
### Standards: Eucleadian Distance, average linkage kernel
#### Note
`dist` calculates distance between rows

```r
# create the distance matrix,
pr.dis <- dist(t(sprDat), method = "euclidean")

# interaction is a cool method :)
prDesc$grp <- with(prDesc, interaction(gType, devStage))
summary(prDesc$grp)
```

```
##        wt.E16     NrlKO.E16         wt.P2      NrlKO.P2         wt.P6 
##             4             3             4             4             4 
##      NrlKO.P6        wt.P10     NrlKO.P10    wt.4_weeks NrlKO.4_weeks 
##             4             4             4             4             4
```



```r
# do hierarchical clustering using different kernels
linkage.types <- c("single", "complete", "average", "ward")
pr.hc <- lapply(linkage.types, function(x) hclust(pr.dis, method = x))

# plot the results
op <- par(mar = c(0, 4, 4, 2), mfrow = c(2, 2))

lapply(seq(pr.hc), function(x) plot(pr.hc[[x]], labels = F, main = linkage.types[[x]], 
    xlab = ""))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6.png) 

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
```



