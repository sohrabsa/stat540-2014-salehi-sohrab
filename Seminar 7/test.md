Title
========================================================

This is an R Markdown document. Markdown is a simple formatting syntax for authoring web pages (click the **Help** toolbar button for more details on using R Markdown).

When you click the **Knit HTML** button a web page will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

$$ Y \sim NB(\mu, \varphi) \text{ with mean } \mu \text{ and dispersion } \varphi$$

$$ f(y; \mu, \phi) = P(Y=y) = \frac{\Gamma(y+\phi^{-1})}{\Gamma(\phi^{-1})\Gamma(y+1)}(\frac{1}{1+\mu \phi})^{\phi^{-1}}(\frac{\mu}{\phi^{-1} + \mu})^y  $$



```r
# Here we're estimating different types of dispersion. According to
# [here](https://stat.ethz.ch/pipermail/bioconductor/2014-January/056975.html),
# this has little to do with the final results.  We can just use y <-
# estimateDisp(y,design) and it will do these steps itself altogether a.
# assuming all genes/tags share a common dispersion

```

# According to [this,](https://stat.ethz.ch/pipermail/bioconductor/2014-January/056975.html) the edgeR pipline has been updated, now you can just use:


# Limma and Voom
## [ref 2: ](http://www.statsci.org/smyth/pubs/VoomPreprint.pdf) explains how a not-count specific method
such as limma could be adapted to count data. Main idea is based on modeling mean-varience relationship.



