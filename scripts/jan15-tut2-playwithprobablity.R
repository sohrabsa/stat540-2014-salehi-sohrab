# second lab

min.y = 0
max.y = 20
num.of.samples = 1000
y = seq(min.y, max.y, length = num.of.samples)

plot(c(0,20), c(0, .5), xlab="y", ylab="f(y)", main="The shaped based realization of Gamma pdf", type="n")

# do it with sapply with more clearly defined variables
shapes=c(1,2,3,5,9)
scales=c(2,2,2,1,.5)
cols=c("red", "green", "blue", "magenta", "black")
plots=sapply(1:5, function(x)  { lines(y, dgamma(y, shape=shapes[x], scale=scales[x]), col=cols[x]  ) } )

# do this with mapply and vage variables
plot(c(0,20), c(0, .5), xlab="y", ylab="f(y)", main="The shaped based realization of Gamma pdf", type="n")
f <- function(col, ...) {
  lines(y, dgamma(y, ...), col = col)
}

plot.status <- mapply(f, shape=shapes, scale=scales, col=cols)