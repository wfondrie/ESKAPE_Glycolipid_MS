set.seed(1)

library("plyr")

a <- seq(from=0, to=3,by=0.1)
b <- rnorm(length(a))

t <- rbind(a,b)

ints <- seq(0,3,by=1)

bins <- cut(a,breaks=ints)

test<- tapply(t[2,], bins, sum)

test<- aaply(t, 1, function(x) tapply(x,bins,sum))
