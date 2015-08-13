 library("plyr")
 library("reshape2")
 library("ggplot2")
 library("scales")

# masses <- as.numeric(colnames(features))
# ints <- seq(400,2400,by=1)
# 
# bins <- cut(masses,breaks=ints)
# 
# binnedFeatures <- aaply(features, 1, function(x) tapply(x,bins,sum))
# binnedFeatures[is.na(binnedFeatures)] <- 0
# 
# dot <- binnedFeatures %*% t(binnedFeatures)
# 
# normDot <- (dot-min(dot))/(max(dot)-min(binnedFeatures))
 
# altFeatures <- intensityMatrix(peaks)
# rownames(altFeatures) <- rownames(features)
# altFeatures[is.na(altFeatures)] <- 0
# 
# dot <- altFeatures %*% t(altFeatures)
 
normFeatures <- aaply(features,1,function(x) x/sqrt(x %*% x))

dot <- normFeatures %*% t(normFeatures)

normDot <- dot
#normDot <- aaply(dot,1,rescale)

df <- as.data.frame(normDot)
df$names <- row.names(df)

dfMelt <- melt(df)

p <- ggplot(dfMelt, aes(variable,names)) + geom_tile(aes(fill=value)) +
  scale_fill_continuous(low=("black"), high="red")
p
