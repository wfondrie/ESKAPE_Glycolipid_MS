library(ggdendro)
library(gridExtra)
library(grid)
library(xlsx)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
source("R/AnalyzeSpectra.R")
set.seed(1234)

features <- AnalyzeSpectra("data/large/", snr = 8)

# Heirarchical cluster of samples ----------------------------------------------
library("pvclust")
pv <- pvclust(t(features), method.hclust = "ward.D2",
              method.dist = "euclidean")
plot(pv, print.num=F)

# Dot Product ------------------------------------------------------------------
# Normalizes all vectors to 1
normFeatures <- aaply(features,1,function(x) x/sqrt(x %*% x))
dot <- normFeatures %*% t(normFeatures)
df <- as.data.frame(dot)

# reorders to match heirarchical clustering
df2 <- df
df2 <- df2[pv$hclust$order, pv$hclust$order]

# Add a column with sample name
df$names <- row.names(df)
df2$names <- row.names(df2)
write.csv(df2[nrow(df2):1, ],"results/largeHeatMap.csv")


# Making a heatmap -------------------------------------------------------------
# Reformats wide dataframe to long form for plotting with ggplot
dfMelt <- melt(df)
df2Melt <- melt(df2)

df2Melt$names <- factor(df2Melt$names, levels = levels(df2Melt$variable), ordered = T)

# The plotting functions!
p <- ggplot(df2Melt, aes(variable,names)) + 
  geom_tile(aes(fill=value), color = "black", size=0.3) +
  scale_fill_gradient2(name = "Spectrum\nSimilarity", 
                       low="white", mid = "red", high="black",
                       midpoint = 0.5) +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(color = "black", size = 6),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 12, color = "black"),
        legend.position = "left") +
  guides(fill = guide_colorbar(ticks = F )) +
  coord_fixed(ratio = 1)
p
ggsave("results/largeHeatmap.png",width=10,height=10,units="in")

dData <- dendro_data(as.dendrogram(pv$hclust))

theme_none <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.title.x = element_text(colour=NA),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  plot.margin = unit(c(0,0,0,0), "cm"),
  panel.margin = unit(c(0,0,0,0), "lines")
)


xDendro <- ggplot(segment(dData)) + 
  geom_segment(aes(x=x,y=y, yend=yend, xend=xend)) +
  theme(axis.title.x = element_blank()) +
  theme_none
xDendro

yDendro <- ggplot(segment(dData)) + 
  geom_segment(aes(x=x,y=y, yend=yend, xend=xend)) +
  theme(axis.title.x = element_blank()) +
  theme_none + coord_flip()
yDendro


pdf("results/fullHeatmap.pdf")

# Arrange dendrograms and spectra so they match, then write to file.
print(xDendro, vp=viewport(0.585, 0.2, x=0.49, y=0.9))
print(yDendro, vp=viewport(0.2, 0.6, x=0.87, y=0.535))
print(p, vp=viewport(0.8, 0.8, x=0.4, y=0.53))
dev.off()


