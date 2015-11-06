library(MALDIquant)
library(MALDIquantForeign)
library(rafalib)

file.list <- list.files("data/small/", pattern=".mzXML$",full.names=T)
spectra <- importMzXml(file.list)
#spectra <- trim(spectra)

#idx <- sample(length(spectra), size=2)
mypar(4,4)
#lapply(spectra,plot,xlim=c(400,2000))
mypar()

spectra2 <- transformIntensity(spectra, method="sqrt")
mypar(4,4)
#lapply(spectra2,plot,xlim=c(400,2000))
mypar()

spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=5)
mypar(4,4)
#lapply(spectra3,plot,xlim=c(400,2000))
mypar()

mypar(4,4)
# lapply(spectra3, function(x) {
#   baseline <- estimateBaseline(x, method="SNIP",iterations=60)
#   plot(x)
#   lines(baseline,col="red",lwd=2)
# })
mypar()

# ## define iteration steps: 25, 50, ..., 100
# iterations <- seq(from=25, to=100, by=25)
# ## define different colors for each step
# col <- rainbow(length(iterations))
# plot(spectra3[[1]])
# ## draw different baseline estimates
# for (i in seq(along=iterations)) {
#   baseline <- estimateBaseline(spectra3[[1]], method="SNIP",
#                                iterations=iterations[i])
#   lines(baseline, col=col[i], lwd=2)
# }
# legend("topright", legend=iterations, col=col, lwd=1)

spectra4 <- removeBaseline(spectra3, method="SNIP", iterations=60)
mypar(4,4)
#lapply(spectra4,plot)

spectra5 <- calibrateIntensity(spectra4,method="TIC")
#lapply(spectra5, plot)

spectra6 <- alignSpectra(spectra5)
#lapply(spectra6, plot)
mypar()

# ## define snrs steps: 1, 1.5, ... 2.5
# snrs <- seq(from=1, to=2.5, by=0.5)
# ## define different colors for each step
# col <- rainbow(length(snrs))
# ## estimate noise
# noise <- estimateNoise(spectra6[[1]],
#                        method="SuperSmoother")
# plot(spectra6[[1]])
# for (i in seq(along=snrs)) {
#   lines(noise[, "mass"],
#         noise[, "intensity"]*snrs[i],
#         col=col[i], lwd=2)
# }
# legend("topright", legend=snrs, col=col, lwd=1)

peaks <- detectPeaks(spectra6, SNR=2.2, halfWindowSize=5, method="SuperSmoother")
mypar(4,4)

# for (i in 1:length(spectra6)) {
#   plot(spectra6[[i]])
#   points(peaks[[i]], col="red", pch=4)
# }


peaks <- binPeaks(peaks, tolerance=0.002)

sample <- sapply(spectra6, function(x) metaData(x)$file)
sample <- gsub(".*\\\\", "", sample)
sample <- gsub(".mzXML$","",sample)
sample <- factor(sample)

features <- intensityMatrix(peaks, spectra6)
rownames(features) <- sample


library("pvclust")
pv <- pvclust(t(features), method.hclust = "ward.D2",
              method.dist = "euclidean")
plot(pv, print.num=F)


library("plyr")
library("reshape2")
library("ggplot2")
library("scales")

normFeatures <- aaply(features,1,function(x) x/sqrt(x %*% x))

dot <- normFeatures %*% t(normFeatures)

df <- as.data.frame(dot)
df$names <- row.names(df)
labels <- gsub(" \\(colistin-resistant\\)","\\*", df$names)

dfMelt <- melt(df)
dfMelt$names <- factor(dfMelt$names, levels = levels(dfMelt$variable))

p <- ggplot(dfMelt, aes(variable,factor(names))) + 
  geom_tile(aes(fill=value), color = "black", size = 1) +
  scale_fill_continuous(name = "Spectrum\nSimilarity", low=("white"), high="darkgreen") +
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.text = element_text(color="black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        text = element_text(size = 18, color = "black")) +
  guides(fill = guide_colorbar(ticks = F )) +
  scale_x_discrete(labels = labels) +
  scale_y_discrete(labels = labels) +
  coord_fixed(ratio = 1)
p
ggsave("results/smallHeatmap_green.png",width=10,height=10,units="in")