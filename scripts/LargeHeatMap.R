library(MALDIquant)
library(MALDIquantForeign)
library(rafalib)
library(ggdendro)
library(gridExtra)
library(grid)
library(xlsx)
library(plyr)
library(reshape2)
library(ggplot2)
library(scales)

#############################################################################
################## Preparing Mass Spectra ###################################
#############################################################################
## Any parts with "mypar(x,y)...mypar()" are just plotting mass spectra. I commented out most of these.


# Read in mzXML files #######################################################
file.list <- list.files("data/large/", pattern=".mzXML$",full.names=T)
spectra <- importMzXml(file.list)

# mypar(4,4)
# lapply(spectra,plot,xlim=c(400,2000))
# mypar()

# Normalize MS intensities ##################################################
spectra2 <- transformIntensity(spectra, method="sqrt")
mypar(4,4)
lapply(spectra2,plot,xlim=c(400,2000))
mypar()

spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=5)
# mypar(4,4)
# lapply(spectra3,plot,xlim=c(400,2000))
# mypar()

# mypar(4,4)
# lapply(spectra3, function(x) {
#   baseline <- estimateBaseline(x, method="SNIP",iterations=60)
#   plot(x)
#   lines(baseline,col="red",lwd=2)
# })
# mypar()

# Remove Baseline ################################################################
spectra4 <- removeBaseline(spectra3, method="SNIP", iterations=60)
# mypar(4,4)
# lapply(spectra4,plot)

spectra5 <- calibrateIntensity(spectra4,method="TIC")
# lapply(spectra5, plot)

# Align MS #######################################################################
spectra6 <- alignSpectra(spectra5,reference = spectra5[[1]])
# lapply(spectra6, plot)
# mypar()

# Pick Peaks #####################################################################
peaks <- detectPeaks(spectra6, SNR=2.2, halfWindowSize=5, method="SuperSmoother")
# mypar(4,4)

# for (i in 1:length(spectra6)) {
#   plot(spectra6[[i]])
#   points(peaks[[i]], col="red", pch=4)
# }


peaks <- binPeaks(peaks, tolerance=0.002) # This bins peaks that are very close together

# Preparing Spectra for dot product ##############################################
# This small section extracts the sample file name and adjusts it to be the sample name.
sample <- sapply(spectra6, function(x) metaData(x)$file)
sample <- gsub(".*\\\\", "", sample)
sample <- gsub(".mzXML$","",sample)
sample <- factor(sample)

# Retrieves identified peaks and intensities as a matrix (row = sample, col = m/z, value = intensity)
features <- intensityMatrix(peaks, spectra6)
rownames(features) <- sample


# Heirarchical cluster of samples
library("pvclust")
pv <- pvclust(t(features), method.hclust = "ward.D2",
              method.dist = "euclidean")
mypar()
plot(pv, print.num=F)



# Dot Product ########################################################################################
# Normalizes all vectors to 1
normFeatures <- aaply(features,1,function(x) x/sqrt(x %*% x))
dot <- normFeatures %*% t(normFeatures)
df <- as.data.frame(dot)

# reorders to match heirarchical clustering
df2 <- df
df2 <- df2[pv$hclust$order, pv$hclust$order]

# The part below would have been used to match a file name to key
### Match file names to Organism name ### Uncomment for once names are updated!
# keyRough <- read.xlsx("data/large/ESKAPE library spectra DB.xlsx",1, stringsAsFactors =F)
# keyLeft <- keyRough[,1:8]
# keyRight <- keyRough[,10:17]
# 
# names(keyLeft)[2:8] <- keyLeft[1, 2:8]
# names(keyRight) <- names(keyLeft)
# key <- merge(keyLeft,keyRight, all = T)
# key <- key[!is.na(key[,1]) & !is.na(key[,5]),]
# 
# newName <- lapply(row.names(df2), function(x) {
#   match <- grepl(paste0("^",x,"$"),key[,"Heat map number"], ignore.case = T)
#   r <- key[match,"Heat map name"] 
#   r
# })
# 
# row.names(df2) <- names(df2) <- newName

#########################################

# Add a column with sample name
df$names <- row.names(df)
df2$names <- row.names(df2)


# Making a heatmap ##############################################################################
# Reformats wide dataframe to long form for plotting with ggplot
dfMelt <- melt(df)
df2Melt <- melt(df2)

df2Melt$names <- factor(df2Melt$names, levels = levels(df2Melt$variable), ordered = T)

# The plotting functions!
p <- ggplot(df2Melt, aes(variable,names)) + 
  geom_tile(aes(fill=value), color = "black", size=0.3) +
  scale_fill_continuous(name = "Spectrum\nSimilarity", low=("white"), high="darkgreen") +
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
ggsave("results/largeHeatmap_green_ordered.png",width=10,height=10,units="in")

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

# Fin


