library(plyr)
library(reshape2)
library(ggplot2)
library(scales)
source("R/AnalyzeSpectra.R")
set.seed(1234)

# Imports spectra, normalizes and performs peak picking. Check 
# <functions/AnalyzeSpectra.R> for details

features <- AnalyzeSpectra("data/small", snr = 8)

library("pvclust")
pv <- pvclust(t(features), method.hclust = "ward.D2",
              method.dist = "euclidean")
plot(pv, print.num=F)

# Dot Product ------------------------------------------------------------------
# Normalizes all vectors to 1
normFeatures <- aaply(features,1,function(x) x/sqrt(x %*% x))

dot <- normFeatures %*% t(normFeatures)

# Correct Names ----------------------------------------------------------------
# Reorder names to ESKAPE order
ord <- unique(c(grep("faecium", rownames(dot), value = T, ignore.case = T),
                grep("aureus", rownames(dot), value = T, ignore.case = T),
                grep("klebsiella", rownames(dot), value = T, ignore.case = T),
                grep("pneumoniae", rownames(dot), value = T, ignore.case = T),
                grep("baumannii", rownames(dot), value = T, ignore.case = T),
                grep("aeruginosa", rownames(dot), value = T, ignore.case = T),
                grep("enterobacter", rownames(dot), value = T, ignore.case = T)))

dot <- dot[ord, ord]

# needs to be flipped around to match heatmap
dfwrite <- as.data.frame(dot[nrow(dot):1, ])
write.csv(dfwrite, "results/smallHeatMap.csv")


# make names good for heatmap
df <- as.data.frame(dot)
df$names <- factor(row.names(df), levels = ord)
labels <- gsub(" \\(colistin-resistant\\)","\\*", df$names)


# Making a heatmap -------------------------------------------------------------
# Reformats wide dataframe to long form for plotting with ggplot
dfMelt <- melt(df)
dfMelt$names <- factor(dfMelt$names)
dfMelt$variable <- factor(dfMelt$variable)
dfMelt$names <- factor(dfMelt$names, levels = ord)
dfMelt$variable <- factor(dfMelt$variable, levels = ord)

# The Plotting
p <- ggplot(dfMelt, aes(variable,names)) + 
  geom_tile(aes(fill=value), color = "black", size = 1) +
  scale_fill_gradient2(name = "Spectrum\nSimilarity", 
                       low="white", mid = "red", high="black",
                       midpoint = 0.5) +
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
ggsave("results/smallHeatmap.png",width=10,height=10,units="in")
ggsave("results/smallHeatmap.pdf", width=10,height=10,units="in", useDingbats = F)
