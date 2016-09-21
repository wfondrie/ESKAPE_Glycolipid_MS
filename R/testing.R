library(MALDIquant)
library(MALDIquantForeign)

file.list <- list.files("data/", pattern=".mzXML$")
spectra <- importMzXml(paste0("Data/",file.list))
#spectra <- trim(spectra, c(100,2000))

idx <- sample(length(spectra), size=2)

spectra2 <- transformIntensity(spectra, method="sqrt")

spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=10)

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

spectra5 <- calibrateIntensity(spectra4,method="TIC")
spectra6 <- alignSpectra(spectra5)

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

peaks <- detectPeaks(spectra6, SNR=2.5, halfWindowSize=10, method="SuperSmoother")
a <- 16
plot(spectra6[[a]])
points(peaks[[a]], col="red", pch=4)

peaks <- binPeaks(peaks)

points(peaks[[a]], col="blue", pch=3)

sample <- sapply(spectra6, function(x) metaData(x)$file)
sample <- gsub(".*\\\\", "", sample)
sample <- gsub(".mzXML$","",sample)
sample <- factor(sample)

features <- intensityMatrix(peaks, spectra6)
rownames(features) <- sample
