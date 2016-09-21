library(MALDIquant)
library(MALDIquantForeign)
file <- "E. faecium FN-45 Isolate-29.1 from JKJ 25C.mzXML"
test <- importMzXml(paste0("Data/",file))
spec <- test[[1]]
plot(spec)

spec2 <- transformIntensity(spec,method="sqrt")
plot(spec2)

spec3 <- smoothIntensity(spec2,method="SavitzkyGolay", halfWindowSize=10)
plot(spec3)

baseline <- estimateBaseline(spec3, method="SNIP", iterations=100)
plot(spec3)
lines(baseline, col="red", lwd=2)

spec4 <- removeBaseline(spec3, method="SNIP", iterations=100)
plot(spec4)

#spec5 <- alignSpectra(spec4, halfWindowSize=20, SNR-2, tolerance=0.02, warpingMethod="lowess")
#plot(spec5)

noise <- estimateNoise(spec4)
plot(spec4)
lines(noise, col="red")
lines(noise[,1],noise[,2]*2,col="blue")
lines(noise[,1],noise[,2]*3,col="green")
lines(noise[,1],noise[,2]*4,col="yellow")

peaks <- detectPeaks(spec4, method="MAD", halfWindowSize=20, SNR=5)
plot(spec4)
points(peaks, col="red", pch=4)

peaks <- binPeaks(peaks, tolerance=0.002)
plot(spec4)
points(peaks, col="red", pch=4)

