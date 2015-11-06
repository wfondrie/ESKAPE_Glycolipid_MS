library(caret)
library(ggplot2)
library(reshape2)
library(plyr)
library(MALDIquant)
library(MALDIquantForeign)
library(rafalib)
library(pROC)
set.seed(1234)

ddP <- 20 # Number of peaks to keep from each spectra

#### Read in Files #############
file.list <- list.files("data/ROC/lib/", pattern=".mzML$",full.names=T)
spectra <- importMzMl(file.list)

mypar(7,5)
lapply(spectra,plot,xlim=c(400,2000))
mypar()

spectra2 <- transformIntensity(spectra, method="sqrt")

spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=5)

spectra4 <- removeBaseline(spectra3, method="SNIP", iterations=60)

spectra5 <- calibrateIntensity(spectra4,method="TIC")

spectra6 <- alignSpectra(spectra5,reference = spectra5[[1]])


peaks <- detectPeaks(spectra6, SNR=3.5, halfWindowSize=5, method="SuperSmoother")

peaks <- binPeaks(peaks, tolerance=0.002)

peaks <- filterPeaks(peaks, minNumber = 2)

mypar(7,5)
for (i in 1:length(spectra6)) {
  plot(spectra6[[i]])
  points(peaks[[i]], col="red", pch=4)
}
mypar()

############################################################

sample <- sapply(spectra6, function(x) metaData(x)$file)
sample <- gsub(".*\\\\", "", sample)
sample <- gsub(" [0-9]+.mzML$","",sample)
sample <- gsub(" ","_", sample)
sample[sample != "Pa"] <- "Not.Pa" 

featureMatrix <- intensityMatrix(peaks, spectra6)
featureMatrix <- aaply(featureMatrix, 1, function(x){
  ord <- order(x, decreasing = T)
  if(length(ord) > ddP){
    ord <- ord[(ddP + 1):length(ord)]
    x[ord] <- 0
  }
  x
})
row.names(featureMatrix) <- sample
featureMatrix <- featureMatrix[,colSums(featureMatrix)!=0]

features <- as.data.frame(featureMatrix)
features$sample <- sample
#features <- features[features$sample != "Ecf",]
features$sample <- factor(features$sample, levels = c("Pa","Not.Pa"))



########### Machine Learning #################
trainIdx <- createDataPartition(features$sample, p = 0.5, list = F, times = 1)
training <- features[trainIdx,]
test <- features[-trainIdx,]

fitCtrl <- trainControl(method = "LOOCV")
fit <- train(sample ~ ., data = training, method = "rf",
             trControl = fitCtrl)
t <- predict(fit, newdata=test)

postResample(t,test$sample)

confusionMatrix(t, test$sample)


########## Adding in silico Samples ################
notPa <- features[features$sample == "Not.Pa", names(features) != "sample"]
Pa <- features[features$sample == "Pa", names(features) != "sample"]


avg <- colMeans(notPa)
sd <- aaply(notPa,2,function(x){
  x <- as.numeric(unname(unlist(x)))
  x <- x[x > 0]
  x <- sd(x)
  x
})
sd[is.na(sd)] <- 0
statNotPa <- data.frame(avg = avg, sd = sd)

avg <- colMeans(Pa)
sd <- aaply(Pa,2,function(x){
  x <- as.numeric(unname(unlist(x)))
  x <- x[x > 0]
  x <- sd(x)
  x
})
sd[is.na(sd)] <- 0
statPa <- data.frame(avg = avg, sd = sd)

############### Making New Spectra #######################################
generateSpectra <- function(realLibrary, libraryStats, n, numPeaks) {
  
  baseSpectra <- sample(1:nrow(realLibrary), n, replace = T)
  
  newSpectra <- adply(baseSpectra, 1, function(x) {
    spec <- realLibrary[x,]
    spec <- adply(spec,1, function(y) {
      t <- y > 0
      y[t] <- y[t] + rnorm(length(y[t]), mean = 0, sd= libraryStats$sd[t]*5)
      y[t] <- rbinom(length(y[t]), 1, 0.50) * y[t]
      y[!t] <- abs(rbinom(length(y[!t]), 1, 0.75) * rnorm(length(y[!t]), 
                                                      mean = 0.01, 
                                                      sd = 0.01))#mean(libraryStats$sd))
      
    })
    spec
  })
  
  newSpectra
}

newNotPa <- generateSpectra(notPa, statNotPa,20, ddP)
newPa <- generateSpectra(Pa, statPa, 20, ddP)
newNotPa$sample <- rep("Not.Pa", nrow(newNotPa))
newPa$ sample <- rep("Pa", nrow(newPa))
newSpec <- merge(newNotPa, newPa, all = T)
newSpec$inSilico <- rep(T, nrow(newSpec))

features$inSilico <- rep(F, nrow(features))

features2 <- merge(features, newSpec, all = T)


newNotPa$num <- seq_along(newNotPa[,1])
newSpec <- melt(newNotPa[,sapply(newNotPa, is.numeric)], id.vars = "num")

pSpec <- ggplot(newSpec, aes(x=variable, ymax=value, ymin=0)) + 
  geom_linerange() +
  facet_wrap(~num, ncol = 4)
#pSpec


newPa$num <- seq_along(newPa[,1])
newSpec2 <- melt(newPa[,sapply(newPa,is.numeric)], id.vars = "num")

pSpec2 <- ggplot(newSpec2, aes(x=variable, ymax=value, ymin=0)) + 
  geom_linerange() +
  facet_wrap(~num, ncol = 4)
#pSpec2

################## ML with in silico library ##############################

featuresM <- features2[, names(features2) != "inSilico"]
trainIdx <- createDataPartition(featuresM$sample, p = 0.5, list = F, times = 1)
training <- featuresM[trainIdx,]
test <- featuresM[-trainIdx,]

fitCtrl <- trainControl(method = "LOOCV")
fit <- train(sample ~ ., data = training, method = "rf",
             trControl = fitCtrl)
prob <- predict(fit,newdata=test, type = "prob")
t$prob <- prob[,"Pa"]
t$class <- predict(fit, newdata=test)

ROC <- roc(response = t$class,
           predictor = t$prob,
           levels = rev(levels(t$class)))

confusionMatrix(t$class, t$class, positive = "Pa")


