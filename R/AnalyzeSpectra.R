## Description
# Takes input mzXML spectra and runs them through a MALDIQuant pipeline for processing
# MALDI spectra. Outputs a peak list for each spectra.

## Input
# Path - a folder containing mzXML spectra
# snr - acceptable SNR ratio for peaks

## Output
# Features - a peak list of detected from the paraeters providedm

AnalyzeSpectra <- function(path,snr){
  
  library(MALDIquant)
  library(MALDIquantForeign)
  library(rafalib)
  
  # Read in mzXML files #######################################################
  file.list <- list.files(path, pattern=".mzXML$",full.names=T)
  spectra <- importMzXml(file.list)

  # Normalize MS intensities ##################################################
  spectra2 <- transformIntensity(spectra, method="sqrt")
  spectra3 <- smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=5)

  # Remove Baseline ############################################################
  spectra4 <- removeBaseline(spectra3, method="SNIP", iterations=60)
  
  spectra5 <- calibrateIntensity(spectra4,method="TIC")
  
  # Align MS ###################################################################
  spectra6 <- alignSpectra(spectra5,reference = spectra5[[1]])
  
  # Pick Peaks #################################################################
  peaks <- detectPeaks(spectra6, SNR=snr, halfWindowSize=5, method="MAD")
  peaks <- binPeaks(peaks, tolerance=0.002) # This bins peaks that are very close together
  
  # Preparing Spectra for dot product ##########################################
  # This small section extracts the sample file name and adjusts it to be the sample name.
  sample <- sapply(spectra6, function(x) metaData(x)$file)
  sample <- gsub(".*\\\\", "", sample)
  sample <- gsub(".mzXML$","",sample)
  sample <- factor(sample)
  
  # Retrieves identified peaks and intensities as a matrix (row = sample, col = m/z, value = intensity)
  features <- intensityMatrix(peaks, spectra6)
  rownames(features) <- sample
  
  return(features)
}