################################################################################
# PARAMETERS ###################################################################
################################################################################
# These are the parameters for the user to set

# Path to directory containing spectra. Use quotes and forward slashes.
path <- file.path("data/copyTest")

# Minimum signal-to-noise ratio to define a peak
snr <- 8 

# Minimum number of peaks for an acceptable spectrum
minPeaks <- 10 

# If TRUE, only considers spectra within the bounds below (use FALSE otherwise)
trimSpec <- TRUE 

lowerBound <- 1100 # Lower end of mass range. Choose 0 to include all
upperBound <- Inf  # Upper end of mass range. Choose Inf to include all


################################################################################
################################################################################
################################################################################

# Check Input and Other Requirements -------------------------------------------
cat("\n Checking Parameters...\n")

if(!file.exists(path)) {
  stop("The input directory does not exist.")
}

if(is.na(snr)) {
  stop("The specified SNR was non-numeric.")
}

if(is.na(minPeaks)) {
  stop("The specified minimum number of peaks was non-numeric.")
}

cat(" Loading Required Packages...\n")
suppressPackageStartupMessages(
  if(!require(MALDIquant) | !require(MALDIquantForeign) | !require(plyr)) {
    cat("The MALDIquant, MALDIquantForeign and plyr packages are necessary to use this script. Please install these packages and try again.\n")
  }
)

cat(" Loading Spectra...\n")

# Read in mzXML files ----------------------------------------------------------
file.list <- list.files(path, pattern=".mzXML$",full.names=T)
if(length(file.list) == 0) {
  stop("No files to be analyzed. Please verify that all input files are mzXML")
} else {
  cat(paste0("  -> ", length(file.list)," spectra to analyze\n"))
}
spectra <- importMzXml(file.list)

# Trim Spectra -----------------------------------------------------------------
if(trimSpec) {
  spectra <- trim(spectra, range = c(lowerBound, upperBound))
}

# Normalize MS intensities -----------------------------------------------------
cat(" Processing Spectra... \n")
spectra2 <- transformIntensity(spectra, method="sqrt")

spectra3 <- suppressWarnings(
  smoothIntensity(spectra2, method="SavitzkyGolay", halfWindowSize=5)
)

# Remove Baseline ------------------------------------------------------------
spectra4 <- removeBaseline(spectra3, method="SNIP", iterations=60)

spectra5 <- calibrateIntensity(spectra4,method="TIC")

# Align MS -------------------------------------------------------------------
spectra6 <- alignSpectra(spectra5,reference = spectra5[[1]])


# Pick Peaks -----------------------------------------------------------------
cat(" Picking Peaks... \n")
peaks <- detectPeaks(spectra6, SNR=snr, halfWindowSize=5, method="MAD")
peaks <- binPeaks(peaks, tolerance=0.002)

pVec <- ldply(peaks, function(x){
  num <- length(x@snr)
  file <- x@metaData$file
  data.frame(file = file, num = num)
})

goodSpec <- pVec[pVec$num >= minPeaks, ]
badSpec <- pVec[pVec$num < minPeaks, ]

cat(paste0("  -> ", nrow(goodSpec), " spectra passed\n"))
cat(paste0("  -> ", nrow(badSpec), " spectra failed\n"))
names(goodSpec) <- names(badSpec) <- c("File", "Peaks")

# Write Results ----------------------------------------------------------------
cat(" Writing Results to Tab-Delimited Files...\n")
cat("  -> ", file.path(path,"passingSpectra.txt"), "\n")
cat("  -> ", file.path(path, "failedSpectra.txt"), "\n")

write.table(goodSpec, file.path(path,"passingSpectra.txt"), quote = F, sep = "\t", row.names = F)
write.table(badSpec, file.path(path,"failedSpectra.txt"), quote = F, sep = "\t", row.names = F)


# Write good files to new directory --------------------------------------------
cat(paste0(" Writing ", nrow(goodSpec), " spectra to \"passed\" directory..."))

dir.create(file.path(path,"passed"), showWarnings = F)

naught <- llply(goodSpec$File, function(spec) {
  file.copy(as.character(spec), 
            file.path(path, "passed", gsub("^.*\\\\", "", spec)))
})

cat("\n Done!\n")