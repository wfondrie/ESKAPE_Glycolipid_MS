# Extracts m/z, peak intensity sample information from a MassPeaks object.

# Input:
# specPeaks - a MassPeaks object

# Output:
# A long data frame containing 3 columns of sample information, an m/z column and
# a relative intensity column. Intensities are relative to the base peak.

extractPeaks <- function(specPeaks) {
  # extract m/z
  mz <- as.numeric(specPeaks@mass)
  
  # extract relative intensity
  relInt <- specPeaks@intensity
  relInt <- relInt/max(relInt) # Relative intensity to base peak
  
  # extract sample data
  fileName <- metaData(specPeaks)$file
  samps <- str_match(fileName, ".*(\\\\|/)(.*).mzXML$")[ , 3] # Extracts sample name
  org <- str_match(fileName, 
                   ".*(\\\\|/)(.*)(\\\\|/).*mzXML$")[ , 3] # Extracts organism
  spec <- metaData(specPeaks)$nu
  
  return(data.frame(org = org, samps = samps, spec = spec, mz = mz, relInt = relInt))
}