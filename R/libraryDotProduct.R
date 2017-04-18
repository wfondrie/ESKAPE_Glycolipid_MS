# Calculates the dot product between data frames containing test spectra
# and library consensus spectra

# Input:
# testDf - Data frame containing the test spectra. Must contain columns "org",
# "samps", "spec", and "test" for unique identifiers and columns "mz", and 
# "relInt" storing spectral information.
#
# consensusDf - Data frame with the same format as testDf

# Output:
# A data frame containing the calculated dot product where rows are 
# test spectra and columns are library spectra

libraryDotProduct <- function(testDf, consensusDf) {
  widePeaks <- testDf   %>%
    full_join(consensusDf) %>% # Combined into one data frame
    mutate(id = ifelse(train,  # Create unique id's for row names
                       paste0(org, "_lib"),
                       paste(org, samps, spec, "test", sep = "_"))) %>%
    select(id, org, samps, spec, mz, relInt) %>% # drop any irrelevant columns
    spread(key = mz, value = relInt, fill = 0) %>% # long to wide format
    ungroup()
  
  wideLib <- widePeaks %>% filter(is.na(samps)) # consensus spectra
  wideTest <- widePeaks %>% filter(!is.na(samps)) # test spectra
  
  # Create matrices to calculate dot product
  # - rows are spectra, columns are m/z, values are intensity
  # - unique id in row names
  libMat <- as.matrix(select(wideLib, -id, -org, -samps, -spec)) 
  row.names(libMat) <- wideLib$id
  
  testMat <- as.matrix(select(wideTest, -id, -org, -samps, -spec))
  row.names(testMat) <- wideTest$id
  
  # Normalize each spectrum to a unit vector
  libMat <- aaply(libMat, 1, function(x) x / sqrt(x %*% x))
  testMat <- aaply(testMat, 1, function(x) x / sqrt(x %*% x))
  
  # Calculate the dot product
  # - results in rows as test spectra, columns as consensus spectra
  # and values being the dot product between them
  dot <- as.data.frame(testMat %*% t(libMat))
  dot$id <- row.names(dot)
  
  # Turn matrix back into data frame, with all info
  r <- wideTest %>%
    select(id, org, samps, spec) %>%
    inner_join(dot) %>%
    select(-id)
  return(r)
}