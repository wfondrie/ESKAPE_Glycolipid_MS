library(MALDIquant)
library(MALDIquantForeign)
library(pROC)
library(plyr)
library(tidyverse)
library(stringr)
library(forcats)
set.seed(321654)

source("R/extractPeaks.R")
source("R/libraryDotProduct.R")

snr <- 8 # minimum SNR for peak picking

# Import Spectra ---------------------------------------------------------------
path <- c("data/fullLib/Klebsiella pneumoniae - res",
          "data/fullLib/Klebsiella pneumoniae - sen")

files <- list.files(path, pattern = ".mzXML$", full.names = T, recursive = T)

spectra <- importMzXml(files, verbose = F)

# Spectra preprocessing. -------------------------------------------------------
# See MALDIquant documentation for details.

spectra <- transformIntensity(spectra, method="sqrt")
spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
spectra <- removeBaseline(spectra, method="SNIP", iterations=60)
spectra <- calibrateIntensity(spectra,method="TIC") # Intensity proportional to TIC

# The organisms of interest
ooi <- c("Klebsiella pneumoniae - res",
         "Klebsiella pneumoniae - sen")

# Pick and Extract Peaks -------------------------------------------------------
peaks <- detectPeaks(spectra, SNR = snr, halfWindowSize = 10, method = "MAD")
peaks <- binPeaks(peaks, tolerance = 0.5)

# Extract peak data to create single data frames
peakDat <- as_tibble(ldply(peaks, extractPeaks))


# Split Data into Library and Test Sets ----------------------------------------
sampSummary <- peakDat %>% 
  group_by(org, samps, spec) %>%
  summarize(numPeaks = length(mz)) %>%
  mutate(ooi = ifelse(org %in% ooi, T, F)) %>%
  # Pick an equal proportion of random samples from the res and sen spectra
  group_by(org) %>%
  mutate(train = ifelse(ooi,
                        as.logical(rbinom(length(samps), 1, p = 0.5)),
                        T))

# Spectra used to build consensus spectrum
peakLib <- sampSummary %>%
  filter(train) %>%
  select(org, samps, spec, train) %>%
  inner_join(peakDat)

# Consensus spectrum is the sum of all training spectra for each org
# Intensities relative to base peak
peakConsensus <- peakLib %>%
  group_by(org, mz, train) %>%
  summarize(totalInt = sum(relInt)) %>%
  group_by(org) %>%
  mutate(relInt = totalInt / max(totalInt))

# The spectra to test against the consensus spectrum
test <- sampSummary %>%
  filter(!train) %>%
  select(org, samps, spec, train) %>%
  inner_join(peakDat)

# Dot product between test and consensus spectra  ------------------------------
peaksDot <- libraryDotProduct(test, peakConsensus)

# Dot Product Score Scatter ----------------------------------------------------
dot <- peaksDot %>% 
  gather(lib, dotP, -org, -samps, -spec) %>%
  mutate(lib = gsub("_lib", "", lib),
         truth = ifelse(org == lib, T, F))

dot %>%
  select(org, samps, spec, lib, dotP) %>%
  spread(lib, dotP) %>%
  select(org, starts_with("Klebsiella pneumoniae")) %>%
  mutate(`Colistin\nResistance` = ifelse(str_detect(org, "- res"), "+", "-")) %>%
  ggplot(aes(x = `Klebsiella pneumoniae - sen`, 
             y = `Klebsiella pneumoniae - res`, 
             color = `Colistin\nResistance`)) +
  geom_point(size = 0.75) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Klebsiella pneumoniae",
       x = "Sensitive Library Dot Product",
       y = "Resistant Library Dot Product") +
  theme_bw() +
  theme(text = element_text(size = 8),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = c(0.0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.key.size = unit(2, "mm"),
        legend.title = element_text(size = 6))

ggsave("results/scatterKp.pdf",  width = 70, height = 70, unit = "mm", useDingbats = F)

# Summed Spectra After Peak Picking --------------------------------------------
# labels for spectra
resLabs <- data.frame(mz = rep(2400, 2), 
                      relInt = c(1, -1), 
                      txt = c("Resistant", "Sensitive"))


# offsets for labels
xoff <- 25
# sum spectra
summedPeaks <- peakDat %>%
  group_by(org, mz) %>%
  summarize(relInt = sum(relInt)) %>%
  group_by(org) %>%
  mutate(relInt = relInt / max(relInt),
         relInt = ifelse(str_detect(org, "res"), relInt, -relInt),
         lab = ifelse(relInt == max(relInt) &
                        relInt > 0, "1840", ""),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 1970 & mz < 1974])] &
                        relInt > 0, "1971", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 1954 & mz < 1958])] &
                        relInt > 0, "1955", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 2078 & mz < 2082])] &
                        relInt > 0, "2079", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 2062 & mz < 2066])] &
                        relInt > 0, "2063", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 1823 & mz < 1827])] &
                        relInt > 0, "1824", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 1375 & mz < 1379])] &
                        relInt > 0, "1376", lab),
         lab = ifelse(mz == mz[relInt == max(relInt[mz > 1403 & mz < 1407])] &
                        relInt > 0, "1404", lab),
         xoffset = ifelse(lab %in% c("1840", "1404", "1971", "2079"), xoff, 0),
         xoffset = ifelse(lab %in% c("1955", "2063", "1376", "1824"), -xoff, xoffset))

labeledPeaks <- summedPeaks %>%
  ungroup %>%
  filter(lab != "") %>%
  select(`Labelled m/z` = lab, `Measured m/z` = mz, `Relative Intensity` = relInt)

write_csv(labeledPeaks, "results/labeledPeaksKp.csv")

# plot mirrored summed spectrum
summedPeaks %>%
  ggplot(aes(x = mz, ymax = relInt, ymin = 0)) +
  geom_linerange(size = 0.25, aes(color = fct_rev(org))) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black"),
        text = element_text(size = 8),
        legend.position = "none") +
  geom_text(aes(x = mz + xoffset, y = relInt, label = lab), 
            size = 2, nudge_y = 0.05) +
  geom_text(data = resLabs, aes(x = mz, y = relInt, label = txt), 
            hjust = "inward", vjust = "inward") +
  labs(title = "Klebsiella pneumoniae",
       x = "m/z") +
  xlim(c(1000, 2400))

ggsave("results/mirrorSpecKp.pdf", width = 140, height = 70, unit = "mm", useDingbats = F)

