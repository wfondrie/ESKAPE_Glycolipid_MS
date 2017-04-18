# Analysis of ESKAPE Pathogen Glycolipids

This repository contains the R scripts used in the analysis seen in "Identification of the ESKAPE pathogens by mass spectrometric analysis of microbial membrane glycolipids" by Leung et al.   


### Files

- **analyzeSpectra.R** - Script that contains the main function used for preprocessing the glycolipid mass spectra analysed in **LargeHeatmap.R** and **SmallHeatMap.R**. 

- **LargeHeatMap.R** - Script that creates the heatmap comparing many glycolipid mass spectra by dot product (Supplementary Figure 5 in paper).  

- **SmallHeatMap.R** - Script that creates the ESKAPE pathogen heatmap (Figure 2 in paper) by comparing glycolipid mass spectra by dot product.  

- **filterSpectra.R** - Script that preprocesses spectra for peak detection, then returns a list of spectra that pass the given criteria. The user can specify a minimum number of mass peaks at a minimum signal-to-noise ratio in either a defined mass rage or the full spectrum.

- **ResistanceAnalysisKp.R** and **ResistanceAnalysisAb** - Scripts that create the figures comparing colistin-sensitive and -susceptible Klebsiella *pneumoniae* and Acinetobacter *baumannii*, respectively. Together they create everything for figures 4 and S5 in the manuscript.  