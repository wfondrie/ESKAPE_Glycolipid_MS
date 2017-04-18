# Analysis of ESKAPE Pathogen Glycolipids

This repository contains the R scripts used in the analysis seen in "Identification of the ESKAPE pathogens by mass spectrometric analysis of microbial membrane glycolipids" by Leung et al.   


### Files
#### Helper functions
- **analyzeSpectra.R** - Contains the main function used for preprocessing the glycolipid mass spectra analysed in **LargeHeatmap.R** and **SmallHeatMap.R**. 

- **extractPeaks.R** - Contains the function used to extract *m/z*, relative intensity, and sample data from a MassPeaks object in the MALDIQuant package into a data frame. This function is used in **ResistanceAnalysisKp.R** and **ResistanceAnalysisAb.R**.  

- **libraryDotProduct.R** - Takes a long-formatted data frame containing one or multiple test spectra and a long-formatted data frame containing one or multiple library spectra and computes the a dot product between each test spectra and library spectra. The result is a wide-formatted data frame with rows containing each test sample and columns that contain the calculated dot product for each library spectrum. This function is used in **ResistanceAnalysisKp.R** and **ResistanceAnalysisAb.R**.  
  


#### Analysis Scripts
- **LargeHeatMap.R** - Script that creates the heatmap comparing many glycolipid mass spectra by dot product (Supplementary Figure 5 in paper).  

- **SmallHeatMap.R** - Script that creates the ESKAPE pathogen heatmap (Figure 2 in paper) by comparing glycolipid mass spectra by dot product.  

- **filterSpectra.R** - Script that preprocesses spectra for peak detection, then returns a list of spectra that pass the given criteria. The user can specify a minimum number of mass peaks at a minimum signal-to-noise ratio in either a defined mass rage or the full spectrum.

- **ResistanceAnalysisKp.R** and **ResistanceAnalysisAb.R** - Scripts that create the figures comparing colistin-sensitive and -susceptible Klebsiella *pneumoniae* and Acinetobacter *baumannii*, respectively. Together they create everything for figures 4 and S5 in the manuscript.  