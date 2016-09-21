# Analysis of ESKAPE Pathogen Glycolipids

This repository contains the R scripts used in the analysis seen in "Identification of the ESKAPE pathogens by mass spectrometric analysis of microbial membrane glycolipids" by Leung et al.   


### Files
- **analyzeSpectra.R** contains the main function used for preprocessing the glycolipid mass spectra analysed in **LargeHeatmap.R** and **SmallHeatMap.R**. 

- **LargeHeatMap.R** Script that creates the heatmap comparing many glycolipid mass spectra by dot product (Supplementary Figure 5 in paper).  

- **SmallHeatMap.R** Script that creates the ESKAPE pathogen heatmap (Figure 2 in paper) by comparing glycolipid mass spectra by dot product.  

- **filterSpectra.R** Script that preprocesses spectra for peak detection, then returns a list of spectra that pass the given criteria. The user can specify a minimum number of mass peaks at a minimum signal-to-noise ratio in either a defined mass rage or the full spectrum.

### Session Information

Below is the session information from the published work analysis.

```R
> devtools::session_info()
Session info -------------------------------------------------------------------------------------
 setting  value                       
 version  R version 3.3.1 (2016-06-21)
 system   x86_64, mingw32             
 ui       RStudio (0.99.893)          
 language (EN)                        
 collate  English_United States.1252  
 tz       America/New_York            
 date     2016-09-21                  

Packages -----------------------------------------------------------------------------------------
 package            * version  date       source        
 base64enc            0.1-3    2015-07-28 CRAN (R 3.2.1)
 colorspace           1.2-6    2015-03-11 CRAN (R 3.2.1)
 devtools             1.12.0   2016-06-24 CRAN (R 3.3.1)
 digest               0.6.10   2016-08-02 CRAN (R 3.2.5)
 evaluate             0.9      2016-04-29 CRAN (R 3.2.5)
 ggdendro           * 0.1-20   2016-04-27 CRAN (R 3.2.5)
 ggplot2            * 2.1.0    2016-03-01 CRAN (R 3.2.5)
 gridExtra          * 2.2.1    2016-02-29 CRAN (R 3.2.5)
 gtable               0.2.0    2016-02-26 CRAN (R 3.2.5)
 htmltools            0.3.5    2016-03-21 CRAN (R 3.2.5)
 magrittr             1.5      2014-11-22 CRAN (R 3.2.1)
 MALDIquant         * 1.15     2016-06-25 CRAN (R 3.2.5)
 MALDIquantForeign  * 0.10     2015-11-01 CRAN (R 3.2.5)
 MASS                 7.3-45   2016-04-21 CRAN (R 3.3.1)
 memoise              1.0.0    2016-01-29 CRAN (R 3.2.3)
 munsell              0.4.3    2016-02-13 CRAN (R 3.2.5)
 plyr               * 1.8.4    2016-06-08 CRAN (R 3.2.5)
 Rcpp                 0.12.6   2016-07-19 CRAN (R 3.2.5)
 readBrukerFlexData   1.8.2    2014-12-16 CRAN (R 3.2.1)
 readMzXmlData        2.8.1    2015-09-16 CRAN (R 3.2.2)
 reshape2           * 1.4.1    2014-12-06 CRAN (R 3.2.1)
 rJava              * 0.9-8    2016-01-07 CRAN (R 3.2.3)
 rmarkdown            1.0      2016-07-08 CRAN (R 3.2.5)
 scales             * 0.4.0    2016-02-26 CRAN (R 3.2.5)
 stringi              1.1.1    2016-05-27 CRAN (R 3.2.5)
 stringr              1.0.0    2015-04-30 CRAN (R 3.2.1)
 withr                1.0.2    2016-06-20 CRAN (R 3.3.1)
 xlsx               * 0.5.7    2014-08-02 CRAN (R 3.2.1)
 xlsxjars           * 0.6.1    2014-08-22 CRAN (R 3.2.1)
 XML                  3.98-1.4 2016-03-01 CRAN (R 3.3.0)
```