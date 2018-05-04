# SMMR_Remontet2018

This repository contains the R-code to perform the analysis of the case study (available in the supplementary material) of the following article:

Flexible and structured survival model for a simultaneous estimation of non-linear and non-proportional effects and complex interactions between continuous variables: performance of this multidimensional penalized spline approach in net survival trend analysis.

Authors: Laurent Remontet, Zoé Uhry, Nadine Bossard, Jean Iwaz, Aurélien Belot, Coraline Danieli, Hadrien Charvat, Laurent Roche and the CENSUR Working Survival Group
Statistical Methods in Medical Research,2018 (in revision).


However:
-	due to copyright issues, we cannot provide the original real dataset. Thus, we uploaded one of the simulated datasets used in the simulation study on cervix uteri cancer data on 10,000 patients. The results may thus differ, to some extent, from those presented in the article.
-	For simplicity the PH model is not presented

This repository is organized as follows:
- /data/: cervix dataset and expected hazard of the general population
- /function/: Two functions to split the data and the modified link function
- ana_cervix.R: main program (load and split the data, fit the model, prediction of excess hazard with CI, and prediction of net survival
- Ana_cervix.lis: Results from ana_cervix.R
- CI_delta_NS.R: prediction of net survival with CI (using Delta method)
- CI_delta_Standardized_NS.R	prediction of standardized net survival with CI
- FigureS11_Trends_Std_NS.R to plot the figureS11
- FigureS12_Trends_NS_by_Age.R to plot the figureS12
- FigureS13_ExcessHazard.R to plot the figureS13
- /res/: figureS11, figureS12, figureS13 in eps format


The following versions were used :
R version 3.1.1 (2014-07-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] relsurv_2.0-9   date_1.2-34     survival_2.37-7 statmod_1.4.21  mgcv_1.8-3      nlme_3.1-117   

loaded via a namespace (and not attached):
[1] grid_3.1.1      lattice_0.20-29 Matrix_1.2-6    tools_3.1.1    
