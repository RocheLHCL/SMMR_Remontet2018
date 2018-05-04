# SMMR_Remontet2018

In this document readme.pdf, we present the R-code that allows reproducing the analysis of the case study (available in the supplementary material) of the following article:
Flexible and structured survival model for a simultaneous estimation of non-linear and non-proportional effects and complex interactions between continuous variables: performance of this multidimensional penalized spline approach in net survival trend analysis.
Laurent Remontet, Zoé Uhry, Nadine Bossard, Jean Iwaz, Aurélien Belot, Coraline Danieli, Hadrien Charvat, Laurent Roche and the CENSUR Working Survival Group, Statistical Methods in Medical Research,2018.


However:
-	due to copyright issues, we cannot provide the original real dataset. So, we provide one of the simulated datasets used in the simulation study on cervix uteri cancer data on 10,000 patients. The results may thus differ, to some extent, from those presented in the article.
-	For simplicity the PH model is not presented


File / directory	Contents
/data/	cervix dataset 
expected hazard of the general population
/function/	Function to split the data
modified link function
ana_cervix.R

= main program	Load and split the data
Model fit
Prediction of excess hazard with CI
Prediction of net survival
Ana_cervix.lis	Results from ana_cervix.R
CI_delta_NS.R	prediction of net survival with CI
CI_delta_Standardized_NS.R	prediction of standardized net survival with CI
FigureS11_Trends_Std_NS.R	plot the figureS11
FigureS12_Trends_NS_by_Age.R	plot the figureS12
FigureS13_ExcessHazard.R	plot the figureS13
/res/	figureS11, figureS12, figureS13 in eps format
CI=confidence intervals

The following versions were used :
R version 3.1.1 (2014-07-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)

attached base packages:
[1] splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] relsurv_2.0-9   date_1.2-34     survival_2.37-7 statmod_1.4.21  mgcv_1.8-3      nlme_3.1-117   

loaded via a namespace (and not attached):
[1] grid_3.1.1      lattice_0.20-29 Matrix_1.2-6    tools_3.1.1    
