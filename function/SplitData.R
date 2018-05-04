
# Function to 'split' the data into subject-bands of follow-up
# 
# This function was written for spliting the simulated datasets in the 5 scenarios 
# Likely, user would have to adapt it for another dataset
# 
# bands: vector of times defining the bands of follow-up for the split. Length of each band must be small for reasonable numerical approximation of the cumulative hazard
# relative: relative=TRUE for excess-hazard modelling and net survival estimation; 
# 			if relative=FALSE, mortality hazard expected mortality hazard are set to 0 

split.data <- function(jeudata,exp.haz, bands, relative )
{
# Here, jeudata corresponds to the cervix dataset. 
# exp.haz: data.frame of the 'expected' mortality hazard of the general population,
# In the cervix example, the expected hazard depends on the year of death ANNEE, age of death AGEX, and sex SEXE
# 

	# age: age at diagnosis; fup: observed time elapsed since diagnosis
	# age.exit.r: age at the end of the follow-up; required for the merge of jeudata and exp.haz
jeudata$age.exit <-  jeudata$age + jeudata$fup
jeudata$age.exit <-  ifelse(jeudata$age.exit >= 99, 99, jeudata$age.exit)
jeudata$age.exit.r <- floor(jeudata$age.exit)

	# ydiag: year of diagnosis
	# year.exit.r: calendar year at the end of the follow-up; required for the merge of jeudata and exp.haz 
	# (e.g., if an individual was diagnosed in 1990 and followed-up fup=4.5 years, we need the expected hazard during the year year.exit.r=1994)
jeudata$year.exit <-  jeudata$ydiag + jeudata$fup
jeudata$year.exit.r <- floor(jeudata$year.exit)

# What is needed is the expected mortality hazard at the time of death (for those who died); this is done here
# merge can of course be made by additional variables if the expected hazard depends on other variables (e.g., deprivation index, geographical area, ...)
jeudata <- merge(jeudata,exp.haz,by.x=c("year.exit.r","age.exit.r","sex"),by.y=c("ANNEE","AGEX","SEXE"))

# the data can now be split, using the lexis() function. The include=list() option must include:
# - the expected mortality hazard 'MUA'
# - the variables used in the analysis (here: adiagc and agec)
splitjeudata          <- lexis(entry = 0, exit = fup, fail = life.status, breaks = bands,
                         include = list(Ident, ydiag, ydiagc, sex, age, agec, age.exit, year.exit, life.status, MUA),
                          data = jeudata)

# tik: length of each subject-band of follow-up						  
splitjeudata$tik      <-  splitjeudata$Exit - splitjeudata$Entry


# intnum: corresponds to the time since diagnosis t, used when calling gam()
# We recall that the 'Poisson approach' is just a numerical trick to compute the cumulative hazard H, required for log-likelihood computation
# Setting intnum to the midpoint of each band of follow-up, we actually approximate this integral H by the point-milieu method,
# a classical method for numerical computation of integrals.
# See Remontet et al Statistics in Medicine 2007 for details     						  
splitjeudata$intnum <- (splitjeudata$Entry + splitjeudata$Exit) / 2
splitjeudata[splitjeudata$Fail == 1, c("intnum")] <- splitjeudata[splitjeudata$Fail == 1, c("Exit")]
 
return(splitjeudata)
}
