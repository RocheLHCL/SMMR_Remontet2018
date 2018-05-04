#-------------------------#
#   ana_cervix.R          |
#-------------------------#


rm(list=ls())
gc()

# this path must be changed 	!!
my.rep<-"S:/etude35/3596_survie_1989_2013/communic/paper_tensor/soumission_SMMR_2018_05_04/github/"




# Functions to split data (lexis() and split.data() functions) and load modified link function, adapted for gam() (package mgcv) 
source(paste(my.rep,"function/Lexis.R",sep=""))
source(paste(my.rep,"function/SplitData.R",sep=""))
source(paste(my.rep,"function/Poisson_gam.R",sep=""))
#-------------------------------------------------------------------------------------------------------------------------------

		# Load the data
		#--------------
	
# Load cervix dataset. This is the first of the 200 simulated datasets of the cervix scenario (N=10000)
load(paste(my.rep,"data/cervix10000_SimulatedDataset.RData",sep=""))
# age: age at diagnosis
# agec: age at diagnosis - 70
# ydiag: year of diagnosis
# ydiagc: year of diagnosis - 2000
# fup : time elapsed since diagnosis (in years) = length of follow-up 
# life.status: 0 for alive at end of follow-up (censored), 1 if dead at time fup

# Load expected hazard of the general population
load(paste(my.rep,"data/ExpectedHazard.RData",sep=""))
# ANNEE: year of death
# AGEX: age at death
# SEXE: 1 for male, 2 for female
# MUA: expected mortality hazard (in deaths per person-year) per year, age, and sex. 
#-------------------------------------------------------------------------------------------------------------------------------


sink(paste(my.rep,"ana_cervix.lis",sep=""))

print(date())
print(sessionInfo())


cervix$status5<-cervix$life.status
cervix$status5[cervix$fup>5]<-0
cervix$time5<-cervix$fup
cervix$time5[cervix$fup>5]<-5

# Deaths (all causes) at 5 years 
cat("\n", paste("N (%) deaths at 5 years: ",sum(cervix$status5)," (",round(100*sum(cervix$status5)/(dim(cervix)[1]),1),")\n",sep=""))

# Persons alive at 5 years 
cat(paste("N (%) alive individuals at 5 years: ",sum(cervix$status5==0)," (",round(100*sum(cervix$status5==0)/(dim(cervix)[1]),1),")\n",sep=""))
#-------------------------

		# split of the data
		#------------------

suppressWarnings(rm(splitdon))
# NOTE: bands should stay rather small so that the numerical integration of the hazard is reasonably accurate
splitdon<-split.data(jeudata=cervix,exp.haz=exp.hazard, bands=c(seq(0,1,by=0.05),seq(1.1,5,by=0.1)), relative=T )
	# dcatt: expected number of deaths due to other causes than cancer (required for the modified link function)
splitdon$dcatt<-splitdon$tik*splitdon$MUA

cat("\n", "---------- Dimension of the splitted dataset -----------","\n")
print(dim(splitdon))

cat("\n", "---------- head of the splitted dataset -----------","\n")
print(head(splitdon))

# example of split of the follow-up for individual Ident=1

cat("\n", "---------- split of the follow-up individual Ident=1 -----------","\n")
print(splitdon[splitdon$Ident==1,])

# Fail: indicator of death WITHIN each subject-band of follow-up; this will be the response variable in the gam() formula
# tik: length of the subject-band of follow-up (used as an offset)
# intnum: the variable to use in the formula of gam() for the time elapsed since diagnosis t  (do NOT use the variable 'Time')

gc()
#-------------------------

		# Knots specification
		#--------------------

# specification of the knots was based on the empirical percentiles observed in the population of patients who died in the article
# 6 knots for time
k.t=quantile(cervix[cervix$status5==1,c("time5")], p=c(0, 0.20, 0.40, 0.60, 0.80, 1))

cat("\n", "---------- Knots of time, age and year of diagnosis -----------","\n")
print(k.t)

# 5 knots for age
k.agec<-quantile(cervix[cervix$status5==1,c("agec")], p=c(0, 0.25, 0.50, 0.75, 1))
print(k.agec)

# 4 knots for the year of diagnosis
k.ydiagc<-quantile(cervix[cervix$status5==1,c("ydiagc")], p=c(0, 0.33, 0.66, 1))
print(k.ydiagc)


#-------------------------------------------------------------------------------------------------------------------------------

		# fit of the model
		#-----------------

library(mgcv)
suppressWarnings(rm(tempmod))

	# log of the excess hazard is modelled here as a tensor product of time, age, and year of diagnosis.

	# Warning: results are not guaranteed at all if scale is not set to 1
time.start<-Sys.time()

# took about 21 min for the fitting with N=10000 patients but 467,986 lines (memory required 4 GB)
# Note: with esophagus, N=10000 patients and nearly 200,000 lines, run time is # 3 minutes

# Run with 8 GB of RAM and processor: Inter(R) Core(TM) i7-4770 CPU @ 3.4 GHz 3.4 GHz)
tempmod<-try(gam(Fail~te(intnum,agec,ydiagc, k=c(6,5,4),bs='cr'),offset=log(tik),scale=1,method="REML",
		family=POISS.RS.SPLIT.R.GAM(MyData.dcatt=splitdon$dcatt),data=splitdon,
		knots=list(intnum=k.t, agec=k.agec, ydiagc=k.ydiagc )))

cat("\n", "---------- Summary of the gam fit and run time -----------","\n")
print(summary(tempmod))
time.end<-Sys.time()
difftime(time.end,time.start,units="min")

#================================================================================

		# Estimate excess hazard
		#-----------------------


	# say we want to estimate excess hazard at time t=1, for an individual aged 60 and diagnosed in 1997

# the design matrix can be obtained using predict.gam() with the type="lpmatrix" option
newdon<-data.frame(intnum=1,ydiag=1997,ydiagc=1997-2000,age=60,agec=60-70)
myDesign<-predict.gam(tempmod,newdata=newdon,type="lpmatrix")
# 6*5*4=120 parameters here, so myDesign should be a matrix with 1 row and 120 columns

cat("\n", "---------- dim(myDesign) -----------","\n")
dim(myDesign)

# estimated parameters
lparam<-tempmod$coefficients
log.excess.hazard<-as.numeric(myDesign%*%lparam)
excess.hazard<-exp(log.excess.hazard) 

# CI is easy to compute here since log.excess.hazard is a linear combination of the parameters
# Estimated covariance matrix: 
varcov<-tempmod$Vp
# Standard errors on the log excess hazard scale
se.logh<-as.numeric(sqrt(myDesign%*%varcov%*%t(myDesign)))

lower<- exp(log.excess.hazard-1.96*se.logh)
upper<- exp(log.excess.hazard+1.96*se.logh)

cat("\n", "---------- excess hazard (t=1,age=60,year=1997) -----------","\n")
cat(paste("Estimate [95% CI]: ",round(excess.hazard,3)," [",round(lower,3),",",round(upper,3),"]\n\n",sep=""))

#---------------------

	# plot of the dynamics of excess hazard  for an individual aged 60 and diagnosed in 1997

ltime<-seq(0,5,0.1)
newdon<-data.frame(intnum=ltime,ydiag=rep(1997,length(ltime)),age=rep(60,length(ltime)))
newdon$ydiagc<-newdon$ydiag-2000
newdon$agec<-newdon$age-70
myDesign<-predict.gam(tempmod,newdata=newdon,type="lpmatrix")
# 6*5*4=120 parameters here, so myDesign should be a matrix with length(ltime) rows and 120 columns
cat("\n", "---------- dim(myDesign) -----------","\n")
dim(myDesign)

excess.hazard<-as.numeric(exp(myDesign%*%tempmod$coefficients))
par(mgp=c(2,0.5,0),mar=c(3.5,3.5,2.5,0))
plot(ltime,excess.hazard,xlim=c(0,5),ylim=c(0,0.3),lwd=1.2,cex.axis=0.85,cex.lab=1,cex.main=1,type="l",
	main="Dynamics of excess hazard for women aged 60 and diagnosed in 1997",xlab="Time since diagnosis",ylab="Excess hazard")
#================================================================================

		# Estimate net survival (NS)
		#---------------------------

	# This time, we want to estimate NS at time t=5, for an individual aged 60 and diagnosed in 1997

# Additional work required because we have to compute the integral of the excess hazard to compute NS estimates
# We use here Gauss-Legendre (GL) quadrature 

# gauss.quad() function in package statmod give us the weights and the nodes for GL
library(statmod)
# We used n=100 nodes in the simulation study, but 20 are actually sufficient for accurate numerical integration in our experience
GL<-gauss.quad(n=100,kind="legendre")

# gauss.quad() returns the nodes and the weights for integration over interval [-1,1]
# we need to 'Rescale' the nodes and the weights to integrate over interval [a,b]; this is done here
Rescale<-function(gl,a,b)
{
	gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2
	gl$weights <- gl$weights*(b-a)/2
	return(gl)
}
gg<-Rescale(gl=GL,a=0,b=5)

# For some function foo(u), integral of foo(u) is computed as the sum of the (gg$weights* foo(gg$nodes)
#
## Simple example to understand GL with foo(u)=u*u:
## foo<-function(u)return(u*u)
## cat(paste("\n\nIntegral of u^2 over [0,1] computed with GL: ",sum(gg$weights*foo(gg$nodes)),"\n",sep=""))
## cat(paste("True value: ",(5^3)/3 - (0^3)/3,"\n\n",sep=""))


# Thus, we need to compute the hazard for each value of node:
ltime<-gg$nodes
newdon<-data.frame(intnum=ltime,ydiag=rep(1997,length(ltime)),age=rep(60,length(ltime)))
newdon$ydiagc<-newdon$ydiag-2000
newdon$agec<-newdon$age-70
design.GL<-predict.gam(tempmod,newdata=newdon,type="lpmatrix")
# 6*5*4=120 parameters here, so design.GL should be a matrix with 100 rows (1 row per node) and 120 columns
cat("\n", "---------- design.GL -----------","\n")
dim(design.GL)
# hazard.GL: vector of hazard at each node value
hazard.GL<-exp(as.numeric(design.GL%*%tempmod$coefficients))

cum.hazard<-sum(gg$weights*hazard.GL)
net.survival<-exp(-cum.hazard)
cat("\n","5-year Net Survival estimate for women aged 60 y diagnosed in 1997: ",round(100*net.survival,1),"%\n",sep="")
#-------------------------

	# Confidence interval of NS estimates

# 2 methods can be used to compute Confidence Interval of NS estimates: Delta method and Monte-Carlo simulation
# Delta method was used in the simulation study because it is much less computationally intensive 
# 
# See file 'CI_delta_NS.R' for the code to obtain CI using the delta method

		
#================================================================================

sink()


