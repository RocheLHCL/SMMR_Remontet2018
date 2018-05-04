#---------------------------------------#
#   CI_delta_Standardized_NS.R          |
#---------------------------------------#

# Code to estimate standardized net survival and its confidence interval with the delta method

# It is recommended to read and understand the code in 'CI_delta_NS.R' before reading the code below
#
# If you are NOT interested in standardized net survival, we suggest you to disregard this program file

# The maths to compute standardized NS and the CIs are not sophisticated at all, but tedious (especially for CI)
# this makes this code hardly readable 

#================================================================================

# Standardized net survival at time t for year of diagnosis y is simply a weighted sum of 'individual' net survival SN(t,y,a) over the age a,
# the weights corresponding to the distribution of age of a reference population
# 
# Classicaly, International Cancer Survival Standards (ICSS) weights are used
# Reference: Corazziari et al. 2004. Standard cancer patient population for age standardising survival ratios. Eur J Cancer. 2004 Oct;40(15):2307-16. 


# run 'ana_cervix.R' to load the cervix data and fit the model (object stored as tempmod)
exists("tempmod")
exists("cervix")

#================================================================================

# creation of the data.frame containting the ICSS weights

icss <-data.frame(
		typepop=rep(c("1","2","3","prostate"),each=5),
		ageclasse=c(rep(c("[15;45[","[45;55[","[55;65[","[65;75[","[75;++["),3),
						c("[15;55[","[55;65[","[65;75[","[75;85[","[85;++[")),
		numclas=rep(rep(c(1:5),4)),
		wstd=c(
			c(0.07,0.12,0.23,0.29,0.29),
			c(0.28,0.17,0.21,0.20,0.14),
			c(0.60,0.10,0.10,0.10,0.10),
			c(0.19,0.23,0.29,0.23,0.06))
	)

	
#================================================================================

# table.age: 1 row per age (rounded down) 
la<-floor(cervix$age)
lacl<-cut(la,breaks=c(0,45,55,65,75,200),labels=c("[15;45[","[45;55[","[55;65[","[65;75[","[75;++["),right=FALSE,include.lowest=TRUE)
table.age<-data.frame(age=as.numeric(names(tapply(la,la,sum))),n=tapply(la,la,length),agecl=tapply(lacl,la,FUN=function(x)return(as.character(x[1]))))

# merge of table age and icss, and computation of weights for each age
# There are 3 standards depending on the cancer site; use typepop=2 for cervical cancer
table.age<-merge(table.age,icss[icss$typepop==2,c("ageclasse","wstd")],by.x="agecl",by.y="ageclasse")
ncl<-data.frame(agecl=names(table(lacl)),ncl=as.numeric(table(lacl)))
table.age<-merge(table.age,ncl,by="agecl")
table.age$wstd2<-table.age$wstd * (table.age$n/table.age$ncl)

# ltemps: vector of time 
ltemps<-c(1,5)
# ly: vector of years of diagnosis
ly<-c(1995,2000,2005)
la<-table.age$age

# table.yat: 1 row per time, age, and year of diagnosis
# -> length(ly)*length(ltemps)*length(la) rows
table.yat<-data.frame(temps=rep(ltemps,each=length(la)*length(ly)),adiag=rep(rep(ly,each=length(la)),length(ltemps)),age=rep(la,length(ltemps)*length(ly)))
table.yat<-merge(table.yat,table.age[,c("age","agecl","wstd2")],by="age")
table.yat$agec<-table.yat$age-70
table.yat$adiagc<-table.yat$adiag-2000
table.yat<-table.yat[order(table.yat$temps,table.yat$adiag,table.yat$age),]
#-------------------------------------------------------------------------------------

# MyTP.gam: get the 'design' matrix for centered age 'agec', centered year of diagnosis 'adiagc', at the time=nodes value
MyTP.gam<-function(modele.gam,nodes,agec,adiagc)
{
	return(predict.gam(modele.gam,newdata=data.frame(intnum=nodes,agec=agec,ydiagc=adiagc),type="lpmatrix"))
}

# surv: Net survival computation
# 	mat: 'design' matrix, as returned by function 'MyTP.gam'
#	param: estimated coefficients 
# 	weights: weights of Gauss-Legendre quadrature
surv<-function(mat,param,weights)
{
	return(exp(-sum(exp(mat%*%as.numeric(param))*weights)))
}

# TxCum: cumulative hazard computation, similar to -log(surv(mat,param,weights))
# 	mat: 'design' matrix, as returned by function 'MyTP.gam'
#	param: estimated coefficients 
# 	weights: weights of Gauss-Legendre quadrature
TxCum<-function(mat,param,weights)
{
	return(sum(exp(mat%*%as.numeric(param))*weights))
}

# TauxX: derivative of the cumulative hazard; required for CI computation with Delta method
# 	mat: 'design' matrix, as returned by function 'MyTP.gam'
#	param: estimated coefficients 
# 	weights: weights of Gauss-Legendre quadrature
TauxX<-function(mat,param,weights,ii)
{
	sum(mat[,ii]*(exp(mat%*%as.numeric(param)))*weights)
}

# gauss.quad() returns the nodes and the weights for integration over interval [-1,1]
# 'Rescale' compute the nodes and the weights to integrate over interval [a,b]
# 	gl: object returned by function gauss.quad of package statmod
Rescale<-function(gl,a,b)
{
	gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2
	gl$weights <- gl$weights*(b-a)/2
	return(gl)
}
#-------------------------------------------------------------------------------------

library(mgcv)
library(statmod)

GL<-gauss.quad(n=100,kind="legendre")

table.yat$surv<-NA
table.yat$surv.PH<-NA
ly<-sort(unique(table.yat$adiag))
ltemps<-sort(unique(table.yat$temps))

# sstd: will store Standardized Net Survival estimation and confidence interval
# one row per time and year of diagnosis; --> data.frame with length(ltemps)*length(ly) rows

sstd<-data.frame(
	temps=rep(ltemps,each=length(ly)),adiag=rep(ly,length(ltemps)),
	sstd=rep(-1,length(ltemps)*length(ly)),se.logH=rep(-1,length(ltemps)*length(ly)))

# lcoef: estimated coefficients
lcoef<-tempmod$coefficients
# number of coefficients (6*5*4=120 here)
ncoef<-length(lcoef)

for(temps in ltemps)
{
	# temps<-ltemps[1]
	suppressWarnings(rm(gg,surv0,tempmat,temp2,tagec,tadiagc,PartialSurv,snstd,partialsn,partialLogT))
	
	# we need to rescale the nodes and weights for each time value 'temps'
	gg <- Rescale(GL,0,temps)
	
	for(indy in 1:length(ly))
	{
		# indy<-1
		y<-ly[indy]
		# temptable : will contain net survival (multiplied by standard weights) by time 'temps', year 'y', and age 
		# 1 row per age (for a given time 'temps' and year 'y')
		temptable<-table.yat[table.yat$adiag==y & table.yat$temps==temps,]
		# PartialSurv: derivative of NS(t,a,y) with respect to each parameter (multiplied by standard weights)
		# 1 row per parameter and 1 column per age
		PartialSurv<-matrix(0,ncoef,dim(temptable)[1])
		
		for(x in 1:(dim(temptable)[1]))
		{
			# x<-1
			tagec<-temptable$agec[x]
			tadiagc<-temptable$adiagc[x]
			
			# tempmat: design matrix for centered age 'tagec', year 'tadiagc', and for each time=node value
			tempmat<-MyTP.gam(modele.gam=tempmod,nodes=gg$nodes,agec=tagec,adiagc=tadiagc)
			surv0<-surv(tempmat,lcoef,gg$weights)
			temp2<-rep(0,length(ncoef))
			for(k in 1:ncoef)
			{
				temp2[k]<- - TauxX(tempmat,lcoef,gg$weights,k)
			}
			PartialSurv[,x]<-temp2*surv0*temptable$wstd2[x]
			temptable$surv[x]<-surv0*temptable$wstd2[x]
		}
		
		# snstd: standardized net survival at time 'temps' and year 'y'
		snstd<-sum(temptable$surv)
		
		# derivative of the standardized net survival with respect to each parameter
		partialsn<-as.matrix(apply(PartialSurv,1,sum),ncoef,1)

		# derivative of log-log(Standardized NS), with respect to each parameter
		partialLogT<-partialsn/(snstd*log(snstd))
		
		# variance of log-log(Standardized NS), as obtained by the Delta method
		varLogT<-t(partialLogT)%*%tempmod$Vp%*%partialLogT
		
		sstd$sstd[sstd$temps==temps & sstd$adiag==y]<-snstd
		sstd$se.logH[sstd$temps==temps & sstd$adiag==y]<-sqrt(varLogT)
	}
}

sstd$lower<-exp(-exp( log(-log(sstd$sstd)) - qnorm(0.025) * sstd$se.logH))
sstd$upper<-exp(-exp( log(-log(sstd$sstd)) + qnorm(0.025) * sstd$se.logH))

sstd[,"Std. NS [95% CI]"]<-paste(	format(round(100*sstd$sstd,1),nsmall=1)," [",
									format(round(100*sstd$lower,1),nsmall=1),",",
									format(round(100*sstd$upper,1),nsmall=1),"]",sep="")


sstd$time<-sstd$temps
sstd$year.diag<-sstd$adiag

print(sstd[,c("time","year.diag","Std. NS [95% CI]")])

# results ===>
#
#
#  time year.diag Std. NS [95% CI]
#1    1      1995 86.7 [85.8,87.5]
#2    1      2000 86.6 [85.9,87.3]
#3    1      2005 86.5 [85.6,87.4]
#4    5      1995 66.2 [64.8,67.6]
#5    5      2000 65.0 [63.9,66.0]
#6    5      2005 63.5 [62.1,64.9]
> 




