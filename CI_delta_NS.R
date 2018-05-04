#--------------------------#
#   CI_delta_NS.R          |
#--------------------------#



# run 'ana_cervix.R' to load the cervix data and fit the model (object stored as tempmod)
exists("tempmod")

# In 'ana_cervix.R', Net Survival (NS) at t=5 for women aged 60 and diagnosed in 1997 was estimated,
#
# Here, we show how to compute CI of NS estimates using the Delta method

#================================================================================


	# NS estimate at time t=5, for an individual aged 60 and diagnosed in 1997

# (NOTE: this is a copy-paste of the last lines of 'ana_cervix.R')


# gauss.quad() function in package statmod give us the weights and the nodes for GL
library(statmod)
# We used n=100 nodes in the simulation study, but 20 are actually sufficient for accurate numerical integration in our experience
GL<-gauss.quad(n=20,kind="legendre")

Rescale<-function(gl,a,b)
{
	gl$nodes <- gl$nodes*(b-a)/2+(a+b)/2
	gl$weights <- gl$weights*(b-a)/2
	return(gl)
}
gg<-Rescale(gl=GL,a=0,b=5)
# For some function foo(u), integral of foo(u) is computed as the sum of the (gg$weights* foo(gg$nodes)

# Thus, we need to compute the hazard for each value of node:
ltime<-gg$nodes
newdon<-data.frame(intnum=ltime,ydiag=rep(1997,length(ltime)),age=rep(60,length(ltime)))
newdon$ydiagc<-newdon$ydiag-2000
newdon$agec<-newdon$age-70
design.GL<-predict.gam(tempmod,newdata=newdon,type="lpmatrix")
# 6*5*4=120 parameters here, so design.GL should be a matrix with 20 rows (1 row per node) and 120 columns
dim(design.GL)
# hazard.GL: vector of hazard at each node value
hazard.GL<-exp(as.numeric(design.GL%*%tempmod$coefficients))
cum.hazard<-sum(gg$weights*hazard.GL)
net.survival<-exp(-cum.hazard)
#================================================================================


# Delta method is applied on the log of the cumulative hazard scale

# let:
#  - H and (d_H/d_beta_i) be the cumulative hazard and its derivative  with respect to the i th parameter beta_i
#  - (d_logH/d_beta_i) be the derivative of the logarithm of the cumulative hazard H  with respect to the i th parameter beta_i
# The derivative of log(H) is then: (d_logH/d_beta_i) = (d_H/d_beta_i) / H
#
# According the Delta-method, the variance of the log of the cumulative hazard scale is the summation over the couples (i,j) of:
# var.logH = SUM { (d_logH/d_beta_i) * (d_logH/d_beta_j) * cov(beta_i,beta_j)}
# var.logH = 1/(H^2) *  SUM { (d_H/d_beta_i) * (d_H/d_beta_j) * cov(beta_i,beta_j)}
#

# Thus, we need: 
# i) the estimated covariance matrix;
varcov<-tempmod$Vp

# ii) to compute H; this was already done above and result is stored in object 'cum.hazard'

# iii) to compute the derivatives of the cumative hazard H with respect to each parameter

#---------------------------------

	# derivative of the cumative hazard H with respect to beta_i: (d_H/d_beta_i)

# (d_H/d_beta_i) can be seen as the integral over [0,t] of the derivative of the instantaneous hazard (with respect to beta_i)
# and GL quadrature is again used to compute this last integral
# 

# the hazard at the qth node value is: exp(design.GL[q,] %*% tempmod$coefficients)
# Thus, with respect to beta_i, the derivative of the hazard at the qth node value is: design.GL[q,i] * exp(design.GL[q,] %*% tempmod$coefficients)
# This can be compactly written as:
deriv.hazard.GL <- design.GL * exp(as.numeric(design.GL%*%tempmod$coefficients))

# Using GL quadrature, (d_H/d_beta_i) is then the weighted sum: sum(deriv.hazard.GL[,i] * gg$weights)
deriv.cum.hazard <- apply(deriv.hazard.GL,2,FUN=function(x,weight)return(sum(x*weight)),weight=gg$weights)

#---------------------------------

# we can now compute the variance of the log of the cumulative hazard H^
lparam<-tempmod$coefficients
tempmat<-matrix(as.numeric(NA),ncol=length(lparam),nrow=length(lparam))
for(i in 1:length(lparam))
{
for(j in 1:length(lparam))
{
	tempmat[i,j]<-varcov[i,j]*deriv.cum.hazard[i]*deriv.cum.hazard[j]
}
}
var.logH<-sum(tempmat,na.rm=FALSE)/(cum.hazard*cum.hazard)

lower.delta<-exp(-exp(log(cum.hazard) - qnorm(0.025)*sqrt(var.logH) ))
upper.delta<-exp(-exp(log(cum.hazard) + qnorm(0.025)*sqrt(var.logH) ))

cat("\n5-year NS estimate [95% CI]: ",round(100*net.survival,1),"% [",round(100*lower.delta,1),",",round(100*upper.delta,1),"]\n",sep="")

# result ==> 5-year NS estimate [95% CI]: 64.4% [62.4,66.3]


#================================================================================




