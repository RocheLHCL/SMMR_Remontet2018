
# WARNING : do NOT change the link="log" option
# MyData.dcatt : vector of expected number of deaths due to other causes than cancer. See SplitData.R for details
# 
# If an error is cast during initialisation, this is likely du to a misspecification of the 'initialize' expression (lines 63-73)
# If so, try to set other constant values in:  mustart <- y + 0.3 (e.g., mustart <- y + 0.5)

POISS.RS.SPLIT.R.GAM <- function(link="log", MyData.dcatt)
{
			# The first lines correspond to usual specifications of a link function for a GLM model
	poissRS <- make.link("log")
	poissRS$linkfun <- function(mu, dcatt=MyData.dcatt)
	{
		log(mu-dcatt)
	}
	poissRS$linkinv <- function(eta, dcatt=MyData.dcatt) exp(eta)+dcatt
	poissRS$mu.eta <- function(eta) exp(eta)


			# The following lines are specific to GAM using mgcv package
			# The algorithm developped by Wood for estimation of GAM is based on 
			# Outer Newton-Raphson iterations combined with Inner P-IRLS iterations
		# This requires computation of 2nd (d2link), 3rd (d3link) and 4th (d4link) derivatives of the link function
		# see ? fix.family.link in mgcv package
	poissRS$d2link <- function(mu, dcatt=MyData.dcatt) -1/(mu-dcatt)^2
	poissRS$d3link <- function(mu, dcatt=MyData.dcatt) 2/(mu-dcatt)^3
	poissRS$d4link <- function(mu, dcatt=MyData.dcatt) -6/(mu-dcatt)^4
	
		# Also requires computation of 1st (dvar), 2nd (d2var), and 3rd (d3var) derivatives of variance
		# see ? fix.family.var in mgcv package
	poissRS$dvar <- function(mu) rep.int(1, length(mu))
	poissRS$d3var <- poissRS$d2var <- function(mu) rep.int(0, length(mu))

		# ls : log-saturated likelihood ; required for REML criterion
		# see ?fix.family.ls in mgcv package
	poissRS$ls <- function(y, w, n, scale)
	{
		res <- rep(0, 3)
		res[1] <- sum(dpois(y, y, log = TRUE) * w)
		res
	}

		# see ?fix.family.qf and ?fix.family.rd for qf and rd
		# Objects qf and rd are not actually required in the estimation process 
	poissRS$qf <- function(p, mu, wt, scale)
	{
		qpois(p, mu)
	}
	poissRS$qf <- function(p, mu, wt, scale)
	{
		qpois(p, mu)
	}
	poissRS$rd <- function(mu, wt, scale)
	{
		rpois(length(mu), mu)
	}
 
	variance <- function(mu) mu
	validmu <- function(mu) all(mu > 0)
	dev.resids <- function(y, mu, wt) 2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
	aic <- function(y, n, mu, wt, dev) -2 * sum(dpois(y, mu, log = TRUE) * wt)
	
	initialize <- expression(
		{
			if (any(y < 0)) stop("negative values not allowed for the Poisson family")
			n <- rep.int(1, nobs)
				# !!
				# we must have all(mustart - MyData.dcatt > 0)=TRUE
				# Thus, if y=0, 0.3 must be > MyData.dcatt
				# change this 0.3 to some convenient value if this is not the case
			mustart <- y + 0.3
			# mustart <- y + 0.5
		})
	# simfun <- function(object, nsim)
	# {
		# wts <- object$prior.weights
		# if (any(wts != 1)) 
			# warning("ignoring prior weights")
		# ftd <- fitted(object)
		# rpois(nsim * length(ftd), ftd)
	# }

	structure(list(family="poissonRS", link="poissonRS.link", linkfun=poissRS$linkfun, linkinv=poissRS$linkinv,
		d2link=poissRS$d2link, d3link=poissRS$d3link, d4link=poissRS$d4link, ls=poissRS$ls, 
		variance=variance, dvar=poissRS$dvar, d2var=poissRS$d2var, d3var=poissRS$d3var,
		dev.resids=dev.resids, aic=aic, qf=poissRS$qf, rd=poissRS$rd,
		mu.eta=poissRS$mu.eta, initialize=initialize, validmu=validmu,
		valideta=poissRS$valideta),class="family")
}

