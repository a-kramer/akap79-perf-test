#!/usr/bin/env Rscript
library(rgsl)
library(uqsa)
library(pbdMPI)

start_time <- Sys.time()
comm  <- 0
pbdMPI::init()
r <- pbdMPI::comm.rank(comm=comm)
cs <- pbdMPI::comm.size(comm=comm)
attr(comm,"rank") <- r
attr(comm,"size") <- cs
beta <- 1.0 #(1.0 - (r/cs))^2   # PT: inverse temperature

      N <- 20000           # default sample size

h <- 0.01                  # step size
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"
sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

for (i in seq_along(ex)){
	t_ <- ex[[i]]$outputTimes
	nt <- length(t_)
	BY <- 10 # take every BYth point
	D <- t(ex[[i]]$outputValues)
	D[is.na(D)] <- 0.0
	SD <- D*0.05+apply(D,1,FUN=max,na.rm=TRUE)*0.05 #t(ex[[i]]$errorValues)
	SD[is.na(SD)] <- Inf
	ex[[i]]$time <- t_[seq(1,nt,by=BY)]
	ex[[i]]$data <- D[,seq(1,nt,by=BY)]
	ex[[i]]$stdv <- SD[,seq(1,nt,by=BY)]
}

parMCMC <- log10(sb$Parameter[["!DefaultValue"]])
stopifnot("!Max" %in% names(sb$Parameter))
stopifnot("!Min" %in% names(sb$Parameter))
mu <- 0.5*(log10(sb$Parameter[["!Max"]])+log10(sb$Parameter[["!Min"]]))
stdv <- 0.5*(log10(sb$Parameter[["!Max"]])-log10(sb$Parameter[["!Min"]]))
stopifnot(length(mu)==length(stdv))
stopifnot(length(parMCMC)==length(mu))

dprior <- dNormalPrior(mean=mu,sd=stdv)
rprior <- rNormalPrior(mean=mu,sd=stdv)

gprior <- \(p) {return(-1.0*(p-parMCMC)/stdv^2)}
## ----simulate-----------------------------------------------------------------
sim <- simfi(ex,modelName,log10ParMap) # or simulator.c

## log-likelihood
llf <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(Reduce(\(a,b) a + b$logLikelihood[1], attr(parMCMC,"simulations"), init = 0.0))
}

## gradient of the log-likelihood
gllf <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(as.numeric(Reduce(\(a,b) a + b$gradLogLikelihood[i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J))
}

## Fisher Information
fi <- function(parMCMC){
	i <- seq_along(parMCMC)
	J <- log10ParMapJac(parMCMC)
	return(t(J) %*% Reduce(\(a,b) a + b$FisherInformation[i,i,1], attr(parMCMC,"simulations"), init = 0.0) %*% J)
}

## ----update-------------------------------------------------------------------
smmala <- mcmcUpdate(
	simulate=sim,
	experiments=ex,
	logLikelihood=llf,
	dprior=dprior,
	gradLogLikelihood=gllf,
	gprior,
	fisherInformation=fi,
	fisherInformationPrior=diag(1.0/stdv^2)
)

## ----mcmc-method--------------------------------------------------------------
MC <- mcmc(smmala)
## ----init---------------------------------------------------------------------
	x <- mcmcInit(
		beta,
		parMCMC,
		simulate=sim,
		logLikelihood=llf,
		dprior,
		gllf,
		gprior,
		fi
	)

## beta,parMCMC,simulate,logLikelihood,dprior,gradLogLikelihood=NULL,gprior=NULL,fisherInformation=NULL
A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.5^4 + a^4))
}

## ----converge-and-adapt-------------------------------------------------------
if (r==0){
	cat(sprintf("%10s  %12s %16s %16s\n","rank","iteration","acceptance","step size"))
	cat(sprintf("%10s  %12s %16s %16s\n","----","---------","----------","---------"))
}

for (j in seq(4)){
	for (i in seq(6)){
		s <- MC(x,100,h)           # evaluate rate of acceptance
		a <- attr(s,"acceptanceRate")
		h <- h*A(a)                # adjust h up or down
		x <- attr(s,"lastPoint")   # start next iteration from last point
		cat(sprintf("%10i  %12i %16.4f %16.4g\n",r,i,a,h)) # cat(sprintf("rank: %02i; iteration: %02i; a: %f; h: %g\n",r,i,a,h))
	}
	## ---- here, the sampling happens
	s <- MC(x,N,h) # the main amount of work is done here
	saveRDS(s,file=sprintf("AKAP79-smmala-sample-%i-rank-%i.RDS",j,r))
	x <- attr(s,"lastPoint")
}
time_ <- difftime(Sys.time(),start_time,units="min")
cat("Finished after: ",time_," minutes\n")

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
