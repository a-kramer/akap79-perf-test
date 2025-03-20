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
a <- commandArgs(trailingOnly=TRUE)

if (length(a)>0){
      N <- as.integer(a[1])
} else {
      N <- 50000           # default sample size
}

h <- 1e-2                  # step size

beta <- 1.0
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"

sb <- readRDS(file="AKAP79-sb.RDS")
ex <- readRDS(file="AKAP79-ex.RDS")

parMCMC <- log10(sb$Parameter[["!DefaultValue"]]) #c(1.225141, -2.667506, -2.937709, -1.45128, -2.746507, -2.114591, -5.737067, -3.583369, -3.476828, -0.9773003, 0.6060973, -2.554933, -2.148655, 1.192213, -3.232335, 0.8150445, 0.9318317, -0.8488824, 0.2815703, -0.7135722, -1.276805, 1.183659, 2.071277, 0.0576175, -0.4907902) #
#
stopifnot("!Max" %in% names(sb$Parameter))
stopifnot("!Min" %in% names(sb$Parameter))
mu <- 0.5*(log10(sb$Parameter[["!Max"]])+log10(sb$Parameter[["!Min"]]))
stdv <- 0.5*(log10(sb$Parameter[["!Max"]])-log10(sb$Parameter[["!Min"]]))
stopifnot(length(mu)==length(stdv))
stopifnot(length(parMCMC)==length(mu))

dprior <- dNormalPrior(mean=mu,sd=stdv)
rprior <- rNormalPrior(mean=mu,sd=stdv)

## ----simulate-----------------------------------------------------------------
sim <- simcf(ex,modelName,log10ParMap) # or simulator.c

## ----likelihood---------------------------------------------------------------
logLH <- function(y,h,stdv,name){
	n <- sum(!is.na(stdv))
	llf_const <- sum(log(stdv),na.rm=TRUE) + 0.5*log(2*pi)*n
	llf_sq <- 0.5*sum(((y - h)/stdv)^2,na.rm=TRUE)
	return(-llf_const-llf_sq)
}

suppressMessages(
	llf <- logLikelihoodFunc(ex,simpleUserLLF=logLH)
)
## ----update-------------------------------------------------------------------
metropolis <- mcmcUpdate(
	simulate=sim,
	ex=ex,
	logLikelihood=llf,
	dprior=dprior)

## ----mcmc-method--------------------------------------------------------------
mh <- mcmc(metropolis)
## ----init---------------------------------------------------------------------
x <- mcmcInit(
	beta,
	parMCMC,
	simulate=sim,
	logLikelihood=llf,
	dprior)

A <- function(a) { # step-size adjuster
    return(0.5 + a^4/(0.25^4 + a^4))
}

## ----converge-and-adapt-------------------------------------------------------
if (r==0){
	cat(sprintf("%10s  %12s %16s %16s\n","rank","iteration","acceptance","step size"))
	cat(sprintf("%10s  %12s %16s %16s\n","----","---------","----------","---------"))
}
for (j in seq(4)){
	for (i in seq(6)){
		s <- mh(x,300,h)           # evaluate rate of acceptance
		a <- attr(s,"acceptanceRate")
		cat(sprintf("%10i  %12i %16.4f %16.4g\n",r,i,a,h))
		h <- h*A(a)                # adjust h up or down
		x <- attr(s,"lastPoint")   # start next iteration from last point
	}
	s <- mh(x,N,h) # the main amount of work is done here
	saveRDS(s,file=sprintf("akap79-metropolis-hastings-sample-%i-rank-%i.RDS",j,r))
	## ---- when all are done, we load the sampled points from the files but only for the right temperature:
	x <- attr(s,"lastPoint")
}
time_ <- difftime(Sys.time(),start_time,units="min")
print(time_)

cat(sprintf("rank %02i of %02i finished with an acceptance rate of %f and swap rate of %f.\n",round(r),round(cs),attr(s,"acceptanceRate"),attr(s,"swapRate")))
finalize()
