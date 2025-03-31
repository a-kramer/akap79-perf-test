#!/usr/bin/env Rscript
require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 80000 # Size of the sub-sample from each chain
npc <- 500 # pre-calibration sample size
# Define ABC-MCMC Settings
delta <- 0.01

sb <- readRDS("AKAP79-sb.RDS")
ex <- readRDS("AKAP79-ex.RDS")

print(comment(sb))
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"

parVal <- SBtabVFGEN::sbtab_quantity(sb$Parameter)

ll <- log10(sb$Parameter[["!Min"]])
ul <- log10(sb$Parameter[["!Max"]])

set.seed(7619201)

distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
  if (all(is.finite(funcSim))){
    distance <- mean(((funcSim-t(dataExpr))/t(dataErr))^2, na.rm=TRUE)
  } else {
    distance <- Inf
  }
  return(distance)
}

start_time = Sys.time()
options(mc.cores = length(ex))
simulate <- simulator.c(ex,modelName,log10ParMap)
objectiveFunction <- makeObjective(ex, modelName, distanceMeasure, parMap=log10ParMap, simulate)

message("- Initial Prior: uniform product distribution")

rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)
dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)

## Run Pre-Calibration Sampling
message("- Precalibration")
nChains <- max(1,parallel::detectCores() %/% length(ex))
cat("number of chains: ",nChains,"\n")

start_time_preCalibration <- Sys.time()
options(mc.cores=parallel::detectCores())
pC <- preCalibration(objectiveFunction, npc, rprior, rep = 5, num=nChains)
cat("\nPreCalibration:")
print(Sys.time()-start_time_preCalibration)

## Get Starting Parameters from Pre-Calibration
options(mc.cores=min(parallel::detectCores(),length(ex)))
for(j in 1 : nChains){
  stopifnot(all(dprior(pC$startPar[,j])>0))
  cat("Chain", j, "\n")
  cat("\tMin distance of starting parameter for chain",j," = ", min(objectiveFunction(pC$startPar[,j])),"\n")
  cat("\tMean distance of starting parameter for chain",j," = ", mean(objectiveFunction(pC$startPar[,j])),"\n")
  cat("\tMax distance of starting parameter for chain",j," = ", max(objectiveFunction(pC$startPar[,j])),"\n")
}

## Run ABC-MCMC Sampling
cat(sprintf("-Running MCMC chains \n"))
start_time_ABC = Sys.time()
cl <- parallel::makeForkCluster(nChains)
clusterExport(cl, c("objectiveFunction", "pC", "ns", "delta", "dprior"))

out_ABCMCMC <- parLapply(cl, 1:nChains, function(j) ABCMCMC(objectiveFunction, pC$startPar[,j], ns, pC$Sigma, delta, dprior))
stopCluster(cl)

ABCMCMCoutput <- do.call(Map, c(rbind,out_ABCMCMC))
ABCMCMCoutput$scores <- as.numeric(t(ABCMCMCoutput$scores))
end_time = Sys.time()
time_ = end_time - start_time_ABC
print(time_)
cat("\nRegularizations:", ABCMCMCoutput$nRegularizations)
cat("\nAcceptance rate:", ABCMCMCoutput$acceptanceRate)

saveRDS(ABCMCMCoutput,file=sprintf("abc-akap79-%i-chains.RDS",nChains))

end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)
