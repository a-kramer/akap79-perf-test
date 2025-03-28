#!/usr/bin/env Rscript
require(rgsl)
require(SBtabVFGEN)
library(uqsa)
library(parallel)

# Define Number of Samples for the Precalibration (npc) and each ABC-MCMC chain (ns)
ns <- 80000 # Size of the sub-sample from each chain
npc <- 50000 # pre-calibration sample size
# Define ABC-MCMC Settings
delta <- 0.01

sb <- readRDS("AKAP79-sb.RDS")
ex <- readRDS("AKAP79-ex.RDS")

print(comment(sb))
modelName <- "AKAP79"
comment(modelName) <- "./AKAP79.so"

parVal <- SBtabVFGEN::sbtab_quantity(sb$Parameter)

# scale to determine prior values
defRange <- 1000

# Define Lower and Upper Limits for logUniform prior distribution for the parameters
ll <- c(parVal[1:19]/defRange, parVal[20]/1.9, parVal[21]/defRange, parVal[22:24]/1.25, parVal[25:26]/1.5, parVal[27]/2)
ul <- c(parVal[1:19]*defRange, parVal[20]*1.9, parVal[21]*defRange, parVal[22:24]*1.25, parVal[25:26]*1.5, parVal[27]*2)
ll = log10(ll) # log10-scale
ul = log10(ul) # log10-scale

set.seed(7619201)

distanceMeasure <- function(funcSim, dataExpr=Inf, dataErr=Inf){
  if (all(is.finite(funcSim))){
    distance <- mean((1/71.67*(funcSim-as.matrix(dataExpr))/as.matrix(dataErr))^2, na.rm=TRUE)
  } else {
    distance <- Inf
  }
  return(distance)
}

start_time = Sys.time()
options(mc.cores = length(ex))
simulate <- simulator.c(ex,modelName,parMap)
objectiveFunction <- makeObjective(ex, modelName, distanceMeasure, parMap, simulate)

message("- Initial Prior: uniform product distribution")

#rprior <- rUniformPrior(ll, ul)
rprior <- rNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5) 
# with this choice of sd, each component has 98.8% probability
# of being in its interval [ll,ul], 
# and the vector has (98.8%)^27 = 71.2% probability of being in the hyperrectangle
# defined by ll and ul

#dprior <- dUniformPrior(ll, ul)
dprior <- dNormalPrior(mean = (ll+ul)/2, sd = (ul-ll)/5)

## Run Pre-Calibration Sampling
message("- Precalibration")
start_time_preCalibration <- Sys.time()
options(mc.cores=parallel::detectCores())
pC <- preCalibration(objectiveFunction, npc, rprior, rep = 5)
cat("\nPreCalibration:")
print(Sys.time()-start_time_preCalibration)

nChains <- parallel::detectCores() %/% length(ex)
## Get Starting Parameters from Pre-Calibration
options(mc.cores=length(ex))
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

# Save Resulting Samples to MATLAB and R files.
saveRDS(ABCMCMCoutput,file=sprintf("abc-akap79-%i-chains.RDS",nChains))

end_time = Sys.time()
time_ = end_time - start_time
cat("\nTotal time:")
print(time_)
