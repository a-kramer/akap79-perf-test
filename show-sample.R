library(uqsa)
library(rgsl)
library(hadron)

if ("f" %in% ls()){
	message("using f variable")
}else {
	f <- "AKAP79-temperature-ordered-smmala-sample-2-for-rank-0.RDS"
}

if (length(f)>1){
	s <- gatherSample(f,1.0)
} else {
	s <- readRDS(file=f)
}

l <- attr(s,"logLikelihood")

plot(l,type='l',xlab='Markov chain index')

res <- hadron::uwerr(data=l,pl=TRUE)
tau <- ceiling(res$tauint + res$dtauint)

cat(sprintf("auto-correltation: %i\n",tau))
hexbin::hexplom(s[seq(1,NROW(s),by=tau),seq(6)])

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
	ex[[i]]$data <- D[,seq(1,nt,by=BY),drop=FALSE]
	ex[[i]]$stdv <- SD[,seq(1,nt,by=BY),drop=FALSE]
}

sim <- simcf(ex,modelName,log10ParMap) # or simulator.c

i <- order(l,decreasing=TRUE)
p <- as.numeric(s[i[1],])

y <- sim(t(s[seq(1,NROW(s),by=tau),]))

plotTimeSeriesBase <- function(simulations, experiments, nmax=NULL, ylimit=NULL){
	nf <- dim(simulations[[1]]$func)[1]
	par(mfrow=c(nf,length(experiments)))
	for (i in seq(length(experiments))){
		time <- experiments[[i]]$time #%otherwise% experiments[[i]]$outputTimes
		n <- dim(simulations[[i]]$func)
		if (is.null(nmax)) nmax=n[3]
		Names <- names(experiments[[i]]$outputValues)
		for (j in seq(n[1])){
			o <- experiments[[i]]$data[j,] #%otherwise% experiments[[i]]$outputValues[[j]]
			e <- experiments[[i]]$stdv[j,] #%otherwise% experiments[[i]]$errorValues[[j]]
			if (is.null(ylimit) || length(ylimit)<j || any(is.na(ylimit[[j]])) || !all(is.finite(ylimit[[j]]))){
				yl <- c(min(o-e),max(o+e))
			} else {
				yl <- ylimit[[j]]
			}
			plot(time,o,type='p',ylim=yl,xlab="t",ylab=Names[j],main=names(experiments)[i])
			arrows(time,o,time,o+e,angle=90,length=0.01)
			arrows(time,o,time,o-e,angle=90,length=0.01)
			for (k in seq(nmax)){
				y<-as.numeric(simulations[[i]]$func[j,,k])
				lines(time,y,col=rgb(0.3,0.6,0.9,0.4));
			}
		}
	}
}

png(file=sub("RDS$","png",f),width=length(ex)*1024,height=NROW(y[[1]]$func)*786,res=150)
plotTimeSeriesBase(y,ex,ylimit=rep(list(c(0,200)),length(ex)))
dev.off()

