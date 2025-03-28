#!/usr/bin/env Rscript
library(uqsa)
library(SBtabVFGEN)
library(ggplot2)

CA <- commandArgs(trailingOnly=TRUE)
if (length(CA)==2){
	files <- commandArgs(trailingOnly=TRUE)
} else {
	files <- c("AKAP79-temperature-ordered-pt-smmala-sample-4-for-rank-0.RDS", "AKAP79-temperature-ordered-pt-smmala-sample-1-for-rank-15.RDS")
}
PREFIX <- uqsa::determinePrefix(files)
print(files)
posterior <- readRDS(files[1])
cat("dim(posterior): ",dim(posterior),"\n")
i <- seq(1,NROW(posterior),by=10)
posterior<-posterior[i,]

prior <- readRDS(files[2])
cat("dim(prior): ",dim(prior),"\n")
i <- seq(1,NROW(prior),by=10)
prior <- prior[i,]

ex <- readRDS("AKAP79-ex.RDS")
sb <- readRDS("AKAP79-sb.RDS")

Names <- rownames(sb$Parameter)
colnames(posterior) <- Names
colnames(prior) <- Names

## ---- R functions for this model:
modelName <- "AKAP79"
C <- NCOL(posterior)
n <- NROW(posterior)
m <- NCOL(posterior)

up <- function(x,y,num=10,subscripts,...){
	m <- mean(x)
	LIM <- c(range(x),range(y))
	dx <- diff(range(x))/2
	if (min(y) < min(x)){
		poly_X <- LIM[c(1,1,2,2)]
		poly_Y <- LIM[c(1,3,3,2)]
	} else {
		poly_X <- LIM[c(1,2,2)]
		poly_Y <- LIM[c(1,1,2)]
	}
	k1 <- MASS::kde2d(head(x,num),head(y,num),n=500,lims=LIM)
	k2 <- MASS::kde2d(tail(x,num),tail(y,num),n=500,lims=LIM)
	C <- contourLines(k2)
	colPatch <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	colContour <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	image(k1$x,k1$y,k1$z,add=TRUE,col=colPatch)
	for (i in seq_along(C)){
		lines(C[[i]],col=colContour[i])
	}
	polygon(poly_X,poly_Y,col=rgb(0.2,1,0.2,0.05))
}

lp <- function(x,y,num=10,...){
	rx <- range(x)
	ry <- range(y)
	C <- cor(head(x,num),head(y,num)) # in [-1,1]
	i <- 1+round((1+C)*50)            # 1+C in [0,2], (1+C)*50 in [0,100], 1+round((1+C)*50) in [1,101]
	col <- colorspace::diverging_hcl(101,rev=FALSE,palette="Blue Red") # maps to possible correlation values
	polygon(c(rx,rev(rx)),rep(ry,each=2),col=col[i])
	text(mean(rx),mean(ry),sprintf("%.3f",C),cex=2+round(abs(C)*2))
}

png(file=sprintf("%s-posterior-all-pairs-and-correlations.png",PREFIX),width=m*360,height=m*360,res=150)
print(pairs(rbind(posterior,prior),upper.panel=up,lower.panel=lp,num=NROW(posterior),main=modelName))
dev.off()

png(file=sprintf("%s-posterior-4-pairs-and-corr.png",PREFIX),width=4*480,height=4*480,res=200)
i <- c(3,6,9,12)
print(pairs(rbind(posterior[,i],prior[,i]),upper.panel=up,lower.panel=lp,num=NROW(posterior),main=modelName))
dev.off()


singlePanel <- function(x,y,num=10,...){
	m <- mean(x)
	LIM <- c(range(x),range(y))
	dx <- diff(range(x))/2
	if (min(y) < min(x)){
		poly_X <- LIM[c(1,1,2,2)]
		poly_Y <- LIM[c(1,3,3,2)]
	} else {
		poly_X <- LIM[c(1,2,2)]
		poly_Y <- LIM[c(1,1,2)]
	}
	k1 <- MASS::kde2d(head(x,num),head(y,num),n=1500,lims=LIM)
	k2 <- MASS::kde2d(tail(x,num),tail(y,num),n=500,lims=LIM)
	C <- contourLines(k2)
	colPatch <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	colContour <- c("#FFFFFF",colorspace::sequential_hcl(length(C),rev=TRUE,palette="Blues 3"))
	image(k1$x,k1$y,k1$z,add=FALSE,col=colPatch,xlim=LIM[1:2],ylim=LIM[3:4],...)
	for (i in seq_along(C)){
		lines(C[[i]],col=colContour[i])
	}
	polygon(poly_X,poly_Y,col=rgb(0.2,1,0.2,0.03),border=NA)
}

hist2 <- function(M,...){
	N <- 10
	L <- c(min(as.numeric(M)),max(as.numeric(M)))
	br <- seq(L[1],L[2],length.out=N+1)
	G <- grey.colors(NCOL(M))
	A <- matrix(NA,NCOL(M),N)
	H <- list()
	for (i in seq(NCOL(M))){
		H[[i]] <- hist(M[i],breaks=br,plot=FALSE)
##		A[]
	}
	barplot(rbind(H[[1]]$counts,H[[2]]$counts),beside=TRUE,names.arg=diff(H[[1]]$breaks)*0.5+head(H[[1]]$breaks,-1),main=colnames(posterior)[8],ylab="count",xlab="log10(par[8])")
legend("topright",legend=c("prior","posterior"),fill=c(G[1],G[2]))
}

png(file=sprintf("%s-posterior-prior-8-9.png",PREFIX),width=3000,height=1000,res=300)
par(mfrow=c(1,3),bty="n")
xl <- paste0("log10(",colnames(posterior)[8],")")
yl <- paste0("log10(",colnames(posterior)[9],")")
#XLim <- (c(min(c(posterior[,8],prior[,8]),max(c(posterior[,8],prior[,8]))))
#YLim <- (c(min(c(posterior[,9],prior[,9])),max(c(posterior[,9],prior[,9]))))
singlePanel(c(posterior[,8],prior[,8]),c(posterior[,9],prior[,9]),NROW(posterior),xlab=xl,ylab=yl,main="kernel density estimate") #,xlim=XLim,ylim=YLim)
H <- list()
H[[1]] <- hist(prior[,8],plot=FALSE)
H[[2]] <- hist(posterior[,8],breaks=H[[1]]$breaks,plot=FALSE)
H[[3]] <- hist(prior[,9],plot=FALSE)
H[[4]] <- hist(posterior[,9],breaks=H[[3]]$breaks,plot=FALSE)
G <- grey.colors(2)
barplot(rbind(H[[1]]$counts,H[[2]]$counts),beside=TRUE,names.arg=diff(H[[1]]$breaks)*0.5+head(H[[1]]$breaks,-1),main=colnames(posterior)[8],ylab="count",xlab="log10(par[8])")
legend("topright",legend=c("prior","posterior"),fill=c(G[1],G[2]))
barplot(rbind(H[[3]]$counts,H[[4]]$counts),beside=TRUE,names.arg=diff(H[[3]]$breaks)*0.5+head(H[[4]]$breaks,-1),main=colnames(posterior)[9],ylab="count",xlab="log10(par[9])")
legend("topright",legend=c("prior","posterior"),fill=c(G[1],G[2]))
dev.off()
