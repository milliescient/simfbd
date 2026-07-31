#!/usr/bin/Rscript

suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(ggplot2)))

args = commandArgs(trailingOnly=TRUE)
dir = args[1]

burnin=2000

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

baseLeft <- function(x, n){
  substr(x, 1, nchar(x)-1)
}

mse = function(f) {
	lm = read.table(paste(dirname(f),"/lambdamu.dat",sep=""),head=T)
	truth = c(lm$lamda,lm$mu)
        df = read.table(f,head=T)[-burnin,]
	df = df[,c(grep("lambda",names(df)),grep("mu",names(df)))]
	err = apply(df, FUN=function(x){ (x - truth)^2 },1)
	mean( apply(err,FUN=mean,1) )
}

medians = function(f) {
	df = read.table(f,head=T)[-burnin,]
        df = df[,c(grep("lambda",names(df)),grep("mu",names(df)))]
	apply(df,FUN=median,2)
}

pyrate.files <- list.files(dir, pattern="pyrate.log", full.names=TRUE, recursive=TRUE)
fbd.files <- paste(dirname(pyrate.files),"/fbd.log",sep="")
bds.files <- paste(dirname(pyrate.files),"/bds.log",sep="")
psi.files <- paste(dirname(pyrate.files),"/psi.dat",sep="")

### Get psi

psi <- as.vector(sapply(psi.files, function(x){ read.table(x,header=F)$V1[1] }))

### Medians

#pdf(paste(dir,"_medians.pdf",sep=""))

#pyrate.medians <- medians(paste(dir,"/pyrate.log",sep=""))
#bds.medians <- medians(paste(dir,"/bds.log",sep=""))

#plot(pyrate.medians,bds.medians)
#dev.off()

### Error plot

pdf("mse.pdf")

fbd.error <- as.vector(sapply(fbd.files, FUN=mse) )
bds.error <- as.vector(sapply(bds.files, FUN=mse) )
diff = bds.error - fbd.error

plot(psi,diff,pch=18,xlab="Fossil Sampling Rate",ylab="Difference in MSE (BDS - FBDRP)")

dev.off()
