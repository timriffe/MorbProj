Dat <- local(get(load("/home/tim/git/HLETTD/Data/ResultsIADL_ADL.Rdata")))
source("/home/tim/git/HLETTD/R/SurfMap.R")


pdf("/home/tim/git/MorbProj/MorbProj/Figures/SurfExample.pdf",width=8,height=4.8)
par(mai=c(.5,.5,.5,1))
SurfMap(Dat$m[["1915"]][["adl5_"]],bg=TRUE)
text(103,17,"proportion",xpd=TRUE)
dev.off()


Mat <- Dat$m[["1915"]][["adl5_"]]
Mat[row(Mat) + col(Mat) < 20] <- NA
jy1 <- rowMeans(Dat$m[["1915"]][["adl5_"]],na.rm=TRUE)
jy2 <- rowMeans(Mat,na.rm=TRUE)

pdf("/home/tim/git/MorbProj/MorbProj/Figures/LinesExample.pdf",width=6,height=4.8)
matplot(0:12, Dat$m[["1915"]][["adl5_"]], type = 'l', col = gray(.5),xlab = "years left, y", ylab = "proportion",lty=1)
lines(0:12, jy1, col = "red")
lines(0:12, jy2, col = "blue")
dev.off()

####
library(HMDHFDplus)
library(reshape2)
library(LexisUtils)
library(magrittr)

Pop       <- readHMDweb("USA","Population",username=us,password=pw)
Deaths    <- readHMDweb("USA","Deaths_1x1",username=us,password=pw)
Deathstri <- readHMDweb("USA","Deaths_lexis",username=us,password=pw)
Mx1x1     <- readHMDweb("USA","Mx_1x1",username=us,password=pw)

# some helper functions:
mx2qx <- function(mx){
	mx / (1 + .5 * mx)
}

qx2lx <- function(qx){
	cumprod(1-c(0,qx))
}

lx2dx <- function(lx){
	-diff(c(lx,0))
}




PM   <- acast(Pop, Age~Year, value.var = "Male2")
cohs <- 2014 - (74:110)
PM  <- PM[as.character(74:110),]
# period Mx from most recent year...
Mx   <- Mx1x1$Male[Mx1x1$Age >= 74 & Mx1x1$Year == 2012]

Mxi <- c(0,Mx)
cohs <- 2014 - (74:110)
ExtrapMat <- matrix(0,nrow=length(Mx),ncol=length(Mx),dimnames=list(74:110,cohs))
i <- 1
for (i in 1:length(Mx)){
	Mxi <- Mxi[-1]
	(Mxi*(.99^(1:(length(Mxi))))) %>% mx2qx %>% qx2lx %>% lx2dx -> dx
	dx[length(Mxi)] <- sum(dx[length(Mxi):length(dx)])
	dx <- dx[1:length(Mxi)]
	ExtrapMat[i:nrow(ExtrapMat),i] <- pop[i] * (dx / sum(dx))
}
ExtrapMat <- ExtrapMat[,ncol(ExtrapMat):1]
ExtrapMat <- ExtrapMat[,-ncol(ExtrapMat)]
ACD <- acast(Deathstri, Age~Cohort,sum, value.var = "Male")
ACD <- ACD[as.character(74:110), ]
ACD <- ACD[,-ncol(ACD)]
ACD <- ACD[,colSums(ACD) > 0]
ACD[,colnames(ExtrapMat)] <- ACD[,colnames(ExtrapMat)] + ExtrapMat

jyextend    <- c(jy1,rep(0,24))
ACD[ACD==0] <- NA
ACD         <- ACD[,!is.na(ACD[1,])]
#LexisMap(ACD,log=FALSE)
ACD[is.na(ACD)] <- 0


JM <- PM * 0

cohs <- as.integer(colnames(ACD))
for (a in rownames(PM)){
	for (yr in colnames(PM)){ 
		# 1) how many will die in 0 yrs?
		coh <- as.integer(yr) - as.integer(a)
		if(coh %in% cohs){
			Dx <- ACD[as.character(as.integer(a):110),as.character(coh)]
			Dx <- Dx / sum(Dx) # future dx 
			Jx <- sum(PM[a,yr] * Dx * jyextend[1:length(Dx)])
			JM[a,yr] <- Jx
		}
	}
}
LexisMap(JM / PM, log = FALSE)


Ja <- JM / PM
Ja[is.nan(Ja)] <- NA


pdf("/home/tim/git/MorbProj/MorbProj/Figures/jprimeaExample.pdf",width=6,height=4.8)
plot(74:110, Ja[,"1970"], 
		type = 'l', col = "red",xlab = "chronological age (a)", ylab = "j'(a)",lty=1,ylim=range(Ja,na.rm=TRUE))
lines(74:110, Ja[,"2010"], col = "blue")
text(80,.14,"1970",col="red")
text(85,.11,"2010",col="blue")
dev.off()


pdf("/home/tim/git/MorbProj/MorbProj/Figures/jprimeaBadExample.pdf",width=6,height=4.8)
plot(74:110, Ja[,"1970"] * PM[,"2010"], 
		type = 'l', col = "red",xlab = "chronological age (a)", ylab = "predicted J(a)",lty=1)
lines(74:110, JM[,"2010"], col = "blue")
text(82,65500,"2010 J(a) prediction from 1970 j'(a)",col="red",pos=4)
text(77,40000,"2010 J(a)",col="blue")
dev.off()

(sum(Ja[,"1970"] * PM[,"2010"],na.rm=TRUE) -
sum(JM[,"1970"],na.rm=TRUE))
# 650000 more predicted...
# 30% more than actual, in this case.
(sum(Ja[,"1970"] * PM[,"2010"],na.rm=TRUE)-
sum(JM[,"2010"],na.rm=TRUE)) / sum(JM[,"2010"],na.rm=TRUE) 


LexisMap(Jmat,log=FALSE)
PM <- PM[75:111,]
LexisMap(AC2AP(Jmat)/AC2AP(Emat),log=FALSE)



