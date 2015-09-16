Dat <- local(get(load("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/Data/ResultsIADL_ADL.Rdata")))
source("/home/tim/git/ATL_SexDecomp/ATL_SexDecomp/R/SurfMap.R")


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






library(HMDHFDplus)
Pop    <- readHMDweb("USA","Population",username=us,password=pw)
Deaths <- readHMDweb("USA","Deaths_1x1",username=us,password=pw)
Deathstri <- readHMDweb("USA","Deaths_lexis",username=us,password=pw)

library(reshape2)
DM <- acast(Deaths,Age~Year,value.var = "Male")

DM[as.character(70:110),"2000"]



library(LexisUtils)
LexisMap(DM,log=FALSE)

cMxm <- c(0.049986, 0.052877, 0.057103, 0.061325, 0.066315, 0.07132, 
		0.081139, 0.088158, 0.096355, 0.104718, 0.113543, 0.122644, 0.133971, 
		0.144406, 0.15709, 0.169257, 0.183678, 0.199131, 0.217271, 0.239125, 
		0.263696, 0.285049)
#dput(cMxf)
cMxf <- c(0.028989, 0.031342, 0.034459, 0.037594, 0.041088, 0.045432, 
		0.05254, 0.057924, 0.064801, 0.071759, 0.079933, 0.088454, 0.09763, 
		0.107266, 0.118591, 0.129736, 0.142486, 0.157165, 0.173349, 0.191498, 
		0.212159, 0.237106)

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

qxm <- mx2qx(cMxm)
qxf <- mx2qx(cMxf)
lxm <- qx2lx(qxm)
lxf <- qx2lx(qxf)
# these could of course be rescaled such that the radix is the 1992 population of people from
# the 1915-1919 cohorts. But I think it's not necessary here.
#plot(74:96,lxm,type='l')
#lines(74:96,lxf,col = "blue")

# dx is the real sugar for working with these data.
plot(lx2dx(lxm),type='l')
lines(lx2dx(lxf),col="red")

############################################################################
# note: the people alive in x are the people that die in x+

# we learn that there is another direction in which we ought to extend
# the data. Since we only have up to age 95, there are people alive in l(x)
# that we don't have in this window of d(x)...

# That last open age group of d(x)
# is therefore present in each preceding age of l(x), and if we want to see the
# chrono consequences of these morbidity patterns then we need to be able to 
# compose an l(x) out of it's various d(x) components, taking into account
# the thanatological trajectories of each d(x). So we need values for those 
# missing d(x), for the ages > 95. I'm inclined to just repeat the experience of those
# that died at age 95, shifting it right. Then the matrix will get bigger.
# another problem will be that we need to extend cMx until the cohort
# is mostly extinct. In that case, how about doing the HMD extendo-O-matic
# thing and borrowing from the preceding cohorts.

# in short, we need to 'close out' both the mortality data and the morbidity data. 
# In the end, this will have little leverage on results. If mortality were
# a lot lower, it would have leverage. Realistically, but the time mortality
# ever does get that low, we'll actually have the data and not need to extrapolate!!


coefm <- lm(log(cMxm)~c(74:95))$coef
coeff <- lm(log(cMxf)~c(74:95))$coef

# predict ages 96-110...
cMxme <- exp(coefm[1] + coefm[2] * 96:110)
cMxfe <- exp(coeff[1] + coeff[2] * 96:110)

# looks good.
#plot(74:95, cMxm,log='y', type = 'l', col = "blue",xlim=c(74,110),ylim=c(.02,1))
#lines(74:95, cMxf,col="red")
#points(96:110,cMxme,pch=19,col="blue" )
#points(96:110,cMxfe,pch=19,col="red" )

# so we'll redo to get the final d(x)
cMxmE <- c(cMxm,cMxme)
cMxfE <- c(cMxf,cMxfe)

lxm <- qx2lx(mx2qx(cMxmE))
lxf <- qx2lx(mx2qx(cMxfE))

dxm <- lx2dx(lxm)
dxf <- lx2dx(lxf)

# groups the last two dx elements so everything both adds and is conformable
dxm <- c(dxm[1:(length(dxm)-2)],sum(rev(dxm)[1:2]))
dxf <- c(dxf[1:(length(dxf)-2)],sum(rev(dxf)[1:2]))
# cool, we're closed out now.
#plot(dxm)

length(74:110)
length(dxm)



PM <- acast(Pop, Age~Year, value.var = "Male2")

cohs <- 2014 - (74:110)
pop  <- PM[as.character(74:110),"2013"]

Mx1x1 <- readHMDweb("USA","Mx_1x1",username=us,password=pw)
Mx <- Mx1x1$Male[Mx1x1$Age >= 74 & Mx1x1$Year == 2012]

library(magrittr)
#Mx %>% mx2qx %>% qx2lx %>% lx2dx -> dx
#dxend <- dx[length(dx)] + dx[length(dx)-1]
#dx <- c(dx[1:36],dxend)
#extrapolating deaths out from 74+ pop
length(2013:(2013+length(Mx)))
cohs <- 2014 - (74:110)
ExtrapMat <- matrix(0,nrow=length(Mx),ncol=length(Mx),dimnames=list(74:110,cohs))
for (i in 1:length(Mx)){
	(Mx*.98^i) %>% mx2qx %>% qx2lx %>% lx2dx -> dx
	ExtrapMat[i:nrow(ExtrapMat),i] <- pop[i] * (dx[i:37]) / sum((dx[i:37]))
}
ExtrapMat <- ExtrapMat[,ncol(ExtrapMat):1]
ExtrapMat <- ExtrapMat[,-ncol(ExtrapMat)]
ACD <- acast(Deathstri, Age~Cohort,sum, value.var = "Male")
ACD <- ACD[as.character(74:110), ]
ACD <- ACD[,-ncol(ACD)]
ACD <- ACD[,colSums(ACD) > 0]
ACD[,colnames(ExtrapMat)] <- ACD[,colnames(ExtrapMat)] + ExtrapMat



jyextend <- c(jy1,rep(0,24))
ACD[ACD==0] <- NA
ACD <- ACD[,!is.na(ACD[1,])]
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

(sum(Ja[,"1970"] * PM[,"2010"],na.rm=TRUE)-
sum(JM[,"2010"],na.rm=TRUE)) / sum(JM[,"2010"],na.rm=TRUE) 


LexisMap(Jmat,log=FALSE)
PM <- PM[75:111,]
LexisMap(AC2AP(Jmat)/AC2AP(Emat),log=FALSE)



