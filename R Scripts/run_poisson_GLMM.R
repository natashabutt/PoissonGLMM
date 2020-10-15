# *** SCRIPT FOR RUNNING POISSON GLMM MODEL IN BUGS ***

# *** Imports and Environment ***

library(R2WinBUGS)
library(dplyr)
library(latex2exp)
library(lattice)

# set working directory
setwd("D:/PoissonGLMM")

# *** Main ***

force_data <- read.csv('Police Data/force_data.csv')%>%
  filter(ethnicity!='None'&ethnicity!='Other')

nForces <- length(unique(force_data$PFA16NM))
nEth <- length(unique(force_data$ethnicity))

datalist <- list(
  nForces = nForces,
  nEth = nEth,
  S = structure(
    .Data = force_data$stops,
    .Dim = c(nEth,nForces)),
  A = structure(
    .Data = force_data$arrests,
    .Dim = c(nEth,nForces))
)

#bugs.data(datalist, dir = getwd(), digits = 5, data.file = "BUGS Models/Random Effects/data.txt")

inits <- function(){
  list(
    sigma_e = runif(1,0.001,100),
    sigma_b = runif(1,0.001,100),
    alpha = rnorm(1),
    Mu = rnorm(nEth-1),
    beta = rnorm(nForces),
    e = structure(
      .Data = rnorm(nEth*nForces),
      .Dim =c(nEth,nForces))
  )}
inits()

#bugs.data(inits(), dir = getwd(), digits = 5, data.file = "BUGS Models/Random Effects/inits.txt")

parameters <- c("alpha","mu","sigma_b","sigma_e")

num_iter = 2000000

model_4.1 <- bugs(model.file= "D:/PoissonGLMM/BUGS Models/Random Effects/Poisson_GLMM.odc",
                     data = datalist,
                     parameters = parameters,
                     inits = inits,
                     debug = FALSE,
                     DIC=FALSE,
                     n.chains = 3,
                     n.iter = num_iter, n.burnin = round(num_iter/3), n.thin = 10,
                     bugs.directory = "C:/Users/nbutt/Desktop/winbugs/WinBUGS14")

print(model_4.1)
with_burn <- window(as.mcmc.list(model_4.1),start=num_iter/2)

write.csv(gelman.diag(with_burn)[1],paste('BUGS Models/Random Effects/All/',as.integer(num_iter),'/rhat.csv',sep=""))
write.csv(HPDinterval(with_burn,prob=0.95)[1],paste('BUGS Models/Random Effects/All/',as.integer(num_iter),'/hpd.csv',sep=""))
write.csv(summary(with_burn)$statistics,paste('BUGS Models/Random Effects/All/',as.integer(num_iter),'/summary.csv',sep=""))
write.csv(effectiveSize(with_burn),paste('BUGS Models/Random Effects/All/',as.integer(num_iter),'/neff.csv',sep=""))

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/alpha.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'alpha'],smooth=FALSE)
densplot(with_burn[,'alpha'],smooth=FALSE)
mtext(TeX('$\\alpha$'), outer = TRUE, cex = 1.5)
dev.off()

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/mu_1.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'mu[1]'],smooth=FALSE)
densplot(with_burn[,'mu[1]'],smooth=FALSE)
mtext(TeX('$\\mu_{asian}$'), outer = TRUE, cex = 1.5)
dev.off()

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/mu_2.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'mu[2]'],smooth=FALSE)
densplot(with_burn[,'mu[2]'],smooth=FALSE)
mtext(TeX('$\\mu_{black}$'), outer = TRUE, cex = 1.5)
dev.off()

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/mu_3.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'mu[3]'],smooth=FALSE)
densplot(with_burn[,'mu[3]'],smooth=FALSE)
mtext(TeX('$\\mu_{mixed}$'), outer = TRUE, cex = 1.5)
dev.off()

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/sigma_b.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'sigma_b'],smooth=FALSE)
densplot(with_burn[,'sigma_b'],smooth=FALSE)
mtext(TeX('$\\sigma_{\\beta}$'), outer = TRUE, cex = 1.5)
dev.off()

jpeg(paste("BUGS Models/Random Effects/All/",as.integer(num_iter),"/sigma_e.jpg",sep=''), width = 450, height = 180)
par(mfrow=c(1,2),oma = c(0, 0, 2, 0),mar=c(4,3,1,2))
traceplot(with_burn[,'sigma_e'],smooth=FALSE)
densplot(with_burn[,'sigma_e'],smooth=FALSE)
mtext(TeX('$\\sigma_{\\epsilon}$'), outer = TRUE, cex = 1.5)
dev.off()
