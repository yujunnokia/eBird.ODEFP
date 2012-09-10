# TODO: Learn the model parameter from the synthetic data 
#
# Author: Jun Yu
# Version: Jan 20, 2012
##################################################################

rm(list=ls())

#setwd("C:/Jun_home/workspace/eBird.ODEFP")
setwd("/Users/yujunnokia/Documents/workspace/eBird.ODEFP")
source("ODEFP.R")
source("ODEFP.synthData.R")

library("lattice")
library("Matrix")
library("glmnet")

######################
# experiment settings
######################
nTrSites <- 500  # number of training sites
nTeSites <- 500  # number of testing sites
nVisits <- 3  # number of visits to each site
nObservers <- 50  # number of observers
nOccCovs <- 5  # number of occupancy covariates
nDetCovs <- 5  # number of detection covariates
nExpCovs <- 5  # number of expertise covariates
nParams <- nOccCovs + nDetCovs * 3 + nExpCovs + 5  # total number of paramters.
nRandomRestarts <- 1  # number of random restarts  

#################
# regularization
#################
regType <- 2 # regularization types: 0 for none, 1 for L1, 2 for L2
lambda <- lambdaO <- lambdaD <- lambdaE <- 0.01  # regularization paramters

#######################
# Set model parameters
#######################
alpha <- c(1,rnorm(nOccCovs)*3)
beta <- c(1,rnorm(nDetCovs, mean = 1.2, sd = 1)*3)
kappa <- c(1,rnorm(nDetCovs, mean = -2.0, sd = 1)*3)
gamma <- c(1,rnorm(nDetCovs, mean = 0.6, sd = 1)*3)
eta <- c(1,rnorm(nDetCovs, mean = -1, sd = 1)*3)
nu <- c(1,rnorm(nExpCovs)*3)
covs <- rnorm(1000*nDetCovs, mean = 0.5, sd = 1)
dim(covs) <- c(1000,nDetCovs)
covs <- cbind(1,covs)
cat("mean Ex true detection prob is",mean(Logistic(covs %*% beta)),"\n",sep=" ")
cat("mean Ex false detection prob is",mean(Logistic(covs %*% kappa)),"\n",sep=" ")
cat("mean No true detection prob is",mean(Logistic(covs %*% gamma)),"\n",sep=" ")
cat("mean No false detection prob is",mean(Logistic(covs %*% eta)),"\n",sep=" ")

########################
# generate testing data
########################
teVisits <- array(nVisits, nTeSites)
teData <- GenerateData(nTeSites,teVisits,nObservers,alpha,beta,gamma,eta,nu)
teDetHists <- teData$detHists
teObservers <- teData$observers
teExpertise <- teData$expertise
teOccCovs <- teData$occCovs
teDetCovs <- teData$detCovs
teExpCovs <- teData$expCovs
teTrueOccs <- teData$trueOccs

############################
# generate training data
############################
trVisits <- array(0, c(nTrSites,1))
for (i in 1:nTrSites) {
    isMultipleVisits <- runif(1) < 0.5
    if (isMultipleVisits == TRUE) {
        trVisits[i] <- round(runif(1, min=2, max=nVisits))
    } else {
        trVisits[i] <- 1
    }
}
trData <- GenerateData(nTrSites,trVisits,nObservers,alpha,beta,gamma,eta,nu)
trDetHists <- trData$detHists
trObservers <- trData$observers
trExpertise <- trData$expertise
trOccCovs <- trData$occCovs
trDetCovs <- trData$detCovs
trExpCovs <- trData$expCovs
trTrueOccs <- trData$trueOccs

#####################
# get Bayes rates
#####################
{
    # get occupancy rate and detection rate
    teOccProb <- array(0,c(nTeSites,1))
    teExTrueDetProb  <- array(0,c(nTeSites,nVisits))
    teExFalseDetProb <- array(0,c(nTeSites,nVisits))
    teNoTrueDetProb  <- array(0,c(nTeSites,nVisits))
    teNoFalseDetProb <- array(0,c(nTeSites,nVisits))
    teExpProb <- array(0,c(nObservers,1))
    tePredExpertise <- array(0,c(nObservers,1))
    predDetHists <- array(0,c(nTeSites,nVisits))
    
    teOccProb <- Logistic(teOccCovs %*% alpha)
    teExpProb <- Logistic(teExpCovs %*% nu)
    tePredExpertise <- round(teExpProb)
    
    for (i in 1:nTeSites) {
        for (t in 1:teVisits[i]) {
            teExTrueDetProb[i,t]  <- Logistic(teDetCovs[i,t,] %*% beta)
            teExFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% kappa)
            teNoTrueDetProb[i,t]  <- Logistic(teDetCovs[i,t,] %*% gamma)
            teNoFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% eta)
            
            if (tePredExpertise[teObservers[i,t]] == 1) {
                if (round(teOccProb[i]) == 1) {
                    predDetHists[i,t] <- round(teExTrueDetProb[i,t])
                } else {
                    predDetHists[i,t] <- round(teExFalseDetProb[i,t])
                }
            } else {
                if (round(teOccProb[i]) == 1) {
                    predDetHists[i,t] <- round(teNoTrueDetProb[i,t])
                } else {
                    predDetHists[i,t] <- round(teNoFalseDetProb[i,t])				
                }
            }
        } # t
    } # i
    bayesOcc <- sum(round(teOccProb) == teTrueOccs)  / nTeSites
    bayesExp <- sum(round(teExpProb) == teExpertise) / nObservers
    bayesDet <- sum(sum(predDetHists == teDetHists)) / (sum(teVisits))
    cat("bayes occupancy rate is ",bayesOcc,"\n")
    cat("bayes expertise rate is ",bayesExp,"\n")
    cat("bayes detection rate is ",bayesDet,"\n")
}

########
# ODEFP
########
{
    # run ODE
    params <- RandomRestartEM(trDetHists,trObservers,trExpertise,trOccCovs,trDetCovs,
            trExpCovs,trVisits,regType,lambdaO,lambdaD,lambdaE,nRandomRestarts)
    alphaODEFP <- params$alpha
    betaODEFP <- params$beta
    kappaODEFP <- params$kappa
    gammaODEFP <- params$gamma
    etaODEFP <- params$eta
    nuODEFP <- params$nu
    
    # get occupancy rate and detection rate
    teOccProb <- array(0,c(nTeSites,1))
    teExTrueDetProb <- array(0,c(nTeSites,nVisits))
    teExFalseDetProb <- array(0,c(nTeSites,nVisits))
    teNoTrueDetProb <- array(0,c(nTeSites,nVisits))
    teNoFalseDetProb <- array(0,c(nTeSites,nVisits))
    teExpProb <- array(0,c(nObservers,1))
    tePredExpertise <- array(0,c(nObservers,1))
    predDetHists <- array(0,c(nTeSites,nVisits))
    
    teOccProb <- Logistic(teOccCovs %*% alphaODEFP)
    
    teExpProb <- Logistic(teExpCovs %*% nuODEFP)
    tePredExpertise <- round(teExpProb)
    
    for (i in 1:nTeSites) {
        for (t in 1:teVisits[i]) {
            teExTrueDetProb[i,t]  <- Logistic(teDetCovs[i,t,] %*% betaODEFP)
            teExFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% kappaODEFP)
            teNoTrueDetProb[i,t]  <- Logistic(teDetCovs[i,t,] %*% gammaODEFP)
            teNoFalseDetProb[i,t] <- Logistic(teDetCovs[i,t,] %*% etaODEFP)
            
            if (tePredExpertise[teObservers[i,t]] == 1)  {
                if (round(teOccProb[i]) == 1) {
                    predDetHists[i,t] <- round(teExTrueDetProb[i,t])
                } else {
                    predDetHists[i,t] <- round(teExFalseDetProb[i,t])
                }
            } else {
                if (round(teOccProb[i]) == 1) {
                    predDetHists[i,t] <- round(teNoTrueDetProb[i,t])
                } else {
                    predDetHists[i,t] <- round(teNoFalseDetProb[i,t])				
                }
            }
        } # t
    } # i
    modelOcc <- sum(round(teOccProb) == teTrueOccs)  / nTeSites
    modelExp <- sum(round(teExpProb) == teExpertise) / nObservers
    modelDet <- sum(sum(predDetHists == teDetHists)) / (sum(teVisits))
    cat("------------------------------\n")
    cat("bayes occupancy rate is ", bayesOcc, "\n")
    cat("model occupancy rate is ", modelOcc, "\n")
    cat("bayes expertise rate is ", bayesExp, "\n")
    cat("model expertise rate is ", modelExp, "\n")
    cat("bayes detection rate is ", bayesDet, "\n")
    cat("model detection rate is ", modelDet, "\n")
    
    # predict Z on test data
    trueOccProb <- array(0,c(nTeSites,1))
    modelOccProb <- array(0,c(nTeSites,1))
    for (i in 1:nTeSites) {
        trueOccProb[i] <- PredictOcc(c(alpha,beta,kappa,gamma,eta,nu),
                teOccCovs[i,],teDetCovs[i,,],teExpCovs,teDetHists[i,],teObservers[i,],teVisits[i]) 
        modelOccProb[i] <- PredictOcc(c(alphaODEFP,betaODEFP,kappaODEFP,gammaODEFP,etaODEFP,nuODEFP),
                teOccCovs[i,],teDetCovs[i,,],teExpCovs,teDetHists[i,],teObservers[i,],teVisits[i]) 
    }
    trueOcc <- sum(round(trueOccProb)==teTrueOccs)/nTeSites
    predOcc <- sum(round(modelOccProb)==teTrueOccs)/nTeSites
    cat("------------------------------\n")
    cat("True occupancy prediction is ",trueOcc,"\n")
    cat("Model occupancy prediction is ",predOcc,"\n")
    
    # predict Y on test data
    trueDetHists <- array(0,c(nTeSites,nVisits))
    modelDetHists <- array(0,c(nTeSites,nVisits))
    for (i in 1:nTeSites) {
        for (t in 1:teVisits[i]) {
            trueDetHists[i,t] <- PredictDet(c(alpha,beta,kappa,gamma,eta,nu),
                    teOccCovs[i,],teDetCovs[i,t,],teExpCovs,teObservers[i,t]) 
            modelDetHists[i,t] <- PredictDet(c(alphaODEFP,betaODEFP,kappaODEFP,gammaODEFP,etaODEFP,nuODEFP),
                    teOccCovs[i,],teDetCovs[i,t,],teExpCovs,teObservers[i,t]) 
        }
    }
    trueDet <- sum(sum(round(trueDetHists) == teDetHists)) / (sum(teVisits))
    predDet <- sum(sum(round(modelDetHists) == teDetHists)) / (sum(teVisits))
    cat("True detection prediction is ",trueDet,"\n")
    cat("Model detection prediction is ",predDet,"\n")
    
    # predict E on test data
    trueExpProb <- array(0,c(nObservers,1))
    modelExpProb <- array(0,c(nObservers,1))
    for (j in 1:nObservers) {
        trueExpProb[j] <- PredictExp(c(alpha,beta,kappa,gamma,eta,nu),
                teOccCovs,teDetCovs,teExpCovs,teDetHists,teVisits,teObservers,j) 
        modelExpProb[j] <- PredictExp(c(alphaODEFP,betaODEFP,kappaODEFP,gammaODEFP,etaODEFP,nuODEFP),
                teOccCovs,teDetCovs,teExpCovs,teDetHists,teVisits,teObservers,j) 
    }
    trueExp <- sum(round(trueExpProb) == teExpertise) / nObservers
    predExp <- sum(round(modelExpProb) == teExpertise) / nObservers
    cat("True expertise prediction is ",trueExp,"\n")
    cat("Model expertise prediction is ",predExp,"\n")
    
    
    # compute MSE
    MSE <- sum(sum((alpha-alphaODEFP)^2) +  sum((beta-betaODEFP)^2) + sum((kappa-kappaODEFP)^2) + 
                    sum((gamma-gammaODEFP)^2) + sum((eta-etaODEFP)^2) +  sum((nu-nuODEFP)^2)) / nParams
    cat("------------------------------\n")
    cat("MSE is",MSE,"\n")
}


##################################
## evaluate on testing data on Y
##################################
## evaluate on testing data
#predictionY <- array(0, c(nTeSites,5))
#for (i in 1:nTeSites) {
#	id <- teObservers[i,1]
#	teOccCov <- teOccCovs[i,]
#	teDetCov <- teDetCovs[i,1,]
#	teExpCov <- teExpCovs[id,]
#	expCov <- Logistic(t(delta)%*%teExpCov)
#	teCov <- c(1,teOccCovs[i,2:(nOccCovs+1)],teDetCovs[i,1,2:(nDetCovs+1)],expCov)
#	teDetCovEOM <- c(teDetCovs[i,1,],expCov)
#	teDetCovEONM <- teDetCovEOM
#	
#	predictionY[i,1] <- teDetHists[i,1]
#	
#	# LRM
#	predictionY[i,2] <- exp(t(logitCoeffY) %*% teCov - logsumexp(0,(t(logitCoeffY) %*% teCov)))
#	
#	# EOM
#	Pz <- exp(t(alphaEOM) %*% teOccCov - logsumexp(0,(t(alphaEOM) %*% teOccCov)))
#	Pd <- exp(t(betaEOM) %*% teDetCovEOM - logsumexp(0,(t(betaEOM) %*% teDetCovEOM)))
#	predictionY[i,3] <- Pz*Pd
#	
#	# EOM.FP
#	Pz <- exp(t(alphaEONM) %*% teOccCov - logsumexp(0,(t(alphaEONM) %*% teOccCov)))
#	Pdg <- exp(t(gammaEONM) %*% teDetCovEONM - logsumexp(0,(t(gammaEONM) %*% teDetCovEONM)))
#	Pde <- exp(t(etaEONM) %*% teDetCovEONM - logsumexp(0,(t(etaEONM) %*% teDetCovEONM)))
#	predictionY[i,4] <- Pz*Pdg+(1-Pz)*Pde
#	
#	# EOEM
#	Pz <- exp(t(alphaEOEM) %*% teOccCov - logsumexp(0,(t(alphaEOEM) %*% teOccCov)))
#	Pdb <- exp(t(betaEOEM) %*% teDetCov - logsumexp(0,(t(betaEOEM) %*% teDetCov)))
#	Pdg <- exp(t(gammaEOEM) %*% teDetCov - logsumexp(0,(t(gammaEOEM) %*% teDetCov)))
#	Pde <- exp(t(etaEOEM) %*% teDetCov - logsumexp(0,(t(etaEOEM) %*% teDetCov)))
#	Pe <- exp(t(nuEOEM) %*% teExpCov - logsumexp(0,(t(nuEOEM) %*% teExpCov)))
#	predictionY[i,5] <- Pe*Pz*Pdb+(1-Pe)*Pz*Pdg+(1-Pe)*(1-Pz)*Pde
#} # i
#
######################################
## generate and save ROC and AUC on Y
######################################
#roc <- genROC(predictionY)
#rocfile <- "../data/synthetic/testEOEM.Y.csv"
#rocColNames <- c("Threshold","LRM_FPR","LRM_TPR","EOM_FPR","EOM_TPR","EOM.FP_FPR","EOM.FP_TPR","EOEM.FP_FPR","EOEM.FP_TPR")
#write.table(roc,file=rocfile,row.names=FALSE,col.names=rocColNames,sep=",")
#auc <- genAUC(rocfile,1)
#print ("EOEM synthetic data results on Y:")
#print ("AUC")
#print (auc)
#
#################################
## evaluate on testing data on Z
#################################
#predictionZ <- array(0, c(nTeSites,5))
#teDetCovEOM <- NULL
#for (i in 1:nTeSites) {
#	id <- teObservers[i,1]
#	teExpCov <- teExpCovs[id,]
#	expCov <- Logistic(t(delta)%*%teExpCov)
#	teCov <- c(1,teOccCovs[i,2:(nOccCovs+1)],teDetCovs[i,1,2:(nDetCovs+1)],expCov,teDetHists[i,1])
#	predictionZ[i,1] <- teTrueOccs[i]
#	teDetCovEOM <- rbind(teDetCovEOM,c(teDetCovs[i,1,],expCov))
#	
#	
#	# LRM
#	predictionZ[i,2] <- exp(t(logitCoeffZ) %*% teCov - logsumexp(0,(t(logitCoeffZ) %*% teCov)))
#}
#dim(teDetCovEOM) <- c(nTeSites,1,nDetCovs+2)
#teDetCovEONM <- teDetCovEOM
#
## evaluate EOM
#source("../eBirdR.EOM/utils.R")
#source("../eBirdR.EOM/logreg.R")
#predictionZ[,3] <- doEstep(c(alphaEOM,betaEOM),teDetHists,teOccCovs,teDetCovEOM,teVisits)
#
## evaluate EOM.FP
#source("../eBirdR.EOM.FP/utils.R")
#source("../eBirdR.EOM.FP/logreg.R")
#predictionZ[,4] <- doEstep(c(alphaEONM,gammaEONM,etaEONM),as.matrix(teDetHists[,1]),teOccCovs,teDetCovEOM,teVisits)
#
## evaluate EOEM
#source("../eBirdR.EOEM/utils.R")
#source("../eBirdR.EOEM/logreg.R")
#predictionZ[,5] <- doEstep(c(alphaEOEM,betaEOEM,gammaEOEM,etaEOEM,nuEOEM),as.matrix(teDetHists[,1]),teObservers,teExpertise,teOccCovs,teDetCovs,teExpCovs,teVisits)
#
######################################
## generate and save ROC and AUC on Z
######################################
#roc <- genROC(predictionZ)
#rocfile <- "../data/synthetic/testEOEM.Z.csv"
#rocColNames <- c("Threshold","LRM_FPR","LRM_TPR","EOM_FPR","EOM_TPR","EOM.FP_FPR","EOM.FP_TPR","EOEM.FP_FPR","EOEM.FP_TPR")
#write.table(roc,file=rocfile,row.names=FALSE,col.names=rocColNames,sep=",")
#auc <- genAUC(rocfile,1)
#print ("EOEM synthetic data results on Z:")
#print ("AUC")
#print (auc)
