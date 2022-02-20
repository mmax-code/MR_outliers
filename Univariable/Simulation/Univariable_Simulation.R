list.of.packages = c("mvtnorm", "caret", "dplyr", "doParallel","RcppEigen",
                     "dplyr","ggplot2","scales")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd(".../Github/Univariable/Simulation")

################################################################################
### Libraries
################################################################################

library(mvtnorm)
library(caret)
library(dplyr)
library(doParallel)
library(RcppEigen)
library(dplyr)
library(ggplot2)
library(scales)

################################################################################
### Helper functions
################################################################################

# get coefficients for creation of summarized data
getcoeffs = function(y,x,name,gv){
  betas = list()
  se = list()
  #p = list()
  for(i in names(gv)){
    reg = summary(fastLm(as.matrix(cbind(rep(1,nrow(x)),get(i,x))),y))
    betas[[i]] = reg$coefficients[2,1]
    se[[i]] = reg$coefficients[2,2]
    #p[[i]] = reg$coefficients[2,4]
  }
  res = cbind(betas,se)#,p)
  colnames(res) = c(paste0("betas_",name),paste0("se_",name))#,paste0("p_",name))
  res
}

# detect outliers
det.outliers = function(q.j,alpha){
  indices = as.factor(pchisq(q.j,1,lower.tail = FALSE) <= alpha/length(q.j))
  vec = rep(FALSE,length(q.j))
  vec[which(indices == TRUE)] = TRUE
  vec = as.factor(vec)
}


# adjusted mr_presso function (original version from https://github.com/rondolab/MR-PRESSO)
mr_presso2 = function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, OUTLIERtest = FALSE, DISTORTIONtest = FALSE, SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL){
  
  if(!is.null(seed))
    set.seed(seed)
  
  if(SignifThreshold > 1)
    stop("The significance threshold cannot be greater than 1")
  
  if(length(BetaExposure) != length(SdExposure))
    stop("BetaExposure and SdExposure must have the same number of elements")
  
  if(class(data)[1] != "data.frame")
    stop("data must be an object of class data.frame, try to rerun MR-PRESSO by conversing data to a data.frame \'data = as.data.frame(data)\'")
  
  # Functions
  "%^%" = function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))
  
  getRSS_LOO = function(BetaOutcome, BetaExposure, data, returnIV){
    dataW = data[, c(BetaOutcome, BetaExposure)] * sqrt(data[, "Weights"])
    X = as.matrix(dataW[ , BetaExposure])
    Y = as.matrix(dataW[ , BetaOutcome])
    CausalEstimate_LOO = sapply(1:nrow(dataW), function(i) {
      (t(X[-i, ]) %*% X[-i, ])%^%(-1) %*% t(X[-i, ]) %*% Y[-i, ]
    })
    
    if(length(BetaExposure) == 1)
      RSS = sum((Y - CausalEstimate_LOO * X)^2, na.rm = TRUE)
    else
      RSS = sum((Y - rowSums(t(CausalEstimate_LOO) * X))^2, na.rm = TRUE)
    
    if(returnIV)
      RSS = list(RSS, CausalEstimate_LOO)
    return(RSS)
  }
  
  getRandomData = function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data){
    mod_IVW = lapply(1:nrow(data), function(i) lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)), weights = Weights, data = data[-i, ]))
    dataRandom = cbind(eval(parse(text = paste0("cbind(", paste0("rnorm(nrow(data), data[, \'", BetaExposure, "\'], data[ ,\'", SdExposure, "\'])", collapse = ", "), ", sapply(1:nrow(data), function(i) rnorm(1, predict(mod_IVW[[i]], newdata = data[i, ]), data[i ,\'", SdOutcome,"\'])))"))), data$Weights)
    colnames(dataRandom) = c(BetaExposure, BetaOutcome, "Weights")
    return(dataRandom)
  }
  
  # 0- Transforming the data + checking number of observations
  data = data[, c(BetaOutcome, BetaExposure, SdOutcome, SdExposure)]
  data = data[rowSums(is.na(data)) == 0, ]
  data[, c(BetaOutcome, BetaExposure)] = data[, c(BetaOutcome, BetaExposure)] * sign(data[, BetaExposure[1]])
  data$Weights = data$Weights = 1/data[, SdOutcome]^2
  
  if(nrow(data) <= length(BetaExposure) + 2)
    stop("Not enough intrumental variables")
  
  if(nrow(data) >= NbDistribution)
    stop("Not enough elements to compute empirical P-values, increase NbDistribution")
  
  # 1- Computing the observed residual sum of squares (RSS)
  RSSobs = getRSS_LOO(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, data = data, returnIV = OUTLIERtest)
  
  # 2- Computing the distribtion of expected residual sum of squares (RSS)
  randomData = replicate(NbDistribution, getRandomData(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, data = data), simplify = FALSE)
  RSSexp = sapply(randomData, getRSS_LOO, BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, returnIV = OUTLIERtest)
  if(OUTLIERtest)
    GlobalTest = list(RSSobs = RSSobs[[1]], Pvalue = sum(RSSexp[1, ] > RSSobs[[1]])/NbDistribution)
  else
    GlobalTest = list(RSSobs = RSSobs[[1]], Pvalue = sum(RSSexp > RSSobs[[1]])/NbDistribution)
  
  # 3- Computing the single IV outlier test
  if(GlobalTest$Pvalue < SignifThreshold & OUTLIERtest){
    OutlierTest = do.call("rbind", lapply(1:nrow(data), function(SNV){
      randomSNP = do.call("rbind", lapply(randomData, function(mat) mat[SNV, ]))
      if(length(BetaExposure) == 1){
        Dif = data[SNV, BetaOutcome] - data[SNV, BetaExposure] * RSSobs[[2]][SNV]
        Exp = randomSNP[, BetaOutcome] - randomSNP[, BetaExposure] * RSSobs[[2]][SNV]
      } else {
        Dif = data[SNV, BetaOutcome] - sum(data[SNV, BetaExposure] * RSSobs[[2]][, SNV])
        Exp = randomSNP[, BetaOutcome] - rowSums(randomSNP[, BetaExposure] * RSSobs[[2]][, SNV])
      }
      pval = sum(Exp^2 > Dif^2)/length(randomData)
      pval = cbind.data.frame(RSSobs = Dif^2, Pvalue = pval)
      return(pval)
    }))
    row.names(OutlierTest) = row.names(data)
    OutlierTest$Pvalue = apply(cbind(OutlierTest$Pvalue*nrow(data), 1), 1, min) # Bonferroni correction
  } else{
    OUTLIERtest = FALSE
  }
  
  # 4- Computing the test of the distortion of the causal estimate
  mod_all = lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data)
  if(DISTORTIONtest & OUTLIERtest){
    getRandomBias = function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, refOutlier){
      indices = c(refOutlier, replicate(nrow(data)-length(refOutlier), sample(setdiff(1:nrow(data), refOutlier))[1]))
      mod_random = lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data[indices[1:(length(indices) - length(refOutlier))], ])
      return(mod_random$coefficients[BetaExposure])
    }
    refOutlier = which(OutlierTest$Pvalue <= SignifThreshold)
    
    if(length(refOutlier) > 0){
      if(length(refOutlier) < nrow(data)){
        BiasExp = replicate(NbDistribution, getRandomBias(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, data = data, refOutlier = refOutlier), simplify = FALSE)
        BiasExp = do.call("rbind", BiasExp)
        
        mod_noOutliers = lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data[-refOutlier, ])
        BiasObs = (mod_all$coefficients[BetaExposure] - mod_noOutliers$coefficients[BetaExposure]) / abs(mod_noOutliers$coefficients[BetaExposure])
        BiasExp = (mod_all$coefficients[BetaExposure] - BiasExp) / abs(BiasExp)
        BiasTest = list(`Outliers Indices` = refOutlier, `Distortion Coefficient` = 100*BiasObs, Pvalue = sum(abs(BiasExp) > abs(BiasObs))/NbDistribution)
      } else {
        BiasTest = list(`Outliers Indices` = "All SNPs considered as outliers", `Distortion Coefficient` = NA, Pvalue = NA)
      }
    } else{
      BiasTest = list(`Outliers Indices` = "No significant outliers", `Distortion Coefficient` = NA, Pvalue = NA)
    }
  }
  
  # 5- Formatting the results
  GlobalTest$Pvalue = ifelse(GlobalTest$Pvalue == 0, paste0("<", 1/NbDistribution), GlobalTest$Pvalue)
  if(OUTLIERtest){
    OutlierTest$Pvalue = replace(OutlierTest$Pvalue, OutlierTest$Pvalue == 0, paste0("<", nrow(data)/NbDistribution))
    if(DISTORTIONtest){
      BiasTest$Pvalue = ifelse(BiasTest$Pvalue == 0, paste0("<", 1/NbDistribution), BiasTest$Pvalue)
      res = list(`Global Test` = GlobalTest, `Outlier Test` = OutlierTest, `Distortion Test` = BiasTest)
    } else {
      res = list(`Global Test` = GlobalTest, `Outlier Test` = OutlierTest)
    }
    if(nrow(data)/NbDistribution > SignifThreshold)
      warning(paste0("Outlier test unstable. The significance threshold of ", SignifThreshold, " for the outlier test is not achievable with only ", NbDistribution, " to compute the null distribution. The current precision is <", nrow(data)/NbDistribution, ". Increase NbDistribution."))
  } else {
    res = list(`Global Test` = GlobalTest)
  }
  
  OriginalMR = cbind.data.frame(BetaExposure, "Raw", summary(mod_all)$coefficients)
  colnames(OriginalMR) = c("Exposure", "MR Analysis", "Causal Estimate", "Sd", "T-stat", "P-value")
  if(exists("mod_noOutliers"))
    OutlierCorrectedMR = cbind.data.frame(BetaExposure, "Outlier-corrected", summary(mod_noOutliers)$coefficients,row.names = NULL)
  else
    OutlierCorrectedMR = cbind.data.frame(BetaExposure, "Outlier-corrected", t(rep(NA, 4)),row.names = NULL)
  colnames(OutlierCorrectedMR) = colnames(OriginalMR)
  MR = rbind.data.frame(OriginalMR, OutlierCorrectedMR)
  row.names(MR) = NULL
  
  res = list(`Main MR results` = MR, `MR-PRESSO results` = res)
  return(res)
}



################################################################################
### Simulation
################################################################################

simulation = function(seed,num.ind = 20000,var.expl1 = c(0.1,0.1), var.expl2 = 0.3, outlier.perc = 0.04,
                      ss.effects = c(1,1,1), alpha = 0.05, outlier.func){
  
  require(mvtnorm)
  require(caret)
  require(dplyr)
  require(doParallel)
  require(RcppEigen)
  require(MVMR)
  
  set.seed(seed)
  # assign("last.warning", NULL, envir = baseenv())
  
  # create gene variants
  gv = as.data.frame(replicate(100,rbinom(num.ind, 2, 0.5)))
  u = as.data.frame(replicate(100,rbinom(num.ind, 2, 0.5)))
  
  # create risk factors 
  beta1 = rnorm(100,1,2)
  beta2 = rnorm(100,1,2)
  
  # fix variance explained in first-stage regression
  var.gv = apply(gv, 2, var)
  var.u = apply(u, 2, var)
  var.eps1 = ((1 - var.expl1[1]) / var.expl1[1]) * (beta1^2 %*% var.gv)
  var.eps2 = ((1 - var.expl1[2]) / var.expl1[2]) * (beta2^2 %*% var.u)
  sigma = var.eps1
  epsilon = rnorm(num.ind, mean = 0, sd = sqrt(sigma))
  epsilon.outlier = rnorm(num.ind, 0,1)
  
  # Outlier creation
  X = as.matrix(gv)
  U = as.matrix(u)
  rho = sample(c(rep(0,ceiling((1-outlier.perc)*100)),outlier.func))
  outlier = X %*% rho + epsilon.outlier
  outlier.indices = as.factor(rho!=0)
  
  if (ss.effects[3] == 0){
    outlier.indices = rep(FALSE,100)
  }
  
  # First-stage regression
  rf1 = X %*% beta1 + epsilon
  rf2 = U %*% beta2 + rnorm(num.ind, 0,sqrt(var.eps2))
  
  RF = as.matrix(cbind(rf1,rf2,outlier))
  
  # Second-Stage regression
  # fix variance explained in second-stage regression
  var.rf = cbind(var(rf1),var(rf2),var(outlier))
  var.eps = ((1 - var.expl2) / var.expl2) * ((ss.effects^2 %*% t(var.rf)) +
                                               2*ss.effects[1]*ss.effects[3]*cov(rf1,outlier) +
                                               2*ss.effects[1]*ss.effects[2]*cov(rf1,rf2)+
                                               2*ss.effects[2]*ss.effects[3]*cov(outlier,rf2)) 
  
  
  
  epsilon = rnorm(num.ind,0,sqrt(var.eps))
  y =  RF %*% ss.effects + epsilon
  
  
  # create summarized data - two-sample summary data
  beta_rf1 = getcoeffs(rf1[1:(num.ind/2)],gv[1:(num.ind/2),],"rf1",gv = gv)
  beta_y = getcoeffs(y[((num.ind/2)+1):num.ind],gv[((num.ind/2)+1):num.ind,],"y",gv = gv)
  coeffs = as.data.frame(cbind(beta_rf1,beta_y))
  
  # transform into numeric
  coeffs[] = lapply(coeffs, function(x) {
    if(!is.numeric(x)) as.numeric(x) else x
  })
  
  # Mendelian Randomization in summarized setting
  lm1 = lm(betas_y ~ betas_rf1 - 1, data = coeffs, weights = se_y^-2)
  mvmr.lm = summary(lm1)
  mvmr.lm
  
  # simple q.values
  q.values = (1/(coeffs$se_y^2))*(coeffs$betas_y - (mvmr.lm$coefficients[1]*coeffs$betas_rf1))^2
  q.outliers = det.outliers(q.values,alpha)
  ind = which(q.outliers == FALSE)
  lm2 = lm(betas_y[ind] ~ betas_rf1[ind] - 1, data = coeffs, weights = se_y[ind]^-2)
  q.lm = summary(lm2)
  
  # Sanderson Test
  sigma.aj = coeffs$se_y^2 + mvmr.lm$coefficients[1]^2*coeffs$se_rf1^2
  q.sanderson = (1/sigma.aj)*(coeffs$betas_y - (mvmr.lm$coefficients[1]*coeffs$betas_rf1))^2
  sanderson.outliers = det.outliers(q.sanderson,alpha)
  ind = which(sanderson.outliers == FALSE)
  lm3 = lm(betas_y[ind] ~ betas_rf1[ind] - 1, data = coeffs, weights = se_y[ind]^-2)
  sanderson.lm = summary(lm3)
  
  # Adjusted q-statistics
  lambda = median(q.values)/(0.456)
  lambda[lambda<1] = 1
  q.adjust = q.values/lambda
  median.outliers = det.outliers(q.adjust,alpha)
  ind = which(median.outliers == FALSE)
  lm4 = lm(betas_y[ind] ~ betas_rf1[ind]  - 1, data = coeffs, weights = se_y[ind]^-2)
  median.lm = summary(lm4)
  
  # MR Presso test
  mr.presso = mr_presso2(BetaOutcome = "betas_y", BetaExposure = c("betas_rf1"),
                         SdOutcome = "se_y", SdExposure = c("se_rf1"), OUTLIERtest = TRUE,
                         DISTORTIONtest = TRUE, data = coeffs, NbDistribution = 1000,  SignifThreshold = alpha)
  # mr.presso = NA
  
  result = list(coeffs,lambda, ss.effects, mvmr.lm,q.values,q.outliers,q.lm,q.sanderson, sanderson.outliers, sanderson.lm,q.adjust,
                median.outliers,median.lm, mr.presso, outlier.indices,warnings(), mvmr.lm$fstatistic,lm1,lm2,lm3,lm4)
  names(result) = c("coeffs","lambda","causal.effects","model.outlier",
                    "q.values", "q.outliers", "q.lm",
                    "q.sanderson","sanderson.outliers", "sanderson.lm",
                    "q.lambda","median.outliers","median.lm",
                    "mr.presso", "outlier.indices","warnings","f.overall","lm.full","lm.q","lm.qsand","lm.qmedian")
  
  rm(list=setdiff(ls(), c("result")))
  result
}


# parallel simulation
t = Sys.time()
cl = makeCluster(50)
registerDoParallel(cl)
x1 = foreach(n = 1:1000) %dopar% simulation(seed = n, num.ind = 500000,var.expl1 = c(0.15,0.15), var.expl2 = 0.5,
                                            outlier.perc = 0.15, ss.effects = c(1,0.5,1),alpha = 0.05,
                                            runif(floor(0.15*100),5,6))


print(difftime(Sys.time(), t, units = 'sec'))

save(x1, file = "Univariable_Outliers.RData")


################################################################################
### Results
################################################################################

load("Univariable_Outliers.RData")


results.confusion = function(object){
  res = list()
  for (i in 1:1000){
    # real outliers
    outliers = as.factor(object[[i]]$outlier.indices)
    levels(outliers) = c(FALSE,TRUE)
    # simple q values
    q.outliers = as.factor(object[[i]]$q.outliers)
    levels(q.outliers) = c(FALSE,TRUE)
    tmp.q = confusionMatrix(data = q.outliers, reference = outliers, positive = "TRUE")
    # sanderson q.values
    sanderson.outliers = as.factor(object[[i]]$sanderson.outliers)
    levels(sanderson.outliers) = c(FALSE,TRUE)
    tmp.sanderson = confusionMatrix(data = sanderson.outliers, reference = outliers, positive = "TRUE")
    # median outliers
    median.outliers = as.factor(object[[i]]$median.outliers)
    levels(median.outliers) = c(FALSE,TRUE)
    tmp.median = confusionMatrix(data = median.outliers, reference = outliers, positive = "TRUE")
    # mr.presso outliers
    if (is.numeric(object[[i]]$mr.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)){
      vec = rep(FALSE,100)
      indices = object[[i]]$mr.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
      vec[indices] = TRUE
      vec = as.factor(vec)
      levels(vec) = c(FALSE,TRUE)
      tmp.presso = confusionMatrix(data = vec, reference = outliers, positive = "TRUE")
    } else{
      vec = rep(FALSE,100)
      vec = as.factor(vec)
      levels(vec) = c(FALSE,TRUE)
      tmp.presso = confusionMatrix(data = vec, reference = outliers, positive = "TRUE")
    }
    tmp = list(tmp.q,tmp.sanderson,tmp.median,tmp.presso)
    names(tmp) = c("conf.q","conf.sanderson","conf.median","conf.presso")
    res[[i]] = tmp
  }
  res
}



aggregate.results = function(confusion.object,method,name){
  overall.sens = NULL
  overall.spec = NULL
  for (i in 1:1000){
    # Sensitivity
    tmp1 = confusion.object[[i]][[method]]$byClass[1]
    # Specificity
    tmp2 = confusion.object[[i]][[method]]$byClass[2]
    overall.sens = c(overall.sens,tmp1)
    overall.spec = c(overall.spec,tmp2)
  }
  res1 = mean(overall.sens,na.rm=TRUE)
  res2 = mean(overall.spec,na.rm=TRUE)
  res.complete = list(res1,res2)
  names(res.complete) = c(paste0("sensitivity ",name), paste0("specificity ",name))
  res.complete
}



result = results.confusion(x1)
q = aggregate.results(result,"conf.q","normal.q")
sanderson = aggregate.results(result,"conf.sanderson","sanderson.q")
median = aggregate.results(result,"conf.median","median.q")
presso = aggregate.results(result,"conf.presso" ,"presso.q")

q;sanderson;median;presso



# Mean Bias
full = list()
lm.q = list()
lm.sanderson = list()
lm.median = list()
lm.presso = list()
for (i in 1:1000){
  full[[i]] = c((x1[[i]]$lm.full$coefficients[1]-1))
  lm.q[[i]] = c((x1[[i]]$lm.q$coefficients[1]-1))
  lm.sanderson[[i]] = c((x1[[i]]$lm.qsand$coefficients[1]-1))
  lm.median[[i]] = c((x1[[i]]$lm.qmedian$coefficients[1]-1))
  lm.presso[[i]] = c((x1[[i]]$mr.presso$`Main MR results`$`Causal Estimate`[2]-1))
}
colMeans(do.call(rbind,full),na.rm = TRUE)
colMeans(do.call(rbind,lm.q),na.rm = TRUE)
colMeans(do.call(rbind,lm.sanderson),na.rm = TRUE)
colMeans(do.call(rbind,lm.median),na.rm = TRUE)
colMeans(do.call(rbind,lm.presso),na.rm = TRUE)


# MSE
full = list()
lm.q = list()
lm.sanderson = list()
lm.median = list()
lm.presso = list()
for (i in 1:1000){
  full[[i]] = c((x1[[i]]$lm.full$coefficients[1]-1)^2) 
  lm.q[[i]] = c((x1[[i]]$lm.q$coefficients[1]-1)^2)
  lm.sanderson[[i]] = c((x1[[i]]$lm.qsand$coefficients[1]-1)^2)
  lm.median[[i]] = c((x1[[i]]$lm.qmedian$coefficients[1]-1)^2)
  lm.presso[[i]] = c((x1[[i]]$mr.presso$`Main MR results`$`Causal Estimate`[1]-1)^2)
}

colMeans(do.call(rbind,full),na.rm = TRUE)
colMeans(do.call(rbind,lm.q),na.rm = TRUE)
colMeans(do.call(rbind,lm.sanderson),na.rm = TRUE)
colMeans(do.call(rbind,lm.median),na.rm = TRUE)
colMeans(do.call(rbind,lm.presso),na.rm = TRUE)


# Average number of detected outliers
outlier.average = function(object){
  res = list()
  for (i in 1:1000){
    # simple q values
    tmp.q = as.numeric(table(object[[i]]$q.outliers)[2])
    # sanderson q.values
    tmp.sanderson = as.numeric(table(object[[i]]$sanderson.outliers)[2])
    # median outliers
    tmp.median = as.numeric(table(object[[i]]$median.outliers)[2])
    # mr.presso outliers
    if (is.numeric(object[[i]]$mr.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)){
      tmp.presso= length(object[[i]]$mr.presso$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)
    } else{
      tmp.presso = 0
    }
    tmp = list(tmp.q,tmp.sanderson,tmp.median,tmp.presso)
    names(tmp) = c("q","sanderson","median","presso")
    res[[i]] = tmp
  }
  res
}

outliers = outlier.average(x1)
tmp = do.call(rbind, lapply(outliers, as.data.frame))
colMeans(tmp)/15
colMeans(tmp)




rm(list = ls())
.rs.restartR()
