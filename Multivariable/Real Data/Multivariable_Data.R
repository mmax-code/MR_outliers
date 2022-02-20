
setwd(".../Github/Multivariable/Real Data/")

################################################################################
# Helper Functions
################################################################################

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



# detect outliers
det.outliers = function(q.j,alpha){
  indices = as.factor(pchisq(q.j,1,lower.tail = FALSE) <= alpha/length(q.j))
  vec = rep(FALSE,length(q.j))
  vec[which(indices == TRUE)] = TRUE
  vec = as.factor(vec)
}


data = read.csv("chd_amd_alz.csv")
data = na.omit(data)
attach(data)


################################################################################
### CHD 
################################################################################

# Full model
lm = lm(chdbeta ~ ldlbeta + hdlbeta + tgbeta - 1, data = data, weights = chdse^-2)
mvmr.lm1 = summary(lm)

# Standard
q.values1 = (1/(chdse^2))*(chdbeta - (mvmr.lm1$coefficients[1]*ldlbeta +
                                       mvmr.lm1$coefficients[2]*hdlbeta + mvmr.lm1$coefficients[3]*tgbeta))^2
tmp = det.outliers(q.values1,0.05)
data$rsid[which(tmp == TRUE)]
ind = which(tmp == FALSE)
lm1 = lm(chdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = chdse[ind]^-2)
mvmr.lm1 = summary(lm1)
mvmr.lm1



# Sanderson 
sigma.aj1 = chdse^2 + mvmr.lm1$coefficients[1]*ldlse^2 + mvmr.lm1$coefficients[2]*hdlse^2 + mvmr.lm1$coefficients[3]*tgse^2
q.sanderson1 = (1/sigma.aj1)*(chdbeta - (mvmr.lm1$coefficients[1]*ldlbeta +
                                         mvmr.lm1$coefficients[2]*hdlbeta + mvmr.lm1$coefficients[3]*tgbeta))^2
tmp = det.outliers(q.sanderson1,0.05)
data$rsid[which(tmp == TRUE)]

ind = which(tmp == FALSE)
lm1 = lm(chdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = chdse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm

# MR-Presso
mr.presso = mr_presso2(BetaOutcome = "chdbeta", BetaExposure = c("ldlbeta", "hdlbeta", "tgbeta"),
                        SdOutcome = "chdse", SdExposure = c("ldlse", "hdlse", "tgse"), OUTLIERtest = TRUE,
                        DISTORTIONtest = TRUE, data = data, NbDistribution = 5000,  SignifThreshold = 0.05)

mr.presso

# Median-adjusted
lambda1 = median(q.values1)/0.456
q.adjust1 = q.values1/lambda1

tmp = det.outliers(q.adjust1,0.05)
data$rsid[which(tmp == TRUE)]

ind = which(tmp == FALSE)
lm1 = lm(chdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = chdse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm



################################################################################
### AMD
################################################################################

# Full model 

lm = lm(amdbeta ~ ldlbeta + hdlbeta + tgbeta - 1, data = data, weights = amdse^-2)
mvmr.lm2 = summary(lm)

# Standard
q.values2 = (1/(amdse^2))*(amdbeta - (mvmr.lm2$coefficients[1]*ldlbeta +
                                        mvmr.lm2$coefficients[2]*hdlbeta + mvmr.lm2$coefficients[3]*tgbeta))^2
tmp = det.outliers(q.values2,0.05)
data$rsid[which(tmp == TRUE)]

ind = which(tmp == FALSE)
lm1 = lm(amdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = amdse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm


# Sanderson 
sigma.aj2 = amdse^2 + mvmr.lm2$coefficients[1]*ldlse^2 + mvmr.lm2$coefficients[2]*hdlse^2 + mvmr.lm2$coefficients[3]*tgse^2
q.sanderson2 = (1/sigma.aj2)*(amdbeta - (mvmr.lm2$coefficients[1]*ldlbeta +
                                           mvmr.lm2$coefficients[2]*hdlbeta + mvmr.lm2$coefficients[3]*tgbeta))^2
tmp = det.outliers(q.sanderson2,0.05)
data$rsid[which(tmp == TRUE)]

ind = which(tmp == FALSE)
lm1 = lm(amdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = amdse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm


# MR-Presso
mr.presso = mr_presso2(BetaOutcome = "amdbeta", BetaExposure = c("ldlbeta", "hdlbeta", "tgbeta"),
                        SdOutcome = "amdse", SdExposure = c("ldlse", "hdlse", "tgse"), OUTLIERtest = TRUE,
                        DISTORTIONtest = TRUE, data = data, NbDistribution = 5000,  SignifThreshold = 0.05)

mr.presso

# Median-adjusted
lambda2 = median(q.values2)/0.456
q.adjust2 = q.values2/lambda2

tmp = det.outliers(q.adjust2,0.05)
data$rsid[which(tmp == TRUE)]
ind = which(tmp == FALSE)
lm1 = lm(amdbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = amdse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm

################################################################################
### ALZ
################################################################################

# Full model


lm = lm(alzbeta ~ ldlbeta + hdlbeta + tgbeta - 1, data = data, weights = alzse^-2)
mvmr.lm3 = summary(lm)
mvmr.lm3

# Standard
q.values3 = (1/(alzse^2))*(alzbeta - (mvmr.lm3$coefficients[1]*ldlbeta +
                                        mvmr.lm3$coefficients[2]*hdlbeta + mvmr.lm3$coefficients[3]*tgbeta))^2
tmp = det.outliers(q.values3,0.05)
data$rsid[which(tmp == TRUE)]
ind = which(tmp == FALSE)
lm1 = lm(alzbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = alzse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm

# Sanderson 
sigma.aj3 = alzse^2 + mvmr.lm3$coefficients[1]*ldlse^2 + mvmr.lm3$coefficients[2]*hdlse^2 + mvmr.lm3$coefficients[3]*tgse^2
q.sanderson3 = (1/sigma.aj3)*(alzbeta - (mvmr.lm3$coefficients[1]*ldlbeta +
                                           mvmr.lm3$coefficients[2]*hdlbeta + mvmr.lm3$coefficients[3]*tgbeta))^2

tmp = det.outliers(q.sanderson3,0.05)
data$rsid[which(tmp == TRUE)]
ind = which(tmp == FALSE)
lm1 = lm(alzbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = alzse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm

# MR-Presso
mr.presso = mr_presso2(BetaOutcome = "alzbeta", BetaExposure = c("ldlbeta", "hdlbeta", "tgbeta"),
                        SdOutcome = "alzse", SdExposure = c("ldlse", "hdlse", "tgse"), OUTLIERtest = TRUE,
                        DISTORTIONtest = TRUE, data = data, NbDistribution = 5000,  SignifThreshold = 0.05)
mr.presso

# Median-adjusted
lambda3 = median(q.values3)/0.456
q.adjust3 = q.values3/lambda3
tmp = det.outliers(q.adjust3,0.05)
data$rsid[which(tmp == TRUE)]
ind = which(tmp == FALSE)
lm1 = lm(alzbeta[ind] ~ ldlbeta[ind] + hdlbeta[ind] + tgbeta[ind] - 1, data = data, weights = alzse[ind]^-2)
mvmr.lm = summary(lm1)
mvmr.lm



