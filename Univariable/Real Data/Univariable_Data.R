setwd(".../Github_MVMR/Univariable/real_data/")

load("VitD.Rdata")

### Helper function

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
    mod_IVW = lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data)
    dataRandom = cbind(eval(parse(text = paste0("cbind(", paste0("rnorm(nrow(data), data[, \'", BetaExposure, "\'], data[ ,\'", SdExposure, "\'])", collapse = ", "), ", rnorm(nrow(data), fitted(mod_IVW), data[ ,\'", SdOutcome,"\']))"))), data$Weights)
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

attach(data_clump)

# Full Model
lm1 = lm(beta_ms ~ beta_vitD - 1, data = data_clump, weights = se_ms^-2)
mvmr.lm = summary(lm1)


# Standard 
q.values = (1/(se_ms^2))*(beta_ms - (mvmr.lm$coefficients[1]*beta_vitD))^2
q.outliers = det.outliers(q.values,0.05)
ind = which(q.outliers == FALSE)
lm2 = lm(beta_ms[ind] ~ beta_vitD[ind] - 1, data = data_clump, weights = se_ms[ind]^-2)
q.lm = summary(lm2)
q.lm

# Sanderson 
sigma.aj = se_ms^2 + mvmr.lm$coefficients[1]*se_vitD^2
q.sanderson = (1/sigma.aj)*(beta_ms - (mvmr.lm$coefficients[1]*beta_vitD))^2

q.outliers = det.outliers(q.sanderson,0.05)
ind = which(q.outliers == FALSE)
lm3 = lm(beta_ms[ind] ~ beta_vitD[ind] - 1, data = data_clump, weights = se_ms[ind]^-2)
sand.lm = summary(lm3)
sand.lm

# Median-adjusted
lambda = median(q.values)/(0.456)
q.adjust = q.values/lambda
q.outliers = det.outliers(q.adjust,0.05)
ind = which(q.outliers == FALSE)
lm4 = lm(beta_ms[ind] ~ beta_vitD[ind] - 1, data = data_clump, weights = se_ms[ind]^-2)
median.lm = summary(lm4)
median.lm

data_clump = as.data.frame(data_clump)
mr.presso = mr_presso2(BetaOutcome = "beta_ms", BetaExposure = "beta_vitD", SdOutcome = "se_ms", SdExposure = "se_vitD", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = data_clump, NbDistribution = 1000,  SignifThreshold = 0.05)
mr.presso