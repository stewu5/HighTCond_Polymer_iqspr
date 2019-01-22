#/ Sample code for generating candidates of high thermal conductivity polymers
#/ Written by Stephen Wu, last update @ 2018.05.26
#/ Note 1: Step 1-5 are commented out due to data closure limitation.
#/         Instead, a pre-trained model is loaded to skip those steps.
#/ Note 2: Here we prepared a fixed set of initial samples due to data closure limitation,
#/         but in our original run, we took random sample from PoLyInfo.

library(iqspr)

### 1. Set working directory for the R codes
# setwd(...)

### 2. load PoLyInfo property data and SMILES strings into "prop" and "smis", respectively
# smis <- ...
# prop <- ...

### 3. set up corresponding fingerprints as input
# fpname <- c(
#   "standard",
#   "extended",
#   "hybridization",
#   "maccs",
#   "circular",
#   "pubchem"
# )  

### 4. get forward model
# qsprpred_env <- QSPRpred()
# qsprpred_env$init_env(smis=smis, prop=prop, v_fnames=fpname)
# qsprpred_env$model_training(model=c("linear_Bayes"),params=NA)

### 5. train ENgram prior model
# engram_poly <- ENgram$new(all.polyinfo.smiles, order = 10)

### [load pre-trained models to skip step 1-5]
load("models.RData")
load("linear_Bayes_params_1.Rda")
load("linear_Bayes_params_2.Rda")

### 6. parameters for iqspr & output variables
para <- list(Tg.min = 200,
             n.sample = 100,
             n.iter = 500,
             n.smis = 1000)
gensmis <- vector("list",length = para$n.iter)
iqspr.pred <- vector("list",length = para$n.iter)

### 7. setup target space
targ.min <- c(para$Tg.min,para$Tg.min+100)
targ.max <- c(para$Tg.min+200,para$Tg.min+100+200)
qsprpred_env$set_target(targ.min,targ.max)

### 8. prepare variables LCP fragment screening
tmp <- read.csv("LC_frags.csv")
frags <- as.character(tmp[,"SMILES"]) #/pick SMILES to match for now

### 9. inverse problem 
# Note: due to data sharing restriction of PoLyInfo, initial structures are selected from the Polymer Genome database
load("initial_SMILES.RData")
smis0 <- sample(smis.uni, para$n.sample)
burnin = 100
temp = c(seq(30, 1, length = burnin), rep(1, para$n.iter - burnin))
for(i in 1:para$n.iter){
  #/run iqspr with best selected molecules from previous step
  smchem <- SmcChem$new(smis = smis0, v_qsprpred = qsprpred_env, v_engram = engram_poly, v_temp = c(temp[i],temp[i]))
  smchem$smcexec(niter = 20, nsteps = max(10,round(para$ngram.O/2)), preorder = 0.2)
  #/get top scored smiles with LC screening
  tmp = smchem$get_hiscores(nsmi = para$n.smis, exsim = 0.7)
  fc <- factor(tmp[,1]) #/make sure only checking unique molecules
  tmp.score <- tapply(tmp[,2], fc, median, na.rm = T)
  #/LC screening
  mols <- parse.smiles(names(tmp.score), kekulise=T) #/was temporarily false before
  mols_op <- tryCatch(lapply(mols,do.aromaticity), error = function(e) NULL)
  mols_op <- tryCatch(lapply(mols,do.typing), error = function(e) NULL)
  mols_op <- tryCatch(lapply(mols,do.isotopes), error = function(e) NULL)
  rm(mols_op)
  LC.uni <- matrix(NA,nrow = length(tmp.score),ncol = length(frags))
  for(k in 1:length(frags)){
    tmp <- which(!sapply(mols,is.na))
    LC.uni[tmp,k] <- matches(frags[[k]],mols[tmp])
  }
  tmp.ind <- which(apply(LC.uni,1,any))
  tmp.smis <- names(tmp.score)
  tmp.score <- as.double(tmp.score)
  tmp.score[tmp.ind] <- tmp.score[tmp.ind] + 10 #/a large constant for sorting
    #/record results and pick the top ones
  tmp <- sort.int(tmp.score,decreasing = T,index.return = T)
  gensmis[[i]] <- list(smis = tmp.smis[tmp$ix], score = tmp.score[tmp$ix])
  if(length(gensmis[[i]]$smis) < para$n.sample){
    smis0 <- sample(gensmis[[i]]$smis, para$n.sample, replace = T)
  }else{
    smis0 <- gensmis[[i]]$smis[1:para$n.sample]
  }
  #/predict properties
  iqspr.pred[[i]] <- tryCatch(qsprpred_env$qspr_predict(gensmis[[i]]$smis), error = function(e) NULL)
  
  cat("\rget iqspr: Finished ", i ,"/",para$n.iter,"     ")   
}

### 10. ordering based on iqspr-Tg values first
smis.all <- unlist(lapply(gensmis[(burnin+1):para$n.iter],function(x) x$smis))
Tg.all <- unlist(lapply(iqspr.pred[(burnin+1):para$n.iter],function(x) x[[1]][1,]))
Tm.all <- unlist(lapply(iqspr.pred[(burnin+1):para$n.iter],function(x) x[[1]][2,]))
Tg.std.all <- unlist(lapply(iqspr.pred[(burnin+1):para$n.iter],function(x) sqrt(x[[2]][1,])))
Tm.std.all <- unlist(lapply(iqspr.pred[(burnin+1):para$n.iter],function(x) sqrt(x[[2]][2,])))
fc <- factor(smis.all)
Tg.uni <- tapply(Tg.all, fc, median, na.rm = TRUE)
smis.uni <- names(Tg.uni)
Tm.uni <- tapply(Tm.all, fc, median, na.rm = TRUE)
Tg.std.uni <- tapply(Tg.std.all, fc, median, na.rm = TRUE)
Tm.std.uni <- tapply(Tm.std.all, fc, median, na.rm = TRUE)

### 11. Post-screening with LCP fragments (separate the calculation for storage concern)
LC.uni <- matrix(NA,nrow = length(smis.uni),ncol = length(frags))
tmp.N <- length(smis.uni)
N.del <- 2000
c0 <- 1
for(j in 1:ceiling(tmp.N/N.del)){
  tmp.ind <- c0:min(c(c0+N.del-1,tmp.N))
  
  mols <- parse.smiles(smis.uni[tmp.ind], kekulise=T) #/was temporarily false before
  mols_op <- tryCatch(lapply(mols,do.aromaticity), error = function(e) NULL)
  mols_op <- tryCatch(lapply(mols,do.typing), error = function(e) NULL)
  mols_op <- tryCatch(lapply(mols,do.isotopes), error = function(e) NULL)
  rm(mols_op)
  
  for(k in 1:length(frags)){
    tmp <- which(!sapply(mols,is.na))
    LC.uni[tmp.ind[tmp],k] <- matches(frags[[k]],mols[tmp])
  }
  
  #/increase counter
  c0 <- c0 + N.del
}

### 12. Summarizing candidates after LCP post-screening
tmp <- apply(LC.uni,1,any)
smis.LC <- smis.uni[tmp]
N.smis <- length(smis.LC)
Tg.LC <- Tg.uni[tmp]
Tg.std.LC <- Tg.std.uni[tmp]
Tm.LC <- Tm.uni[tmp]
Tm.std.LC <- Tm.std.uni[tmp]
