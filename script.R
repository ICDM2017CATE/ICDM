# Setting No.1 in Susan's PNAS paper
library(ggplot2)
library(ggthemes)
library(reshape2)
theme_Publication <- function(base_size=25, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.6, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))

}
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",
                 manual_pal(values = c("#1664d9","red","orange","green4","purple2","#1f9eb3",
                                       "#d93572","#f781bf","#e41a1c")), ...)

}
estimateVariance <- function(model, train, test){
  # model <- prune_CT
  # train <- dataTrain[[1]]
  # test <- dataTest[[1]]
  train$leaves <- predict(model, newdata=train, type = 'vector')
  test$leaves <- predict(model, newdata=test, type = 'vector')

  train$leavesf <- factor(round(train$leaves,4))
  test$leavesf <- factor(round(test$leaves,4))

  if (length(levels(train$leavesf)) == 1){ #length(levels(xxx)) gets the number of leaves in the tree

    modelTrain <- lm(y~w, data=train)
    modelTest <- lm(y~w, data=train)
    return (c(0,0))
  } else{

    modelTrain <- lm(y ~ leavesf + leavesf * w - w - 1, data=train) # the formula does not contain w and intercept
    modelTest <- lm(y ~ leavesf + leavesf * w - w - 1, data=test)

    # extract the coefficient vectors which are the leaf treatment effects
    coefnumh <- length(coef(modelTest))
    coefnuml <- length(coef(modelTest))/2 + 1

    Train.coeftr <- coef(modelTrain)[coefnuml:coefnumh]
    Test.coeftr <- coef(modelTest)[coefnuml:coefnumh]

    # calculate leaf probabilities

    leafprobTrain <- tapply(train$y,list(train$leavesf),length)
    leafprobTest <- tapply(test$y,list(test$leavesf),length)
    leafprob <- (leafprobTrain + leafprobTest) / (nrow(train) + nrow(test))

    # leafprob <- (leafprobEst + leafprobTrain + leafprobTest)/(nrow(dataEst) + nrow(dataTrain) + nrow(dataTest))

    #calculate variance of estimated treatment effects--typically this is higher in the training set, since there is overfitting there
    Train.coefvar <- sum(leafprob * Train.coeftr^2)-(sum(leafprob*Train.coeftr)^2)
    # Train.coefvar <- sum(leafprob * levels(dataTrain[[i]]$leavesf)^2) - (sum(leafprob * dataTrain[[i]]$leaves)^2)
    Test.coefvar <- sum(leafprob * Test.coeftr^2)-(sum(leafprob*Test.coeftr)^2)

    # print("Variance of estimated treatment effects: Train, Test Sets")
    # print("Typically train has higher var--more extreme estimates--due to overfitting")
    return(c(Train.coefvar, Test.coefvar))
  }
}

rep <- 1; # number of experiment repititions
max_depth <- 15
dataList <- list(); dataTrain <- list(); dataTest <- list(); formula <- list()
# generate data
n <- 500000
min_size <- 500
for (i in 1:rep){
  # p is number of variables,
  # outcome is whether outcome is affected by covariates (except treatment effect)
  # pt is the number of variables that affects treatment effect
  # py is the number of variables that affects outcome but not treatment effect (if outcome = T)
  # dataList[[i]] <- susan1(n, p = 2, treatsize = .5)
  # dataList[[i]] <- susan11(n, p = 20, pt = 4, py = 4, treatsize = 1)
  # dataList[[i]] <- susan2(n, p = 20, outcome = T, pt = 4, py = 4, treatsize = .5)
  # dataList[[i]] <- constant(n, p = 20, outcome = T, pt = 4, py = 4, treatsize = .5)
  dataList[[i]] <- susan3(n, p = 50, outcome = T, pt = 4, py = 4, treatsize = 1)
  # dataList[[i]] <- nonlinear1(n, p = 20, outcome = T, pt = 4, py = 4, treatsize = .5)
  data <- dataList[[i]][[1]]
  formula[[i]] <- dataList[[i]][[2]]
  dataTrain[[i]] <- data[1:(n/2),]
  dataTest[[i]] <- data[(n/2+1):n,]
}
# construct model from train and predict on test
CTList<-list(); MHList<-list(); MOList<-list(); TSTATSList <- list(); TOTList <- list(); FITList<- list();
pb <- txtProgressBar(min = 0, max = rep, initial = 0, style = 3)
for (i in 1:rep){
  setTxtProgressBar(pb, i)
  t<-proc.time()
  treeCT<- causalTree(as.formula(paste("y~",formula[[i]])),
                      data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                      split.Rule="CT", split.Honest=F, cv.option="CT", minsize = min_size,
                      cv.alpha = 1, xval=0, cp=0,
                      propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  t<-proc.time()
  treeMH<- causalTree(as.formula(paste("y~",formula[[i]])),
                      data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                      split.Rule="user2", split.Honest=F, cv.option="CT", minsize = min_size,
                      cv.alpha = 1, xval=0, cp=0,
                      propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  t=proc.time()
  treeMO<- causalTree(as.formula(paste("y~",formula[[i]])),
                      data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                      split.Rule="user3", split.Honest=F, cv.option="CT", minsize = min_size,
                      cv.alpha = 1, xval=0, cp=0,
                      propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  t = proc.time()
  treeTOT<- causalTree(as.formula(paste("y~",formula[[i]])),
                       data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                       split.Rule="TOT", split.Honest=F, cv.option="CT", minsize = min_size,
                       cv.alpha = 1, xval=0, cp=0,
                       propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  t = proc.time()
  treeTSTATS<- causalTree(as.formula(paste("y~",formula[[i]])),
                          data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                          split.Rule="tstats", split.Honest=F, cv.option="CT", minsize = min_size,
                          cv.alpha = 1, xval=0, cp=0,
                          propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  t = proc.time()
  treeFIT<- causalTree(as.formula(paste("y~",formula[[i]])),
                       data = dataTrain[[i]], treatment = dataTrain[[i]]$w,
                       split.Rule="fit", split.Honest=F, cv.option="CT", minsize = min_size,
                       cv.alpha = 1, xval=0, cp=0,
                       propensity = rep(0.5, length(dataTrain[[i]]$w)))
  proc.time()-t
  # p<-ncol(dataTrain[[i]])-3
  # forest <- causalForest(as.formula(paste("y~",formula[[i]])), data=dataTrain[[i]], treatment=dataTrain[[i]]$w,
  #                    split.Rule="CT", split.Honest=T,  split.Bucket=F, bucketNum = 5,
  #                    bucketMax = 100, cv.option="CT", cv.Honest=T, minsize = 1L,
  #                    split.alpha = 0.5, cv.alpha = 0.5,
  #                    sample.size.total = floor(nrow(dataTrain[[i]])/2), sample.size.train.frac = .5,
  #                    mtry = ceiling(ncol(dataTrain[[i]])/3), nodesize = 1, num.trees = 500,
  #                    ncolx = p, ncov_sample = floor(p/3)
  # )
  # predictForest <- predict(forest, newdata=dataTest[[i]], type="vector")
  # MSEtau_infeasForest[i] <- mean((predictForest - dataTest[[i]]$tau_true)^2, na.rm = F)

  cpCT <- treeCT$cptable
  cpMH <- treeMH$cptable
  cpMO <- treeMO$cptable
  cpTOT <- treeTOT$cptable
  cpTSTATS <- treeTSTATS$cptable
  cpFIT <- treeFIT$cptable
  predictCT <- list(); predictMH <- list(); predictMO <- list(); predictTOT <- list(); predictTSTATS <- list(); predictFIT<-list();
  for (j in 1:max_depth){
    if (nrow(cpCT) >= j){
      prune_CT <- prune(treeCT, cpCT[j,1])
    }
    else{
      prune_CT <- prune(treeCT, cpCT[nrow(cpCT),1])
    }
    if (nrow(cpMH) >= j){
      prune_MH <- prune(treeMH, cpMH[j,1])
    }
    else{
      prune_MH <- prune(treeMH, cpMH[nrow(treeMH),1])
    }
    if (nrow(cpMO) >= j){
      prune_MO <- prune(treeMO, cpMO[j,1])
    }
    else{
      prune_MO <- prune(treeMO, cpMO[nrow(treeMO),1])
    }
    if (nrow(cpTOT) >= j){
      prune_TOT <- prune(treeTOT, cpTOT[j,1])
    }
    else{
      prune_TOT <- prune(treeTOT, cpTOT[nrow(treeTOT),1])
    }
    if (nrow(cpTSTATS) >= j){
      prune_TSTATS <- prune(treeTSTATS, cpTSTATS[j,1])
    }
    else{
      prune_TSTATS <- prune(treeTSTATS, cpTSTATS[nrow(treeTSTATS),1])
    }
    if (nrow(cpFIT) >= j){
      prune_FIT <- prune(treeFIT, cpFIT[j,1])
    }
    else{
      prune_FIT <- prune(treeFIT, cpFIT[nrow(treeFIT),1])
    }
    predictCT[[j]] = predict(prune_CT, newdata = dataTest[[i]], type = "vector")
    predictMH[[j]] = predict(prune_MH, newdata = dataTest[[i]], type = "vector")
    predictMO[[j]] = predict(prune_MO, newdata = dataTest[[i]], type = "vector")
    predictTOT[[j]] = predict(prune_TOT, newdata = dataTest[[i]], type = "vector")
    predictTSTATS[[j]] = predict(prune_TSTATS, newdata = dataTest[[i]], type = "vector")
    predictFIT[[j]] = predict(prune_FIT, newdata = dataTest[[i]], type = "vector")
  }
  CTList[[i]] <- predictCT
  MHList[[i]] <- predictMH
  MOList[[i]] <- predictMO
  TOTList[[i]] <-predictTOT
  TSTATSList[[i]] <-predictTSTATS
  FITList[[i]] <- predictFIT
}

# compute RMSE and wRMSE
w_MSEtau_infeasCT <- matrix(0, rep, max_depth); w_MSEtau_infeasMH <- matrix(0, rep, max_depth); w_MSEtau_infeasMO <- matrix(0, rep, max_depth); w_MSEtau_infeasTOT <- matrix(0, rep, max_depth); w_MSEtau_infeasTSTATS <- matrix(0, rep, max_depth); w_MSEtau_infeasFIT <- matrix(0, rep, max_depth)
w1 <- .2; w2 <- 1
MSEtau_infeasCT <- matrix(0, rep, max_depth); MSEtau_infeasMH <- matrix(0, rep, max_depth); MSEtau_infeasTOT <- matrix(0, rep, max_depth); MSEtau_infeasTSTATS <- matrix(0, rep, max_depth); MSEtau_infeasFIT <- matrix(0, rep, max_depth); MSEtau_infeasMO <- matrix(0, rep, max_depth); MSEtau_infeasMT <- matrix(0, rep, max_depth)
MAEtau_infeasCT <- matrix(0, rep, max_depth); MAEtau_infeasMH <- matrix(0, rep, max_depth); MAEtau_infeasTOT <- matrix(0, rep, max_depth); MAEtau_infeasTSTATS <- matrix(0, rep, max_depth); MAEtau_infeasFIT <- matrix(0, rep, max_depth); MAEtau_infeasMO <- matrix(0, rep, max_depth); MAEtau_infeasMT <- matrix(0, rep, max_depth)
for (i in 1:rep){
  for (j in 1:max_depth){
    w_MSEtau_infeasCT[i,j] = sqrt(weighted.mean((CTList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                           ifelse(sign(CTList[[i]][[j]] * dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    w_MSEtau_infeasMH[i,j] = sqrt(weighted.mean((MHList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                            ifelse(sign(MHList[[i]][[j]] * dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    w_MSEtau_infeasMO[i,j] = sqrt(weighted.mean((MOList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                           ifelse(sign(MOList[[i]][[j]] * dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    w_MSEtau_infeasTOT[i,j] = sqrt(weighted.mean((TOTList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                            ifelse(sign(TOTList[[i]][[j]]* dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    w_MSEtau_infeasTSTATS[i,j] = sqrt(weighted.mean((TSTATSList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                            ifelse(sign(TSTATSList[[i]][[j]] * dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    w_MSEtau_infeasFIT[i,j] = sqrt(weighted.mean((FITList[[i]][[j]] - dataTest[[i]]$tau_true)^2,
                                            ifelse(sign(FITList[[i]][[j]] * dataTest[[i]]$tau_true)==1,w1,w2), na.rm = F))
    MSEtau_infeasCT[i,j] = sqrt(mean((CTList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))
    MSEtau_infeasMH[i,j] = sqrt(mean((MHList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))
    MSEtau_infeasMO[i,j] = sqrt(mean((MOList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))
    MSEtau_infeasTOT[i,j] = sqrt(mean((TOTList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))
    MSEtau_infeasTSTATS[i,j] = sqrt(mean((TSTATSList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))
    MSEtau_infeasFIT[i,j] = sqrt(mean((FITList[[i]][[j]] - dataTest[[i]]$tau_true)^2, na.rm = F))

    MAEtau_infeasCT[i,j] = sum(abs(CTList[[i]][[j]] - dataTest[[i]]$tau_true))
    MAEtau_infeasMH[i,j] = sum(abs(MHList[[i]][[j]] - dataTest[[i]]$tau_true))
    MAEtau_infeasMO[i,j] = sum(abs(MOList[[i]][[j]] - dataTest[[i]]$tau_true))
    MAEtau_infeasTOT[i,j] = sum(abs(TOTList[[i]][[j]] - dataTest[[i]]$tau_true))
    MAEtau_infeasTSTATS[i,j] = sum(abs(TSTATSList[[i]][[j]] - dataTest[[i]]$tau_true))
    MAEtau_infeasFIT[i,j] = sum(abs(FITList[[i]][[j]] - dataTest[[i]]$tau_true))


  }
}
mean_infeasCT<-colMeans(MSEtau_infeasCT); mean_infeasMH<-colMeans(MSEtau_infeasMH); mean_infeasMO<-colMeans(MSEtau_infeasMO); mean_infeasTOT<-colMeans(MSEtau_infeasTOT); mean_infeasTSTATS<-colMeans(MSEtau_infeasTSTATS); mean_infeasFIT<-colMeans(MSEtau_infeasFIT)
w_mean_infeasCT<-colMeans(w_MSEtau_infeasCT); w_mean_infeasMH<-colMeans(w_MSEtau_infeasMH); w_mean_infeasMO<-colMeans(w_MSEtau_infeasMO); w_mean_infeasTOT<-colMeans(w_MSEtau_infeasTOT); w_mean_infeasTSTATS<-colMeans(w_MSEtau_infeasTSTATS); w_mean_infeasFIT<-colMeans(w_MSEtau_infeasFIT)
mean_maeCT<-colMeans(MAEtau_infeasCT); mean_maeMH <- colMeans(MAEtau_infeasMH); mean_maeMO<-colMeans(MAEtau_infeasMO); mean_maeTOT<-colMeans(MAEtau_infeasTOT);mean_maeTSTATS<-colMeans(MAEtau_infeasTSTATS); mean_maeFIT<-colMeans(MAEtau_infeasFIT)
#show MAE
df<- data.frame(1:max_depth, mean_maeFIT, mean_maeCT, mean_maeTOT,
                mean_maeMH, mean_maeTSTATS, mean_maeMO)
colnames(df) <- c("n_split", "RT", "CT", "TOT", "MH", "TSTATS", "MO")
Molten <- melt(df, id.vars = "n_split")
ggplot(Molten, aes(x = n_split, y = value, colour = variable)) +
  labs(colour = "Criterion") + xlab("Number of splits") + ylab("MAE") + geom_point(size=2) +
  geom_line(size=1.2) + labs(title=paste("Setting 1, n = ",n)) + theme_Publication() +scale_colour_Publication()

#show wRMSE
df<- data.frame(1:max_depth, w_mean_infeasFIT, w_mean_infeasCT, w_mean_infeasTOT,
                w_mean_infeasMH, w_mean_infeasTSTATS, w_mean_infeasMO)
colnames(df) <- c("n_split", "RT", "CT", "TOT", "MH", "TSTATS", "MO")
Molten <- melt(df, id.vars = "n_split")
ggplot(Molten, aes(x = n_split, y = value, colour = variable)) +
  labs(colour = "Criterion") + xlab("Number of splits") + ylab("wRMSE") + geom_point(size=2) +
  geom_line(size=1.2, position =  position_dodge(0.5)) + labs(title=paste("Setting 1, n = ",n)) + theme_Publication() +scale_colour_Publication()

#show RMSE
df2<- data.frame(1:max_depth, mean_infeasFIT, mean_infeasCT, mean_infeasTOT,
                mean_infeasMH, mean_infeasTSTATS, mean_infeasMO)
colnames(df2) <- c("n_split", "RT", "CT", "TOT", "MH", "TSTATS", "MO")
Molten2 <- melt(df2, id.vars = "n_split")
# pdf("s4n2000.pdf")
# par(mai=c(0.1,0.1,0.1,0.1))
ggplot(Molten2, aes(x = n_split, y = value, colour = variable)) +
  labs(colour = "Criterion") + xlab("Number of splits") + ylab("RMSE") + geom_point(size=2) +
  geom_line(size=1.2, position = position_dodge(0.5)) + labs(title=paste("Setting 1, n = ",n)) + theme_Publication() +scale_colour_Publication()
# dev.off()


