
setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")

kdg=read.table(file = "kdg_selected_pff_cen.txt", header=TRUE, sep= "\t", dec = ".") #, na.strings = ".")
head(kdg)
library(caret)
library(gbm)
library(dismo)

#Option1. GBR (Elith et al., 2008) 
##firs step of variables selection, grid of possible parameters, gbm
tree.complexity=c(1,2,3,4,5,6,7,8,9,10)
learning.rate=c(0.001, 0.005, 0.01)
bag.fraction=c(0.3,0.4,0.5,0.6)
param=expand.grid(tree.complexity, learning.rate, bag.fraction)
colnames(param) <- c("tree.complexity", "learning.rate", "bag.fraction")


optimbrt2 <- apply(param,1, function(x){
                        mod <- gbm.step(data=kdg,
                        gbm.x = 4:21,
                        gbm.y = 22,
                        family = "gaussian",
                        tree.complexity = x[1],
                        learning.rate = x[2],
                        bag.fraction = x[3])
                        
                        MSE=sqrt(mean((mod$fit-kdg$LN_KdGlifo)^2))
                        
                        list("Parm"=c(x,"MSE"=MSE, "deviance.mean"=mod$cv.statistics$deviance.mean, 
                                      "corr"=mod$self.statistics$correlation))#, "Modelo"=mod)
                        
})

optim2 <- do.call(rbind,lapply(optimbrt2,"[[",1))

dim(param)
dim(optim2)
colnames(optim2) <- c("tc", "lr","bf","MSE","deviance","corr")

write.table(optim2, file="optim.txt")

which.min(optim2[,"MSE"])

optimbrt2[82]

which.min(optim2[,"deviance"])

optimbrt2[35]

which.max(optim2[,"corr"])

optimbrt2[35]

do.call(data.frame, optimbrt2) 

optim2[35,4]/mean(kdg$LN_KdGlifo)*100
##optim model
mod_optim <- gbm.step(data=kdg,
                gbm.x = 4:21,
                gbm.y = 22,
                family = "gaussian",
                tree.complexity = optim2[which.min(optim2[,"deviance"]),1],
                learning.rate = optim2[which.min(optim2[,"deviance"]),2],
                bag.fraction = optim2[which.min(optim2[,"deviance"]),3])

sqrt(mean((mod_optim$residuals)^2))/mean(kdg$LN_KdGlifo)*100
mod_optim$contributions

class(kdg$Kppm)

mod_optim_simplify <- gbm.simplify(mod_optim, n.folds=10 ,n.drops = 20)
mod_optim_simplify$final.drops
mod_optim_simplify$pred.list$preds.10 
mod_optim_simplify$final.drops

kdg1 <-kdg[,c("X.Al","ARC","X.Fe","Mn","pH","ARENA","CC","CO","LN_KdGlifo")]

class(kdg1)

####Option 2 selection variables with GBM with caret
gbmGrid <-  expand.grid(interaction.depth = c(1:6), 
                        n.trees = (10:30)*50, 
                        shrinkage = c(0.1,0.01, 0.001),
                        n.minobsinnode = c(5,4,3,2))
control<- trainControl(method="repeatedcv", number=5,repeats =5)

gbm_optim_caret<- train(x=kdg1[,4:14],y=kdg1[,15], 
              method = "gbm",
              trControl = control,
              verbose = FALSE,
              #importance=T,
              metric="RMSE",
              tuneGrid = gbmGrid)


# Linear model 2 step selection covariates over_fit
kdg_ov_fit <- cbind(kdg,"ARC2"=kdg$ARC^2,"pH2"=kdg$pH^2,"CO2"=kdg$CO^2,"X.Fe2"=kdg$X.Fe^2,"CC2" =kdg$CC^2,
                    "X.Al2" =kdg$X.Al^2, "Mn2"=kdg$Mn^2)
head(kdg_ov_fit)
getwd()
write.table(as.data.frame(kdg_ov_fit), file="kdg_ov_fit.txt")
library(nlme)

lm_ov_fit<- gls(LN_KdGlifo ~1+ARC+pH+CO+X.Fe+CC+X.Al+Mn+ARC2+pH2+CO2+X.Fe2+CC2+X.Al2+Mn2+Xt+Yt+Xt^2+Yt^2
              ,data= kdg_ov_fit)
summary(lm_ov_fit)

library(automap)
library(gstat)

coordinates(kdg_ov_fit)<-~Xt+Yt

variogram_mod_ov_fit_auto <- autofitVariogram(lm_ov_fit$residuals~1,kdg_ov_fit, c("Sph", "Exp", "Gau", "Mat", "Ste"),cressie=T)
plot(variogram_mod_ov_fit_auto)
variogram_mod_ov_fit_auto$var_model

krig_res <- krige(residuals(lm_ov_fit)~1, kdg_ov_fit, newdata=kdg_ov_fit, variogram_mod_ov_fit_auto$var_model,nmax=25)

lm_ov_fit_pred <- lm_ov_fit$fitted + krig_res$var1.pred

plot(lm_ov_fit_pred,kdg$LN_KdGlifo)

# All Subsets covariates in Regression
library(leaps)
leaps<-regsubsets(LN_KdGlifo ~1+ARC+pH+CO+X.Fe+CC+X.Al+Mn+ARC2+pH2+CO2+X.Fe2+CC2+X.Al2+Mn2, 
                  method=c("exhaustive"), data=kdg_ov_fit)
# view results 
summary(leaps)
# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
