####**Paralelization to make a leave one out validation**
#source("https://bioconductor.org/biocLite.R")
#biocLite() 
#biocLite("BiocParallel")

library(BiocParallel)
i=1

inlakdpred <-function(i){
  
  library(devtools)
  library(INLA)
  library(faraway)
  library(gridExtra)
  library(brinla)
  library(nlme)
  library(geoR)
  library(geoRglm)
  library(sp)
  #Loading data
  
  setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")
  kdg1 <-read.table(file = "kdg_selected_pff_cen_inla.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")
  kdg <- kdg1
  
  head(kdg)
 
  #structure residual table
  tablepred <- matrix(NA, nrow = nrow(kdg), ncol = 3)
  
  colnames(tablepred)<- c("inlaspde")
  #test-train set
  train     <- kdg[-i, ]
  train1    <- kdg[-i, ]
  test     <- kdg[i, ]
  test1     <- kdg[i, ]
  #defining missing value
  iNA <- as.data.frame(cbind(test[,c(1:12)],"LN_KdGlifo"=NA))
   
  testi <- as.data.frame(rbind(train,iNA))
  
  ###building mesh
  loc.obs <- cbind(testi$Xt, testi$Yt)
  
  mesh2 <- inla.mesh.2d(loc.obs, cutoff = 200,
                        max.edge = 20000)

  node <- mesh1$idx$loc

  #spde model
  spde1085 <- inla.spde2.matern(mesh=mesh1, alpha=1.085)
  try({
    #define formula
    formula_spde1085 <-LN_KdGlifo ~1+ CENc_.Al+CENc_pH+CENc_ARENA+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH + f(node, model = spde1085, diagonal = 1e-6)
    #inla model
    inla_pred_spde1085 <- inla(formula_spde1085, family = 'gaussian', data = kdg
                           ,control.predictor = list(link = 1, compute = TRUE)
                           #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
    )

    tablepred[i, "inlaspdem1"] <- (inla_pred_spde1085$summary.fitted.values$mean[89])
  })

  try({
    #result table
    tabla=tablepred[!apply(tablepred,1, function(X){all(is.na(X))}),]
    return(tabla)
    return(tablepred)
  })
  
}

#predictive inla values
results_inla<-do.call("rbind",bplapply(1:89, inlakdpred, BPPARAM=SnowParam(workers=2, progressbar=TRUE, type="SOCK")))
#Prediction errors RMSPE
apply(results_inla_fran, 2 ,function (x) {(sqrt(mean((x-kdg$LN_KdGlifo)^2)))/ mean(kdg$LN_KdGlifo)*100})

