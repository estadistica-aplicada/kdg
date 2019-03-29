library(devtools)
library(INLA)
library(faraway)
library(gridExtra)
library(brinla)
library(nlme)
library(geoR)
library(geoRglm)
library(sp)
library(fields)

setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")

kdg <-read.table(file = "kdg_selected_pff_cen_inla.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")
#alpha=2
#maxedge=200000

##defining mesh and smothness parameter
param <- expand.grid(alpha=c(0.5, 1, 1.5, 2), maxedge=c(100000,10000,2000))
x <- param[1,]

mapinlaoptim <-apply(param, 1, function (x){ 
  
  loc.obs <- cbind(kdg$Xt, kdg$Yt)
  mesh <- inla.mesh.2d(loc.obs, cutoff = 200,
                       max.edge = x["maxedge"] )
  node <- mesh$idx$loc
  spde<- inla.spde2.matern(mesh=mesh, alpha=0.5)#x["alpha"])
  formula <-LN_KdGlifo ~1+ f(node, model = spde, diagonal = 1e-6)
  inla_spde <- inla(formula, family = 'gaussian', data = kdg
                    ,control.predictor = list(link = 1, compute = TRUE)
                    #,control.compute = list(cpo=TRUE, dic=TRUE, waic=TRUE)
  )
  
  list(x,"hyperpar"=inla_spde$summary.hyperpar)
  
}
)
####mesh
spdeparam <- read.table(file = "clipboard", sep = "\t", header=TRUE)

a <- as.data.frame(do.call(rbind, mapinlaoptim))
b <- do.call(rbind, a[,2])
names<- do.call(rbind, a[,1])
write.table(b, file="b.txt", sep = "\t", row.names = T)

spdeparam=read.table(file="spdeparam.txt", header = T)
### range y sigma2 by Bangliardo Cammeletti et al pag 207
range <- sqrt(8/exp(spdeparam$Theta1))
sigma2 <- 1/(4*pi*((exp(spdeparam$Theta1))^2)*((exp(spdeparam$Theta2))^2))
rsv <- range/sigma2

spdeparam=cbind(spdeparam,range,sigma2,rsv)

##selected parameters
meshparam <- spdeparam[which.max(spdeparam$rsv),]


##mapping
setwd("C:/Users/franc/Dropbox/Franca/Doctorado/Kd/Glifo")

kdg1 <-read.table(file = "kdg_selected_pff_cen_inla.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")

grid2 <-read.table(file = "grid2.txt", header=TRUE, dec = "." ,sep= "\t",  na.strings = ".")#resto del muestreo
grid2 <- cbind(grid2[, c("Xt", "Yt")],scale(grid2[,c("pH","ARENA")],center = T,scale = F),rep(NA,nrow(grid2)), 
               scale((grid2[,c("pH","ARENA")]^2), center = T,scale = F),rep(NA,nrow(grid2)))

colnames(grid2) <- c("Xt","Yt","CENc_pH","CENc_ARENA","CENc_.Al", "CENc_POT_pH", "CENc_POT_ARENA", "CENc_POT_.Al")


grid1 <- cbind(as.data.frame(read.table(file = "mapa_gab.csv", header=TRUE, sep= ",",  na.strings = ".")))
grid1 <- grid1[,c("X","Y")]
colnames(grid1) <- c("Xt","Yt")
grid <- rbind(grid1[,c("Xt","Yt")])
plot(grid, main="prediction points", cex=0.3,asp=1)

####

loc.obs <- cbind(kdg1$Xt, kdg1$Yt)
A1 <- inla.spde.make.A(mesh, loc=loc.obs)

stk <- inla.stack(data=list(resp=kdg1$LN_KdGlifo),A=list(A1,1), 
                  effects=list(i=1:spde$n.spde,
                          data.frame(Intercept=1,
                                    CENc_pH=kdg1$CENc_pH,
                                    CENc_ARENA= kdg1$CENc_ARENA,
                                    CENc_.Al=kdg1$CENc_.Al, 
                                    CENc_POT_pH=kdg1$CENc_POT_pH, 
                                    CENc_POT_ARENA=kdg1$CENc_POT_ARENA, 
                                    CENc_POT_.Al=kdg1$CENc_POT_.Al)),
                                    tag='est')

res <- inla(resp ~ 0 + Intercept+CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH + f(i, model=spde),
             data=inla.stack.data(stk),
             control.predictor=list(A=inla.stack.A(stk)))

project <- inla.mesh.projector(mesh, loc = cbind(grid$Xt, grid$Yt)) 

spa.mean <- inla.mesh.project(project, res$summary.ran$i$mean) 

##building stacks
stkgrid <- inla.stack(data=list(resp=NA), A=list(project$proj$A,1),
                      effects=list(i=1:spde$n.spde,
                                   data.frame(Intercept=1,
                                              CENc_pH=rep(NA,nrow(grid)),
                                              CENc_ARENA= rep(NA,nrow(grid)),
                                              CENc_.Al=rep(NA,nrow(grid)), 
                                              CENc_POT_pH=rep(NA,nrow(grid)), 
                                              CENc_POT_ARENA=rep(NA,nrow(grid)), 
                                              CENc_POT_.Al=rep(NA,nrow(grid)))),
                      tag='prd.grd')
stk.all <- inla.stack(stk, stkgrid)

##res2 with parameters which maxmize rsv

res3 <- inla(resp ~ 0 + Intercept+CENc_.Al+CENc_POT_.Al+CENc_POT_ARENA+CENc_pH+CENc_POT_pH + f(i, model=spde),
              data=inla.stack.data(stk.all),
              control.predictor=list(A=inla.stack.A(stk.all),
                                     compute=TRUE), quantiles=NULL,
              control.results=list(return.marginals.random=FALSE,
                                   return.marginals.predictor=FALSE))

igr <- inla.stack.index(stk.all, 'prd.grd')$data

###Ploting
res2_table <-cbind(project$loc, as.data.frame(res2$summary.fitt[igr,1]), as.data.frame(res2$summary.fitt[igr,2]))

colnames(res2_table) <- c("Xt","Yt","LnKdGlifo.pred","sdLnKdGlifo")

quilt.plot(res2_table[,1], res2_table[,2], exp(res2_table[,3]), nx = 150, ny = 150, main="Glyphosate Kd",xlab="Longitude", ylab="Latitude", asp=1)

res3_table <-cbind(project$loc, as.data.frame(exp(res3$summary.fitt[igr,1])), as.data.frame(exp(res3$summary.fitt[igr,2])))

quilt.plot(res2_table[,1], res2_table[,2], exp(res2_table[,4]), nx = 150, ny = 150, main="Ln_Kd Glyphosate",xlab="Longitude", ylab="Latitude", asp=1)


