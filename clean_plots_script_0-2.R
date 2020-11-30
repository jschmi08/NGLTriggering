library(ROCR)
library(ggplot2)
require(martrixStats)
boxcox <- function(x,l){
  if (l == 0) {
    return(log(x))
  } else {
    return((x^(l)-1)/l)
  }
}
revboxcox <- function(x,l){
  if (l == 0 ){
    return(exp(x))
  } else {
    return((x*l+1)^(1/l))
  }
}

make_chart <- function(){
  xgrid <- seq(boxcox(x = min(datan$qc1_mean), l = l_qc1),boxcox(x = max(datan$qc1_mean), l = l_qc1),0.5) 
  ygrid <- seq(boxcox(x = min(datan$CSR_mean), l = l_csr),boxcox(x = max(datan$CSR_mean), l = l_csr),0.1) 
  data.fit <- expand.grid(qc1 = xgrid, CSR = ygrid) 
  for (i in 1:length(data.fit$qc1)){
    data.fit$z[i] <- mean(1/(1+exp(-(b0+b1*data.fit$qc1[i]+b2*data.fit$CSR[i]))))
    data.fit$qc1[i] <- revboxcox(data.fit$qc1[i],l_qc1)
    data.fit$CSR[i] <- revboxcox(data.fit$CSR[i], l_csr)
  }
  ggplot() + theme_classic()+
    geom_point(data = datan[datan$liq == "Yes",], aes(x = qc1_mean, y = CSR_mean, shape = "moss_yes"), size = 1) +
    geom_point(data = datan[datan$liq == "No",], aes(x = qc1_mean, y = CSR_mean, shape = "moss_no"), size = 1) +
    geom_point(data = nz[nz$Liq == "Yes",], aes(x = qc1, y = CSR, shape = "nz_yes"), size = 1) +
    geom_point(data = nz[nz$Liq == "No",], aes(x = qc1, y = CSR, shape = "nz_no"), size = 1) +
    scale_shape_manual( values = c(1,19,2,17)) +
    stat_contour(data = data.fit,color = "black", size = .8, aes(x = qc1, y = CSR, z = z), 
                 breaks = c(0.95, 0.8, 0.5, 0.2, 0.05)) +
    xlab(expression(bold(q["c,1"] (MPa)))) + ylab("CSR") +
    theme(legend.position = "none",
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title =element_text(size = 12, hjust = 0.5))
}

make_performance_chart <- function(measure, x.measure, auc.measure, auc.label){
  par(mar=c(3,3,3,1), mgp=c(2,.5,0), tck=-.01)
  rocc <- performance(pred, measure = measure, x.measure = x.measure) 
  plot(rocc, colorize = F,xaxs = "r", xlim = c(0,1),
       yaxs = "r", ylim = c(0,1))
  if (x.measure == "fpr") {
    abline(a=0, b=1, lty=2, lwd=2, col="black") #45 degree line represents performance of a 'coin flip' model
    abline(h = 1, lty = 3, lwd=2, col="black")
    abline(v = 0, lty = 3, lwd=2, col="black")
  }
  if (x.measure == "rec") {
    abline(h = 1, lty = 3, lwd=2, col="black")
    abline(v = 1, lty = 3, lwd=2, col="black")
  }
  
  auc.perf <- performance(pred, measure = auc.measure)
  text(0.75,0.4,paste(auc.label, " = "))
  text(0.75,0.3,labels = round(auc.perf@y.values[[1]],3))
  
}



datan <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/mosdAT.csv")
names(datan)[1] <- "liq"
datan$liq <- as.factor(datan$liq)
l_qc1 <- 1.6
l_csr <- -0.6
datan$CSR_star <- boxcox(datan$CSR_mean, l_csr)
datan$qc1_star <- boxcox(datan$qc1_mean, l_qc1)

m2 <- glm(liq ~ qc1_star + CSR_star, family = "binomial", data = datan)
bees <- coefficients(m2)

b0 <- bees[1]
b1 <- bees[2]
b2 <- bees[3]

#vizualize model results
plotA <- make_chart()
plotA

sig <- sqrt(data.fit$z*(1-data.fit$z))
mean(sig)

nz <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/nzevents.csv")
for (i in 1:length(nz$qc1)){
  nz$p[i] <- 1/(1+exp(-(b0+b1*boxcox(nz$qc1,l_qc1)[i]+b2*boxcox(nz$CSR,l_csr)[i])))
}
pred <- prediction(nz$p, nz$Liq) 

nz$LiqBinary <- as.numeric(as.factor(nz$Liq))-1
brier <- mean((nz$p-nz$LiqBinary)^2)
log <- mean(-log(1-abs(nz$p-nz$LiqBinary)))
brier
log

#ROC curve
make_performance_chart(measure = "tpr", x.measure = "fpr", auc.measure = "auc", auc.label = "AUC")

#pr curve
make_performance_chart(measure = "prec", x.measure = "rec", auc.measure = "aucpr", auc.label = "AUC-PR")


