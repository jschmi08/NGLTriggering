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

make_summary_histogram <- function(data, name){
  ggplot(data = as.data.frame(data), aes(x = data)) + geom_histogram(aes(y=..density..), color = "black", fill = "white") + theme_classic() + geom_density() +
    geom_vline(aes(xintercept = quantile(data, 0.5)), linetype = "dashed") + geom_vline(aes(xintercept = quantile(data, 0.05)), linetype = "dashed") + 
    geom_vline(aes(xintercept = quantile(data, 0.95)), linetype = "dashed") + ggtitle(paste("Mean = ", round(mean(data),2),"SD = ", round(sd(data),3))) +
    xlab(name) + ylab("Density")+
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.1))) +
    theme(legend.position = "none",
          axis.text=element_text(size=10),
          axis.title=element_text(size=12),
          plot.title =element_text(size = 12, hjust = 0.5, face = "bold"))
  
}  

datan <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/mosdAT.csv")
names(datan)[1] <- "liq"
l_qc1 <- 1.6
l_csr <- -0.6

draws <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/fmodel2.2p10.csv") #make sure you have the right stanfit
b0 <- draws$mu_b0
b1 <- draws$mu_b1
b2 <- draws$mu_b2

#visualize results 


#visualize results 
plotA <- make_chart()




nz <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/nzevents.csv")
for (i in 1:length(nz$qc1)){
  nz$p[i] <- mean(1/(1+exp(-(b0+b1*boxcox(nz$qc1,l_qc1)[i]+b2*boxcox(nz$CSR,l_csr)[i]))))
}




#precrec!

myRoc <- evalmod(scores = nz$p, labels = nz$Liq)
myAUC <- auc(evalmod(scores = nz$p, labels = nz$Liq)) 
plotB <- autoplot(myRoc, "ROC") + xlab("False Positive Rate") + ylab("True Positive Rate") + scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        plot.title = element_blank()) +
  annotate('text',x=0.5, y=0.1, size = 3.5, label=paste("AUC = ", round(myAUC$aucs[1],3)))




#intercept histogram
plotD <- make_summary_histogram(b0, "Intercept")

#qc1 slope histogram
plotE <- make_summary_histogram(b1, expression(q["c,1"] ~ "Slope"))

#CSR Slope histogram
plotF <- make_summary_histogram(b2, "CSR Slope")

#viz results beta

grid.newpage()
vp1 <- viewport(x = 0, y = 0, 
                height = 0.7, width = 0.55,
                just = c("left", "bottom"),
                name = "lower left")
pushViewport(vp1)
grid.rect()
print(plotA, newpage = FALSE)
popViewport()

vp2.5 <- viewport(x = 0.85, y = 0.2, 
                  height = 0.3, width = 0.3,
                  just = c("right", "bottom"),
                  name = "lower right")
pushViewport(vp2.5)
print(plotB, newpage = FALSE)
popViewport()

vp3 <- viewport(x = 0.95, y = 0.5, 
                height = 0.25, width = 0.45,
                just = c("right", "bottom"),
                name = "mid right")
pushViewport(vp3)
print(plotD, newpage = FALSE)
popViewport()

vp4 <- viewport(x = 0.95, y = 1, 
                height = 0.25, width = 0.45,
                just = c("right", "top"),
                name = "upper right")
pushViewport(vp4)
print(plotE, newpage = FALSE)
popViewport()


vp5 <- viewport(x = 0.05, y = 1, 
                height = 0.25, width = 0.45,
                just = c("left", "top"),
                name = "upper left")
pushViewport(vp5)
print(plotF, newpage = FALSE)
popViewport()

nz$LiqBinary <- as.numeric(as.factor(nz$Liq))-1

brier <- mean((nz$LiqBinary-nz$p)^2)
log <- mean(-log(1-abs(nz$LiqBinary-nz$p)))
brier
log



#fake ROC curve
moss <- datan

for (i in 1:length(moss$liq)){
  moss$p[i] <- mean(1/(1+exp(-(b0+b1*boxcox(moss$qc1_mean,l_qc1)[i]+b2*boxcox(moss$CSR_mean,l_csr)[i]))))
}


pred <- prediction(moss$p, moss$liq) 
rocc <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(rocc, colorize = F)
abline(a=0, b=1, lty=2, lwd=2, col="black") #45 degree line represents performance of a 'coin flip' model
abline(h = 1, lty = 3, lwd=2, col="black")
abline(v = 0, lty = 3, lwd=2, col="black")
auc.perf <- performance(pred, measure = "auc")
text(0.75,0.4, "AUC = ")
text(0.75,0.3,labels = round(auc.perf@y.values[[1]],3))
