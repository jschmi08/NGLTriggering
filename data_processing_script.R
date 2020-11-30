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

#moss et al dataset

datan <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/mosdAT.csv")
names(datan)[1] <- "liq"
datan$liq <- factor(datan$liq, ordered = TRUE, levels = c("Yes", "No"))
datan$event <- as.factor(datan$event)
summary(datan)

nz <- read.csv("/Users/JSchmidt/OneDrive - Yeh & Associates/Documents/Personal Library/Paper/nzevents.csv")
nz$Liq <- factor(nz$Liq, ordered = TRUE, levels = c("Yes", "No"))
summary(nz)
