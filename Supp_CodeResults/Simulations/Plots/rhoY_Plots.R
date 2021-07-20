setwd("/Volumes/blue/amolstad/amolstad/MTeQTL/Simulations/Results")
ours.R2 <- NULL
Hu.R2 <- NULL
lasso.R2 <- NULL
GS.lasso.R2 <- NULL
GS.Hu.R2 <- NULL
k20.Hu.R2 <- NULL
ours.TPR <- NULL
Hu.TPR <- NULL
lasso.TPR <- NULL
GS.lasso.TPR <- NULL
GS.Hu.TPR <- NULL
k20.Hu.TPR <- NULL
ours.MS <- NULL
Hu.MS <- NULL
lasso.MS <- NULL
GS.lasso.MS <- NULL
GS.Hu.MS <- NULL
k20.Hu.MS <- NULL
ours.timing <- NULL
k <- 1
rhoY <- NULL

for (ll in c(0.0, 0.1, 0.3, 0.5, 0.7)) {

  for (j in 1:500) {

    temp <- readRDS(paste("rhoYY_",100*ll,"_R2_10_eQTLs_15_", j, ".RDS", sep=""))

    ours.R2 <- c(ours.R2, mean(temp$ours.R2.test))
    Hu.R2 <- c(Hu.R2, mean(temp$Hu.R2.test))  
    lasso.R2 <- c(lasso.R2, mean(temp$Lasso.R2.test))  
    GS.lasso.R2 <- c(GS.lasso.R2, mean(temp$GS.Lasso.R2.test))
    GS.Hu.R2 <- c(GS.Hu.R2, mean(temp$GS.Hu.R2.test))
    k20.Hu.R2 <- c(k20.Hu.R2, mean(temp$k20.Hu.R2.test))

    ours.TPR <- c(ours.TPR, temp$ours.TPR)
    Hu.TPR <- c(Hu.TPR, temp$Hu.TPR)
    lasso.TPR <- c(lasso.TPR, temp$Lasso.TPR)
    GS.lasso.TPR <- c(GS.lasso.TPR, temp$GS.Lasso.TPR)
    GS.Hu.TPR <- c(GS.Hu.TPR, temp$GS.Hu.TPR)
    k20.Hu.TPR <- c(k20.Hu.TPR, temp$k20.Hu.TPR)
    
    ours.MS <- c(ours.MS, temp$ours.MS)
    Hu.MS <- c(Hu.MS, temp$Hu.MS)
    lasso.MS <- c(lasso.MS, temp$Lasso.MS)
    GS.lasso.MS <- c(GS.lasso.MS, temp$GS.Lasso.MS)
    GS.Hu.MS <- c(GS.Hu.MS, temp$GS.Hu.MS)
    k20.Hu.MS <- c(k20.Hu.MS, temp$k20.Hu.MS)
    ours.timing <- c(ours.timing, temp$ours.time[3])

    rhoY[k] <- ll
    k <- k + 1

  }
}



library(ggplot2)
library(cowplot)
library(viridis)
Method <- rep(c("Cov-MT", "MT", "EN", "GS-EN", "GS-MT"
  "KNN(20)-MT"), each=length(rhoY))
Method <- factor(Method, levels=c("Cov-MT", "MT", "EN", "GS-EN", "GS-MT",
                                  "KNN(20)-MT"), ordered=TRUE)

dat <- data.frame("Method" = Method, 
  "R2" = c(ours.R2, Hu.R2, lasso.R2, GS.lasso.R2, GS.Hu.R2,
    k20.Hu.R2), 
  "rhoY" = rep(rhoY, 6))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}



df2 <- data_summary(dat, varname="R2", 
                    groupnames=c("Method", "rhoY"))
p1 <- ggplot(df2, aes(x=as.factor(rhoY), y=R2, group=Method, color=Method, linetype=Method)) +
  ylab(expression(paste("Average test set ", R^2, sep=""))) + 
  geom_line(position=position_dodge(0.2), size=.75) + xlab(expression(rhoY)) + 
  geom_point(position=position_dodge(0.2), show.legend=FALSE)+
  geom_errorbar(aes(ymin=R2-2*se, ymax=R2+2*se), width=0, linetype="solid", 
                 position=position_dodge(0.2), show.legend=FALSE) + 
   theme_classic() + scale_color_viridis(discrete=TRUE, option="A", begin=0, end=.8)


dat <- data.frame("Method" = Method, 
  "TPR" = c(ours.TPR, Hu.TPR, lasso.TPR, GS.lasso.TPR, GS.Hu.TPR, 
   # k5.Hu.TPR, k10.Hu.TPR, 
    k20.Hu.TPR), 
  "rhoY" = rep(rhoY, 6))

df2 <- data_summary(dat, varname="TPR", 
                    groupnames=c("Method", "rhoY"))

p2 <- ggplot(df2, aes(x=as.factor(rhoY), y=TPR, group=Method, color=Method, linetype=Method)) + 
  geom_line(position=position_dodge(0.2), size=.75) + ylab("LD-adjusted true positive rate") + 
  geom_point(position=position_dodge(0.2), show.legend=FALSE)+ xlab(expression(rhoY)) + 
  geom_errorbar(aes(ymin=TPR-2*se, ymax=TPR+2*se), width=0, linetype="solid", 
                 position=position_dodge(0.2), show.legend=FALSE) + 
   theme_classic() + scale_color_viridis(discrete=TRUE, option="A", begin=0, end=.8)


dat <- data.frame("Method" = Method, 
  "FNR" = c(ours.MS, Hu.MS, lasso.MS, GS.lasso.MS, GS.Hu.MS, 
    k20.Hu.MS)/34162, 
  "rhoY" = rep(rhoY, 6))

df2 <- data_summary(dat, varname="FNR", 
                    groupnames=c("Method", "rhoY"))

p3 <- ggplot(df2, aes(x=as.factor(rhoY), y=FNR, group=Method, color=Method, linetype=Method)) + 
  geom_line(position=position_dodge(0.2), size=.75) + ylab("Model size") + 
  geom_point(position=position_dodge(0.2), show.legend=FALSE) + xlab(expression(rhoY)) + 
  geom_errorbar(aes(ymin=FNR-2*se, ymax=FNR+2*se), width=0, linetype="solid", 
                 position=position_dodge(0.2), show.legend=FALSE) + 
   theme_classic() + scale_color_viridis(discrete=TRUE, option="A", begin=0, end=.8)


pdf("/Volumes/blue/amolstad/amolstad/MTeQTL/Simulations/rhoYY_Simulations.pdf", width=11, height=3.2)
plot_grid(p1+ theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), rel_widths = c(1,1, 1), ncol=3)
dev.off()



dat <- data.frame( 
  "Timing" = ours.timing,
  "rhoY" = rhoY)

mat1 <- matrix(unlist(tapply(dat$Timing/60, as.factor(dat$rhoY), summary)), nrow=5, byrow=TRUE)

