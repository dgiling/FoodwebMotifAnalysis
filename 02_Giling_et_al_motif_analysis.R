
#======================================
# 02 Analyse motif counts and z-scores
#    Giling et al. 
#======================================

# This code analyses the effect of sown plant 
# species richness on the motif results


# 1. Required libraries and definitions
#---------------------------------------

library(tidyverse)
library(lme4)
library(lmerTest)

add.alpha <- function(col, alpha=1)
{
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

col.list <- c(rgb(0,158,115,max=255), rgb(0,114,178,max=255), rgb(230,159,0,max=255),  rgb(204,121,167,max=255), "grey50")



# 2. Import results and perform calculations
#---------------------------------------------

# Import mean results of full network (i.e. grounded and free-floating motifs)
# and consumer subweb (i.e. free-floating motifs only)
full <- read.csv("Giling_et_all_motif_data_full_network.csv") # full network
sub <- read.csv("Giling_et_all_motif_data_sub_network.csv") # consumer sub web

# Select desired scenario
full <- filter(full, alpha==10000 & beta==10000)
sub <- filter(sub, alpha==10000 & beta==10000)

# Calculations and transformations
full$other.mean = apply(full[,c("s3.mean", paste0("d",1:8,".mean"))],1,sum)
sub$other.mean = apply(sub[,c("s3.mean", paste0("d",1:8,".mean"))],1,sum)

full$log.s1 <- log10(full$s1.mean +1)
full$log.s2 <- log10(full$s2.mean +1)
full$log.s4 <- log10(full$s4.mean +1)
full$log.s5 <- log10(full$s5.mean +1)
full$log.other <- log10(full$other.mean +1)

sub$log.s1 <- log10(sub$s1.mean +1)
sub$log.s2 <- log10(sub$s2.mean +1)
sub$log.s4 <- log10(sub$s4.mean +1)
sub$log.s5 <- log10(sub$s5.mean +1)
sub$log.other <- log10(sub$other.mean +1)

full$block <- as.factor(substr(full$plot,1,2))
full$log2.psr <- log2(full$sowndiv)
sub$block <- as.factor(substr(sub$plot,1,2))
sub$log2.psr <- log2(sub$sowndiv)

# Create empty dataframe for predictions
divseq <- seq(min(full$log2.psr), max(full$log2.psr), length.out = 100)
newdata <- data.frame(log2.psr = divseq)


# 3. Analyse motif frequencies (Fig. 2)
#-----------------------------------------

# setup for plotting and modelling
results <- data.frame(name=NA, NumDF=NA, DenDF=NA, F=NA, beta.mean=NA, beta.sd=NA, beta.p=NA, int.mean=NA, int.sd=NA, int.p=NA)

# define function for plotting and analysis of full web and consumer sub web
plot.freq <- function(y, yaxt=TRUE, xaxt=TRUE, col=1, ylow=log10(50), yhigh=log10(10000), mod=TRUE) {
  
    full$y <- full[,y]
    sub$y <- sub[,y]
   
    plot(y ~ log2.psr, type='n', data=full, ylim=c(ylow,yhigh), xlim=(c(0,6)),
         xaxt='n', yaxt='n', las=2, ylab = NA, xlab = NA)
    if(xaxt==TRUE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=c(1,2,4,8,16,60))
    if(xaxt==FALSE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=NA)
    if(yaxt==TRUE) axis(2, at=log10(c(1,10,100,1000,10000)+1), lab=c(1,10,100,1000,10000), las=1)
    if(yaxt==FALSE) axis(2, at=log10(c(1,10,100,1000,10000)+1), lab=NA, las=1)
    
    # add points for full web
    points(y ~ log2.psr, data=full, pch=16, col=add.alpha(col.list[col],0.15))
    lmm <- NULL
    lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=full)
    #plot(lmm)
    lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
    results <- rbind(results, c(y, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
    newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
    if(mod==TRUE) points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col=col.list[col])
    
    # add points for consumer sub web
    points(y ~ log2.psr, data=sub, pch=21, col=add.alpha(col.list[col],0.15))
    lmm <- NULL
    lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=sub)
    #plot(lmm)
    lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
    results <- rbind(results, c(y, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
    newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
    if(mod==TRUE) points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col=col.list[col], lty=2)

    return(results)
}

# run analysis and print figure
par(mfrow=c(3,2), oma=c(3,3,1,0), mar=c(1,1,2,1))
results <- plot.freq(y="log.s1", TRUE, FALSE, 1)
results <- plot.freq("log.s2", FALSE, FALSE, 2)
results <- plot.freq("log.s4", TRUE, FALSE, 3)
results <- plot.freq("log.s5", FALSE, TRUE, 4)
results <- plot.freq("log.other", TRUE, TRUE, 5, ylow=0,yhigh=log10(500), mod=FALSE)


# 4. Analyse normalisedz-scores (Fig. 3)
#-----------------------------------------

# function for plotting and analysis
plot.zscores <- function(y, yaxt=TRUE, xaxt=TRUE, col=1) {
  
  full$y <- full[,y]
  sub$y <- sub[,y]
  
  plot(y ~ log2.psr, type='n', data=full, ylim=c(-0.8,0.8), xlim=(c(0,6)),
       xaxt='n', yaxt='n', las=2, ylab = NA, xlab = NA)
  abline(h=0, lty=2)
  if(xaxt==TRUE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=c(1,2,4,8,16,60), cex.axis=0.8)
  if(xaxt==FALSE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=NA)
  if(yaxt==TRUE) axis(2, at=seq(-0.8,0.8,0.4), las=1, cex.axis=0.8)
  if(yaxt==FALSE) axis(2, at=seq(-0.8,0.8,0.4), lab=NA, las=1)
  
  # add points for full web
  points(y ~ log2.psr, data=full, pch=16, col=add.alpha(col.list[col],0.15), cex=0.9)
  lmm <- NULL
  lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=full)
  #plot(lmm)
  lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
  results <- rbind(results, c(y, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
  newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
  points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col=col.list[col])
  
  # add points for consumer sub web
  points(y ~ log2.psr, data=sub, pch=21, col=add.alpha(col.list[col],0.15), cex=0.9)
  lmm <- NULL
  lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=sub)
  #plot(lmm)
  lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
  results <- rbind(results, c(y, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
  newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
  points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col=col.list[col], lty=2)
  
  return(results)
}

# run analysis and print figure
par(mfrow=c(2,2), oma=c(2,2,1,0), mar=c(1,1,2,1))
results <- plot.zscores("s1.z.mean", TRUE, FALSE, 1)
results <- plot.zscores("s2.z.mean", FALSE, FALSE, 2)
results <- plot.zscores("s4.z.mean", TRUE, TRUE, 3)
results <- plot.zscores("s5.z.mean", FALSE, TRUE, 4)






