

#==============================================
# 03 Motif analysis - Supplementary analyses
#    Giling et al.
#==============================================

# This code performs the analyses reported in
# the Supplementary Information


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



# 2. Import results
#---------------------

# Import mean results of full network (i.e. grounded and free-floating motifs)
# and results of complexity scenario (Supplemental Text 2)
full <- read.csv("Giling_et_all_motif_data_full_network.csv") 
comp <- read.csv("Giling_et_all_motif_data_complexity_scenario.csv")

# transformations
full$block <- as.factor(substr(full$plot,1,2))
full$log2.psr <- log2(full$sowndiv)
comp$block <- as.factor(substr(comp$plot,1,2))
comp$log2.psr <- log2(comp$sowndiv)



# 3. Supplemental Note 1: Sensitivity to alpha/beta 
#---------------------------------------------------

# 3.1 - Omnivore generality

ggplot(full, aes(y=omn.gen.mean, x=sowndiv)) + 
  geom_point(alpha=0.25) +
  geom_smooth(method="lm", se=F, col="black") + 
  #geom_errorbar(aes(x=sowndiv, ymin=omn.gen.mean-omn.gen.sd, ymax=omn.gen.mean+omn.gen.sd), width=0) +
  scale_x_continuous(trans="log2", breaks=c(1,2,4,8,16,60), limits = c(1,60)) +
  xlab("Plant species richness") +
  ylab("Omnivore generality") +
  facet_wrap(~ alpha + beta, labeller = "label_both") +
  theme_bw() + 
  theme(#panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))


# 3.2 - Motif normalised z-scores: s1, s2, s4, s5

  ggplot(full, aes(y=s1.z.mean, x=sowndiv)) + 
    geom_point(alpha=0.25) +
    geom_abline(intercept = 0, slope = 0, linetype="dashed") +
    geom_smooth(method="lm", se=F, col="black") + 
    #geom_errorbar(aes(x=sowndiv, ymin=s1.z.mean-s1.z.sd, ymax=s1.z.mean+s1.z.sd), width=0) +
    scale_x_continuous(trans="log2", breaks=c(1,2,4,8,16,60), limits = c(1,60)) +
    xlab("Plant species richness") +
    ylab("s1 normalised z-score") +
    facet_wrap(~ alpha + beta, labeller = "label_both") +
    theme_bw() + 
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black"))

  ggplot(full, aes(y=s2.z.mean, x=sowndiv)) + 
    geom_point(alpha=0.25) +
    geom_abline(intercept = 0, slope = 0, linetype="dashed") +
    geom_smooth(method="lm", se=F, col="black") + 
    #geom_errorbar(aes(x=sowndiv, ymin=s2.z.mean-s2.z.sd, ymax=s2.z.mean+s2.z.sd), width=0) +
    scale_x_continuous(trans="log2", breaks=c(1,2,4,8,16,60), limits = c(1,60)) +
    xlab("Plant species richness") +
    ylab("s2 normalised z-score") +
    facet_wrap(~ alpha + beta, labeller = "label_both") +
    theme_bw() + 
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black"))

  ggplot(full, aes(y=s4.z.mean, x=sowndiv)) + 
    geom_point(alpha=0.25) +
    geom_abline(intercept = 0, slope = 0, linetype="dashed") +
    geom_smooth(method="lm", se=F, col="black") + 
    #geom_errorbar(aes(x=sowndiv, ymin=s4.z.mean-s4.z.sd, ymax=s4.z.mean+s4.z.sd), width=0) +
    scale_x_continuous(trans="log2", breaks=c(1,2,4,8,16,60), limits = c(1,60)) +
    xlab("Plant species richness") +
    ylab("s4 normalised z-score") +
    facet_wrap(~ alpha + beta, labeller = "label_both") +
    theme_bw() + 
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black"))

  ggplot(full, aes(y=s5.z.mean, x=sowndiv)) + 
    geom_point(alpha=0.25) +
    geom_abline(intercept = 0, slope = 0, linetype="dashed") +
    geom_smooth(method="lm", se=F, col="black") + 
    #geom_errorbar(aes(x=sowndiv, ymin=s5.z.mean-s5.z.sd, ymax=s5.z.mean+s5.z.sd), width=0) +
    scale_x_continuous(trans="log2", breaks=c(1,2,4,8,16,60), limits = c(1,60)) +
    xlab("Plant species richness") +
    ylab("s5 normalised z-score") +
    facet_wrap(~ alpha + beta, labeller = "label_both") +
    theme_bw() + 
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black"))

  
  
# 4. Plot general properties (Fig S1)
# ------------------------------------

# Select scenario
full <- filter(full, alpha==10000 & beta==10000)

# 4.1 Proportions of trophic groups in consumers
full = mutate(full, n.cons = n.herb.mean + n.det.mean + n.det.mean + n.pred.mean)
full$Detritivores = full$n.det.mean / full$n.cons
full$Herbivores = full$n.herb.mean / full$n.cons
full$Omnivores = full$n.omn.mean / full$n.cons
full$Predators = full$n.pred.mean / full$n.cons

# setup for predictions
divseq <- seq(min(full$log2.psr), max(full$log2.psr), length.out = 100)
newdata <- data.frame(log2.psr = divseq)
results <- data.frame(name=NA, NumDF=NA, DenDF=NA, F=NA, beta.mean=NA, beta.sd=NA, beta.p=NA, int.mean=NA, int.sd=NA, int.p=NA)

# function for plotting and analysis
plot.prop <- function(ytext, yaxt=TRUE, xaxt=TRUE) {

  full$y <- full[,ytext]
  
  plot(y ~ log2.psr, type='n', data=full, ylim=c(0,0.8), xlim=(c(0,6)),
       xaxt='n', yaxt='n', las=2, ylab = NA, xlab = NA, main=ytext)
  if(xaxt==TRUE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=c(1,2,4,8,16,60))
  if(xaxt==FALSE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=NA)
  if(yaxt==TRUE) axis(2, at=seq(0,0.8,0.2), lab=seq(0,0.8,0.2), las=1)
  if(yaxt==FALSE) axis(2, at=seq(0,0.8,0.2), lab=NA, las=1)
  
  points(y ~ log2.psr, data=full, pch=16, col=add.alpha("black",0.1))
  lmm <- NULL
  lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=full)
    # nb: this is not a binomial model because the mean numbers are not integers
  
  #plot(lmm)
  lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
  results <- rbind(results, c(ytext, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
  newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
  points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col="black")
  
  return(results)
}

# run analysis and print figure
par(mfrow=c(1,4), oma=c(3,3,1,0), mar=c(1,1,2,1))
results <- plot.prop("Detritivores", TRUE, TRUE)
results <- plot.prop("Herbivores", FALSE, TRUE)
results <- plot.prop("Omnivores", FALSE, TRUE)
results <- plot.prop("Predators", FALSE, TRUE)



# 5. Supplemental Text 2: Habitat complexity scenario 
#-------------------------------------------------------- 
  
# Same as Fig 3 but for complexity scenario results
  
comp <- filter(comp, alpha==10000 & beta==10000)

# setup
col.list <- c(rgb(0,158,115,max=255), rgb(0,114,178,max=255), rgb(230,159,0,max=255),  rgb(204,121,167,max=255), "grey50")

# dataframe for predicting
divseq <- seq(min(full$log2.psr), max(full$log2.psr), length.out = 100)
newdata <- data.frame(log2.psr = divseq)
results <- data.frame(name=NA, NumDF=NA, DenDF=NA, F=NA, beta.mean=NA, beta.sd=NA, beta.p=NA, int.mean=NA, int.sd=NA, int.p=NA)
  
# function for plotting and analysis
plot.zscores <- function(y, yaxt=TRUE, xaxt=TRUE, col=1) {
  
  comp$y <- comp[,y]

  plot(y ~ log2.psr, type='n', data=comp, ylim=c(-0.8,0.8), xlim=(c(0,6)),
       xaxt='n', yaxt='n', las=2, ylab = NA, xlab = NA)
  abline(h=0, lty=2)
  if(xaxt==TRUE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=c(1,2,4,8,16,60), cex.axis=0.8)
  if(xaxt==FALSE)axis(1, at=log2(c(1,2,4,8,16,60)), lab=NA)
  if(yaxt==TRUE) axis(2, at=seq(-0.8,0.8,0.4), las=1, cex.axis=0.8)
  if(yaxt==FALSE) axis(2, at=seq(-0.8,0.8,0.4), lab=NA, las=1)
  
  points(y ~ log2.psr, data=comp, pch=16, col=add.alpha(col.list[col],0.15))
  lmm <- NULL
  lmm <-lmer(y ~ log2.psr + (1|block/plot) + (1|time.period), data=comp)
  #plot(lmm)
  lmma <- anova(lmm); lmmc <- summary(lmm)$coefficients
  results <- rbind(results, c(y, lmma$NumDF, lmma$DenDF, lmma$`F value`, lmmc[2,1], lmmc[2,2], lmmc[2,5], lmmc[1,1], lmmc[1,2], lmmc[1,5]))
  newdata$ypred <- predict(lmm, newdata, re.form=NA) # unconditional (level 0) - for re see https://www.rdocumentation.org/packages/lme4/versions/1.1-15/topics/predict.merMod 
  points(newdata$ypred ~ newdata$log2.psr, typ='l', lwd=2, col=col.list[col])
  
  return(results)
}

# run analysis and print figure
par(mfrow=c(2,2), oma=c(2,2,1,0), mar=c(1,1,2,1))
results <- plot.zscores("s1.z.mean", TRUE, FALSE, 1)
results <- plot.zscores("s2.z.mean", FALSE, FALSE, 2)
results <- plot.zscores("s4.z.mean", TRUE, TRUE, 3)
results <- plot.zscores("s5.z.mean", FALSE, TRUE, 4)


  









