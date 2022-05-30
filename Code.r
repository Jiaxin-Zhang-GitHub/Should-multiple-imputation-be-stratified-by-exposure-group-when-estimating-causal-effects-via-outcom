#################################################################################
#  DataCleaning.r                                                               # 
#  VAHCS - Data Cleaning                                                        #
#################################################################################

# Exposure:     thcw2to6any  label: thc weekly:any w2-6
# Outcome:      cmdtot7      label: cis:total score (14 items)
# Confounders:  alchirsk2to6any  asb2to6any  par_div  par_inchs    

rm(list=ls())
library(haven)
library(tidyr)
library(magrittr) 

# Import the original dataset from stata
setwd("~/R/Project No.1")
dat <- read_dta("VAHCS.Raw.dta") 

# Revise alchirsk2to6any and cmd2to6any; keep relative variables; keep female records; drop missing values in thcw2to6any, asb2to6any and par_div
dat[dat[, 8] == 0, 8] <- apply(dat[dat[, 8] == 0, c(3:7)], 1, sum)
dat[dat[, 24] == 0, 24] <- apply(dat[dat[, 24] == 0, c(19:23)], 1, sum)
dat <- dat %>% subset.data.frame(., sex == 1) %>% 
  .[,c("age_2", "par_inchs", "par_div", "asb2to6any", "cmd2to6any", "alchirsk2to6any", "thcw2to6any", "cmdtot7")] %>%
  subset.data.frame(., complete.cases(thcw2to6any, asb2to6any, par_div, par_inchs))

# Check variables
summary(dat)

# Standarise the log-transformed age_2
dat$age_2 <- log(dat$age_2) %>% `-`(mean(., na.rm = T)) %>% `/`(sd(., na.rm = T))

# Symmetrical and standarise cmdtot7
dat$zcmd <- log(dat$cmdtot7+1) %>% `-`(mean(., na.rm = T)) %>% `/`(sd(., na.rm = T))

# Factorization
dat <- dat[,-8]
dat[2:7] <- lapply(dat[2:7], as.factor)

# Rename
names(dat)[1:8] <- c("age","pari","pard","asb","cmd","alc","thc","zcmd")

# Save the cleaned dataset
write.csv(dat,"VAHCS.CleanData.csv", row.names = FALSE)
save(dat, file = "VAHCS.CleanData.Rda")


#################################################################################
#  TableOne.r                                                                   # 
#  VAHCS - Table One                                                            #
#################################################################################
rm(list=ls())
library(haven)
library('tableone')
load("VAHCS.CleanData.Rda")

# Table One
factorVars <- c("alc", "asb", "cmd", "pard", "pari")
vars <- c("age", "pari", "pard", "asb", "cmd", "alc", "zcmd")
tableOne <- CreateTableOne(vars = vars, strata = "thc", data = dat, factorVars = factorVars, includeNA = F, test = F)
VAHCS.table1 <- print(tableOne, missing = T, explain = T, printToggle = T)

# Save to a CSV file
write.csv(VAHCS.table1, file = "VAHCS-Table1.csv")


#################################################################################
#  CaseStudy.r                                                                  # 
#  VAHCS - Case Study                                                           #
#################################################################################
rm(list=ls())
library(mvnmle)
library(BaylorEdPsych)
library('MASS')
library(ggplot2)
library("mice")
library(car)
library(boot)
load("VAHCS.CleanData.Rda")
set.seed(1957427)

# Kernel Density Plot
ggplot(dat, aes(age)) + geom_density(kernel = "gaussian") +
  scale_color_nejm() + theme_classic()
ggplot(dat, aes(zcmd)) + geom_density(kernel = "gaussian") +
  scale_color_nejm() + theme_classic()


## Method I: Complete-case Analysis
CCA <- lm(zcmd ~ pari + pard + asb + cmd + alc + thc, data = dat)
CCA.ACE <- as.data.frame(t(round(cbind(summary(CCA)[["coefficients"]],confint(CCA))[7,-3],3))) %>% `names<-`(c("estimate","std.error","p.value","2.5 %","97.5 %"))


# Multicollinearity: vif score less than 1.5, provide validation for imputation
vif(CCA) 
# Calculate the optimal number of multiple imputations
round(sum(complete.cases(dat))/nrow(dat)*100)


## Method II: Whole-cohort Multiple Imputation
# Imputating
MI.NI <- mice(data = dat, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm"), vis = "monotone", print = F)
# Diagnostic checking
summary(mice::complete(MI.NI)) 
# Convergence check
plot(MI.NI)
# Imputated values range check
stripplot(MI.NI)
# Pooling
MI.NI.fit <- with(MI.NI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.NI.ACE <- round(summary(pool(MI.NI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)


## Method III: Exposure-outcome Multiple Imputation
# Manipulate predicted matrix
dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
predmat <- make.predictorMatrix(dat)
predmat[c("thc","zcmd"), "thc.zcmd"] <- 0
# Imputating
MI.EO <- mice(dat, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*zcmd)"), vis = "monotone", predictorMatrix = predmat, print = F)
# Diagnostic check
summary(mice::complete(MI.EO)) 
# Convergence check
plot(MI.EO)
# Imputated values range check
stripplot(MI.EO)
# Pooling
MI.EO.fit <- with(MI.EO, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.EO.ACE <- round(summary(pool(MI.EO.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)


## Method IV: Exposure-confounder Multiple Imputation
load("VAHCS.CleanData.Rda")
# Manipulate predicted matrix
dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
predmat <- make.predictorMatrix(dat)
predmat[c("cmd","thc"), "thc.cmd"] <- 0
# Imputating
MI.EC <- mice(dat, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))"), vis = "monotone", predictorMatrix = predmat, print = F)
# Diagnostic check
summary(mice::complete(MI.EI)) 
# Convergence check
plot(MI.EC)
# Imputated values range check
stripplot(MI.EC)
# Pooling
MI.EC.fit <- with(MI.EC, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.EC.ACE <- round(summary(pool(MI.EC.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)


## Method V: Exposure-confounder interaction and exposure-outcome interaction Multiple Imputation 
load("VAHCS.CleanData.Rda")
# Manipulate predicted matrix
dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
predmat <- make.predictorMatrix(dat)
predmat["cmd", "thc.cmd"] <- 0
predmat["zcmd", "thc.zcmd"] <- 0
# Imputating
MI.EOC <- mice(dat, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))","~I(as.numeric(as.character(thc))*zcmd)"), vis = "monotone", predictorMatrix = predmat, print = F)
# Diagnostic check
summary(mice::complete(MI.EOC)) 
# Convergence check
plot(MI.EOC)
# Imputated values range check
stripplot(MI.EOC)
# Pooling
MI.EOC.fit <- with(MI.EOC, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.EOC.ACE <- round(summary(pool(MI.EOC.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)


## Method VI: Exposure-incomplete Multiple Imputation
load("VAHCS.CleanData.Rda")
# Manipulate predicted matrix
dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
dat$thc.alc <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$alc)))
dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
predmat <- make.predictorMatrix(dat)
predmat["cmd", "thc.cmd"] <- 0
predmat["alc", "thc.alc"] <- 0
predmat["zcmd", "thc.zcmd"] <- 0
# Imputating
MI.EI <- mice(dat, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))","~I(as.numeric(as.character(thc))*as.numeric(as.character(alc)))","~I(as.numeric(as.character(thc))*zcmd)"), vis = "monotone", predictorMatrix = predmat, print = F)
# Diagnostic check
summary(mice::complete(MI.EI)) 
# Convergence check
plot(MI.EI)
# Imputated values range check
stripplot(MI.EI)
# Pooling
MI.EI.fit <- with(MI.EI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.EI.ACE <- round(summary(pool(MI.EI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)


## Method VII: Exposure-grouping Multiple Imputation
load("VAHCS.CleanData.Rda")
dat0 <- subset(dat, thc == 0)
dat1 <- subset(dat, thc == 1)
# Manipulate predicted matrix
predmat <- mice(data = dat, maxit = 0)[["predictorMatrix"]]  
predmat[,"thc"] <- 0
# Imputating
MI.EG.0 <- mice(data = dat0, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm"), vis = "monotone", predictorMatrix = predmat, print = F)
MI.EG.1 <- mice(data = dat1, m = 100, maxit = 10, method = c("norm","","","","logreg","logreg","","norm"), vis = "monotone", predictorMatrix = predmat, print = F)
# Merging
MI.EG <- rbind(MI.EG.0, MI.EG.1)
# Diagnostic check
summary(mice::complete(MI.EG)) 
# Convergence check
plot(MI.EG)
# Imputated values range check
stripplot(MI.EG)
# Pooling
MI.EG.fit <- with(MI.EG, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
MI.EG.ACE <- round(summary(pool(MI.EG.fit), conf.int = TRUE, conf.level = 0.95)[7,-c(1,4:5)], 3)


## Comparing
Case <- rbind(CCA.ACE, MI.NI.ACE, MI.EO.ACE, MI.EC.ACE, MI.EOC.ACE, MI.EI.ACE, MI.EG.ACE) %>% 
  cbind(., c("CCA","MI-NI","MI-E×O","MI-E×C","MI-E×OC","MI-E×I","MI-EG")) %>% 
  `colnames<-`(c("Estimate", "Std.Err", "p", "CI.low", "CI.upp", "Method"))
Case$Method <- factor(Case$Method, levels = c("CCA","MI-NI","MI-E×O","MI-E×C","MI-E×OC","MI-E×I","MI-EG"))
save(Case, file = "CaseStudy.Rda")

rm(list=setdiff(ls(), "Case"))

Case$Method <- c("CCA","NI","E×O","E×C","E×OC","E×I","EG")
Case$Method <- factor(Case$Method, levels = c("CCA","NI","E×O","E×C","E×OC","E×I","EG"))

Case.plot <- function(result){
  ggplot(result, aes(x = Method, y = Estimate, group = Method, color = Method)) + 
    geom_point(position = position_dodge(width = 0.5), size = 1.5) + 
    geom_errorbar(aes(ymin = CI.low, ymax = CI.upp), width = 0.3, position = position_dodge(width = 0.5), size = 0.5) +
    scale_color_nejm(drop = F) + theme_classic() + 
    theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 15), legend.position = "none") +
    labs(x = "Method", y = "Mean Difference") 
}

Case %>% Case.plot(.) 


#################################################################################
#  source.r                                                                     # 
#  functions and parameters for coef+thc+ACE+missing.r                          #
#################################################################################
library(boot) 
library(nnet) 
library('MASS')
library("mice")
library(car)
library(magrittr) 
library(dplyr)
library(plyr)
library(ggrepel)
library(tidyr)
library(purrr)

load("VAHCS.CleanData.Rda")
load("simu.coef.Rda")
load("simu.coef.10.Rda")
load("simu.coef.30.Rda")
load("simu.coef.50.Rda")
pari.p <- mean(as.numeric(dat$pari)) - 1

# function: get coeffivicients for generation
coeff <- function(i){
  return(t(as.matrix(coef(i))))
}
#function: generate exposure and missingness idicator for alc
gen.thcMalc <- function(n, seed, coef.thc) {
  set.seed(seed)
  unit1 <- rep(1, n)
  sim <- data.frame(unit1)
  sim$age <- rnorm(n, mean = 0, sd = 1)
  sim$pari <- rbinom(n, 1, pari.p)
  sim$pard <- rbinom(n, 1, inv.logit(coef.thc$b.pard[!is.na(coef.thc$b.pard)] %*% t(sim)))
  sim$asb <- rbinom(n, 1, inv.logit(coef.thc$b.asb[!is.na(coef.thc$b.asb)] %*% t(sim)))
  sim$cmd <- rbinom(n, 1, inv.logit(coef.thc$b.cmd[!is.na(coef.thc$b.cmd)] %*% t(sim)))
  sim$alc <- rbinom(n, 1, inv.logit(coef.thc$b.alc[!is.na(coef.thc$b.alc)] %*% t(sim)))
  sim$thc <- rbinom(n, 1, inv.logit(coef.thc$b.thc[!is.na(coef.thc$b.thc)] %*% t(sim)))
  sim$Malc <- rbinom(n, 1, inv.logit(coef.thc$b.Malc[!is.na(coef.thc$b.Malc)] %*% t(sim)[c(1,2,8),]))
  sim <- sim[,-1]
  sim[2:6] <- lapply(sim[2:6], as.factor)
  return(sim)
}
#function: generate missingness indicators
gen.M <- function(n, seed, coef.thc) {
  set.seed(seed)
  unit1 <- rep(1, n)
  sim <- data.frame(unit1)
  sim$age <- rnorm(n, mean = 0, sd = 1)
  sim$pari <- rbinom(n, 1, pari.p)
  sim$pard <- rbinom(n, 1, inv.logit(coef.thc$b.pard[!is.na(coef.thc$b.pard)] %*% t(sim)))
  sim$asb <- rbinom(n, 1, inv.logit(coef.thc$b.asb[!is.na(coef.thc$b.asb)] %*% t(sim)))
  sim$cmd <- rbinom(n, 1, inv.logit(coef.thc$b.cmd[!is.na(coef.thc$b.cmd)] %*% t(sim)))
  sim$alc <- rbinom(n, 1, inv.logit(coef.thc$b.alc[!is.na(coef.thc$b.alc)] %*% t(sim)))
  sim$thc <- rbinom(n, 1, inv.logit(coef.thc$b.thc[!is.na(coef.thc$b.thc)] %*% t(sim)))
  # Missingness patten b: three incomplete variables
  sim$Malc <- rbinom(n, 1, inv.logit(coef.thc$b.Malc[!is.na(coef.thc$b.Malc)] %*% t(sim)[c(1,2,8),]))
  
  sim$Mcmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.A[!is.na(coef.thc$b.Mcmd.A)] %*% t(sim)[c(1,2,9,8),]))
  sim$Mcmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.B[!is.na(coef.thc$b.Mcmd.B)] %*% t(sim)[c(1,2,9,8,6),]))
  sim$Mcmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.C[!is.na(coef.thc$b.Mcmd.C)] %*% rbind(t(sim)[c(1,2,9,8,6),], sim$thc * sim$cmd)))
  
  sim$b.b.Mzcmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.A[!is.na(coef.thc$b.b.Mzcmd.A)] %*% t(sim)[c(1,2,9,10,8),]))
  sim$b.b.Mzcmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.B[!is.na(coef.thc$b.b.Mzcmd.B)] %*% t(sim)[c(1,2,9,11,8,6),]))
  sim$b.b.Mzcmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.C[!is.na(coef.thc$b.b.Mzcmd.C)] %*% rbind(t(sim)[c(1,2,9,12,8,6),], sim$thc * sim$cmd)))
  # Missingness patten a: one incomplete variable
  sim$b.a.Mzcmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.A[!is.na(coef.thc$b.a.Mzcmd.A)] %*% t(sim)[c(1,2,8),]))
  sim$b.a.Mzcmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.B[!is.na(coef.thc$b.a.Mzcmd.B)] %*% t(sim)[c(1,2,8,6),]))
  sim$b.a.Mzcmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.C[!is.na(coef.thc$b.a.Mzcmd.C)] %*% rbind(t(sim)[c(1,2,8,6),], sim$thc * sim$cmd)))
  
  sim <- sim[,-1]
  sim[2:6] <- lapply(sim[2:6], as.factor)
  return(sim)
}
# function: generate outcomes
gen.zcmd <- function(n, seed, coef.thc) {
  set.seed(seed)
  unit1 <- rep(1, n)
  sim <- data.frame(unit1)
  sim$age <- rnorm(n, mean = 0, sd = 1)
  sim$pari <- rbinom(n, 1, pari.p)
  sim$pard <- rbinom(n, 1, inv.logit(coef.thc$b.pard[!is.na(coef.thc$b.pard)] %*% t(sim)))
  sim$asb <- rbinom(n, 1, inv.logit(coef.thc$b.asb[!is.na(coef.thc$b.asb)] %*% t(sim)))
  sim$cmd <- rbinom(n, 1, inv.logit(coef.thc$b.cmd[!is.na(coef.thc$b.cmd)] %*% t(sim)))
  sim$alc <- rbinom(n, 1, inv.logit(coef.thc$b.alc[!is.na(coef.thc$b.alc)] %*% t(sim)))
  sim$thc <- rbinom(n, 1, inv.logit(coef.thc$b.thc[!is.na(coef.thc$b.thc)] %*% t(sim)))
  # Outcome Scenario I: no interaction
  sim$I <- rnorm(n, mean = (coef.thc$b.I[!is.na(coef.thc$b.I)]) %*% t(sim)[-c(2),], sd = 1) 
  # Outcome Scenario II: interaction between confounder 'cmd' and the exposure 'thc'
  sim$II <- rnorm(n, mean = coef.thc$b.II[!is.na(coef.thc$b.II)] %*% rbind(t(sim)[-c(2,9),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario III: doubled interaction between confounder 'cmd' and the exposure 'thc'
  sim$III <- rnorm(n, mean = coef.thc$b.III[!is.na(coef.thc$b.III)] %*% rbind(t(sim)[-c(2,9:10),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario IV: tripled interaction between confounder 'cmd' and the exposure 'thc'
  sim$IV <- rnorm(n, mean = coef.thc$b.IV[!is.na(coef.thc$b.IV)] %*% rbind(t(sim)[-c(2,9:11),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario V: positive interaction between confounder 'cmd' and the exposure 'thc'
  sim$V <- rnorm(n, mean = coef.thc$b.IV[!is.na(coef.thc$b.IV)] %*% rbind(t(sim)[-c(2,9:12),], sim$thc * sim$cmd), sd = 1)
  # # Outcome Scenario VI: doubled positive interaction between confounder 'cmd' and the exposure 'thc'
  sim$VI <- rnorm(n, mean = coef.thc$b.V[!is.na(coef.thc$b.V)] %*% rbind(t(sim)[-c(2,9:13),], sim$thc * sim$cmd), sd = 1)
  # # Outcome Scenario VII: tripled interaction between confounder 'cmd' and the exposure 'thc'
  sim$VII <- rnorm(n, mean = coef.thc$b.VII[!is.na(coef.thc$b.VII)] %*% rbind(t(sim)[-c(2,9:14),], sim$thc * sim$cmd), sd = 1)
  sim$I %<>% `-`(mean(sim$I)); sim$II %<>% `-`(mean(sim$II)); sim$III %<>% `-`(mean(sim$III)); sim$IV %<>% `-`(mean(sim$IV)); sim$V %<>% `-`(mean(sim$V)); sim$VI %<>% `-`(mean(sim$VI)); sim$VII %<>% `-`(mean(sim$VII))
  sim <- sim[,-1]
  sim[2:7] <- lapply(sim[2:7], as.factor)
  return(sim)
}
# function: convert numeric value to NA
na.n <- function(variable, indicator) {
  variable <- ifelse(indicator == 1, NA, variable)
  return(variable)
}
# function: convert factor value to NA
na.f <- function(variable, indicator) {
  variable <- ifelse(indicator == 1, NA, variable) %>% `-`(1) %>% as.factor(.)
  return(variable)
}
# function: generate complete dataset
gen <- function(n, seed, coef.thc) {
  set.seed(seed)
  unit1 <- rep(1, n)
  sim <- data.frame(unit1)
  sim$age <- rnorm(n, 0, 1)
  sim$pari <- rbinom(n, 1, pari.p)
  sim$pard <- rbinom(n, 1, inv.logit(coef.thc$b.pard[!is.na(coef.thc$b.pard)] %*% t(sim)))
  sim$asb <- rbinom(n, 1, inv.logit(coef.thc$b.asb[!is.na(coef.thc$b.asb)] %*% t(sim)))
  sim$cmd <- rbinom(n, 1, inv.logit(coef.thc$b.cmd[!is.na(coef.thc$b.cmd)] %*% t(sim)))
  sim$alc <- rbinom(n, 1, inv.logit(coef.thc$b.alc[!is.na(coef.thc$b.alc)] %*% t(sim)))
  # Exposure proportation adjusted
  sim$thc <- rbinom(n, 1, inv.logit(coef.thc$b.thc[!is.na(coef.thc$b.thc)] %*% t(sim)))
  
  # Outcome scenarios
  # Outcome Scenario 1: no interaction
  sim$I <- rnorm(n, mean = (coef.thc$b.I[!is.na(coef.thc$b.I)]) %*% t(sim)[c(1,3:8),], sd = 1) 
  # Outcome Scenario 2: interaction with ratio -0.25
  sim$II <- rnorm(n, mean = coef.thc$b.II[!is.na(coef.thc$b.II)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario 3: interaction with ratio -0.5
  sim$III <- rnorm(n, mean = coef.thc$b.III[!is.na(coef.thc$b.III)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario 4: interaction with ratio -0.75
  sim$IV <- rnorm(n, mean = coef.thc$b.IV[!is.na(coef.thc$b.IV)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario 5: interaction with ratio 0.25
  sim$V <- rnorm(n, mean = coef.thc$b.V[!is.na(coef.thc$b.V)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario 6: interaction with ratio 0.5
  sim$VI <- rnorm(n, mean = coef.thc$b.VI[!is.na(coef.thc$b.VI)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  # Outcome Scenario 7: interaction with ratio 0.75
  sim$VII <- rnorm(n, mean = coef.thc$b.VII[!is.na(coef.thc$b.VII)] %*% rbind(t(sim)[c(1,3:8),], sim$thc * sim$cmd), sd = 1)
  sim$I %<>% `-`(mean(sim$I)); sim$II %<>% `-`(mean(sim$II)); sim$III %<>% `-`(mean(sim$III)); sim$IV %<>% `-`(mean(sim$IV)); sim$V %<>% `-`(mean(sim$V)); sim$VI %<>% `-`(mean(sim$VI)); sim$VII %<>% `-`(mean(sim$VII))
  
  # Missingness scenarios
  # a. zcmd incomplete only
  sim$Ma.zcmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.A[!is.na(coef.thc$b.a.Mzcmd.A)] %*% t(sim)[c(1,2,8),]))
  sim$Ma.zcmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.B[!is.na(coef.thc$b.a.Mzcmd.B)] %*% t(sim)[c(1,2,8,6),]))
  sim$Ma.zcmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.a.Mzcmd.C[!is.na(coef.thc$b.a.Mzcmd.C)] %*% rbind(t(sim)[c(1,2,8,6),], sim$thc * sim$cmd)))
  # b. alc, cmd and zcmd incomplete
  sim$Mb.alc <- rbinom(n, 1, inv.logit(coef.thc$b.Malc[!is.na(coef.thc$b.Malc)] %*% t(sim)[c(1,2,8),]))
  
  sim$Mb.cmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.A[!is.na(coef.thc$b.Mcmd.A)] %*% t(sim)[c(1,2,19,8),]))
  sim$Mb.cmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.B[!is.na(coef.thc$b.Mcmd.B)] %*% t(sim)[c(1,2,19,8,6),]))
  sim$Mb.cmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.Mcmd.C[!is.na(coef.thc$b.Mcmd.C)] %*% rbind(t(sim)[c(1,2,19,8,6),], sim$thc * sim$cmd)))
  
  sim$Mb.zcmd.A <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.A[!is.na(coef.thc$b.b.Mzcmd.A)] %*% t(sim)[c(1,2,19,20,8),]))
  sim$Mb.zcmd.B <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.B[!is.na(coef.thc$b.b.Mzcmd.B)] %*% t(sim)[c(1,2,19,21,8,6),]))
  sim$Mb.zcmd.C <- rbinom(n, 1, inv.logit(coef.thc$b.b.Mzcmd.C[!is.na(coef.thc$b.b.Mzcmd.C)] %*% rbind(t(sim)[c(1,2,19,22,8,6),], sim$thc * sim$cmd)))
  
  sim <- sim[,-1]
  sim[2:7] <- lapply(sim[2:7], as.factor)
  sim[15:24] <- lapply(sim[15:24], as.factor)
  return(sim)
}
# funtion: resort into incomplete datasets
resort <- function(i, gen.list){
  lapply(c(1:n.sim), function(j){
    comp <- gen.list[[j]]
    # Outcome I, Missing pattern a, Missingness secnario A-C
    {if(i == 1) {cbind(comp[,1:7], na.n(comp[,8], comp[,15]))}
      else if(i == 2) {cbind(comp[,1:7], na.n(comp[,8], comp[,16]))}
      else if(i == 3) {cbind(comp[,1:7], na.n(comp[,8], comp[,17]))}
      # Outcome II, Missing pattern a, Missingness secnario A-C
      else if(i == 4) {cbind(comp[,1:7], na.n(comp[,9], comp[,15]))}
      else if(i == 5) {cbind(comp[,1:7], na.n(comp[,9], comp[,16]))}
      else if(i == 6) {cbind(comp[,1:7], na.n(comp[,9], comp[,17]))}
      # Outcome III, Missing pattern a, Missingness secnario A-C
      else if(i == 7) {cbind(comp[,1:7], na.n(comp[,10], comp[,15]))}
      else if(i == 8) {cbind(comp[,1:7], na.n(comp[,10], comp[,16]))}
      else if(i == 9) {cbind(comp[,1:7], na.n(comp[,10], comp[,17]))}
      # Outcome IV, Missing pattern a, Missingness secnario A-C
      else if(i == 10) {cbind(comp[,1:7], na.n(comp[,11], comp[,15]))}
      else if(i == 11) {cbind(comp[,1:7], na.n(comp[,11], comp[,16]))}
      else if(i == 12) {cbind(comp[,1:7], na.n(comp[,11], comp[,17]))}
      # Outcome V, Missing pattern a, Missingness secnario A-C
      else if(i == 13) {cbind(comp[,1:7], na.n(comp[,12], comp[,15]))}
      else if(i == 14) {cbind(comp[,1:7], na.n(comp[,12], comp[,16]))}
      else if(i == 15) {cbind(comp[,1:7], na.n(comp[,12], comp[,17]))}
      # Outcome VI, Missing pattern a, Missingness secnario A-C
      else if(i == 16) {cbind(comp[,1:7], na.n(comp[,13], comp[,15]))}
      else if(i == 17) {cbind(comp[,1:7], na.n(comp[,13], comp[,16]))}
      else if(i == 18) {cbind(comp[,1:7], na.n(comp[,13], comp[,17]))}
      # Outcome VII, Missing pattern a, Missingness secnario A-C
      else if(i == 19) {cbind(comp[,1:7], na.n(comp[,14], comp[,15]))}
      else if(i == 20) {cbind(comp[,1:7], na.n(comp[,14], comp[,16]))}
      else if(i == 21) {cbind(comp[,1:7], na.n(comp[,14], comp[,17]))}
      
      # Outcome I, Missing pattern b, Missingness secnario A-C
      else if(i == 22) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,8], comp[,22]))}
      else if(i == 23) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,8], comp[,23]))}
      else if(i == 24) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,8], comp[,24]))}
      # Outcome II, Missing pattern b, Missingness secnario A-C
      else if(i == 25) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,9], comp[,22]))}
      else if(i == 26) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,9], comp[,23]))}
      else if(i == 27) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,9], comp[,24]))}
      # Outcome III, Missing pattern b, Missingness secnario A-C
      else if(i == 28) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,10], comp[,22]))}
      else if(i == 29) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,10], comp[,23]))}
      else if(i == 30) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,10], comp[,24]))}
      # Outcome IV, Missing pattern b, Missingness secnario A-C
      else if(i == 31) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,11], comp[,22]))}
      else if(i == 32) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,11], comp[,23]))}
      else if(i == 33) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,11], comp[,24]))}
      # Outcome V, Missing pattern b, Missingness secnario A-C
      else if(i == 34) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,12], comp[,22]))}
      else if(i == 35) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,12], comp[,23]))}
      else if(i == 36) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,12], comp[,24]))} 
      # Outcome VI, Missing pattern b, Missingness secnario A-C
      else if(i == 37) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,13], comp[,22]))}
      else if(i == 38) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,13], comp[,23]))}
      else if(i == 39) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,13], comp[,24]))} 
      # Outcome VII, Missing pattern b, Missingness secnario A-C
      else if(i == 40) {cbind(comp[,1:4], na.f(comp[,5], comp[,19]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,14], comp[,22]))}
      else if(i == 41) {cbind(comp[,1:4], na.f(comp[,5], comp[,20]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,14], comp[,23]))}
      else if(i == 42) {cbind(comp[,1:4], na.f(comp[,5], comp[,21]), na.f(comp[,6], comp[,18]), comp[,7], na.n(comp[,14], comp[,24]))} } %>% 
      `names<-`(c("age","pari","pard","asb","cmd","alc","thc","zcmd"))
  }) 
}


#################################################################################
#  coef - save regression coefficients and missingness coefficients             #
#################################################################################
rm(list=ls())
setwd("~/R/Project No.1")
load("VAHCS.CleanData.Rda")

# function: read the regression coefficients
coeff <- function(i){
  return(t(as.matrix(coef(i))))
}

## Regression coefficients from VAHCS
b.pard <- coeff(glm(pard ~ age + pari, data = dat, family = binomial(link = "logit")))
b.asb <- coeff(glm(asb ~ age + pari + pard, data = dat, family = binomial(link = "logit")))
b.cmd <- coeff(glm(cmd ~ age + pari + pard + asb, data = dat, family = binomial(link = "logit")))
b.alc <- coeff(glm(alc ~ age + pari + pard + asb + cmd, data = dat, family = binomial(link = "logit")))
b.thc <- coeff(glm(thc ~ age + pari + pard + asb + cmd + alc, data = dat, family = binomial(link = "logit")))
# Outcome Scenario I: 'zcmd' depends on all confounders and the exposure, without interaction
b.I <- coeff(glm(zcmd ~ pari + pard + asb + cmd + alc + thc, data = dat, family = gaussian(link = "identity")))
# Outcome Scenario II: 'zcmd' depends on all confounders, the exposure, and interaction between 'thc' and 'cmd'
b.II <- coeff(glm(zcmd ~ pari + pard + asb + cmd + alc + thc + thc : cmd, data = dat, family = gaussian(link = "identity")))
# Outcome Scenario III: 'zcmd' depends on all confounders, the exposure, and doubled interaction between 'thc' and 'cmd'
b.III <- c(b.II[1:7], b.II[8]*2) 
# Outcome Scenario IV: 'zcmd' depends on all confounders, the exposure, and tripled interaction between 'thc' and 'cmd'
b.IV <- c(b.II[1:7], b.II[8]*3) 
# Outcome Scenario V: 'zcmd' depends on all confounders, the exposure, and positive interaction between 'thc' and 'cmd'
b.V <- c(b.II[1:7], -b.II[8])
# Outcome Scenario VI: 'zcmd' depends on all confounders, the exposure, and doubled positive interaction between 'thc' and 'cmd'
b.VI <- c(b.II[1:7], -b.II[8]*2) 
# Outcome Scenario VII: 'zcmd' depends on all confounders, the exposure, and tripled positive interaction between 'thc' and 'cmd'
b.VII <- c(b.II[1:7], -b.II[8]*3) 

# Missingness indicators
dat$Malc <- as.factor(is.na(dat$alc))
dat$Mcmd <- as.factor(is.na(dat$cmd))
dat$Mzcmd <- as.factor(is.na(dat$zcmd))

# One incomeplete variable: zcmd, under three missingness scenarios A-C:
b.a.Mzcmd.A <- cbind(coeff(glm(Mzcmd ~ age, data = dat, family = binomial(link = "logit"))), log(3))
b.a.Mzcmd.B <- cbind(b.a.Mzcmd.A, log(3))
b.a.Mzcmd.C <- cbind(b.a.Mzcmd.A, log(3), log(2))

# Three incomeplete variables: alc, cmd, zcmd, under three missingness scenarios A-C:
b.Malc <- coeff(glm(Malc ~ age + thc, data = dat, family = binomial(link = "logit")))

b.Mcmd.A <- cbind(coeff(glm(Mcmd ~ age + Malc, data = dat, family = binomial(link = "logit"))), log(3))
b.Mcmd.B <- cbind(b.Mcmd.A, log(3))
b.Mcmd.C <- cbind(b.Mcmd.A, log(3), log(2))

b.b.Mzcmd.A <- cbind(coeff(glm(Mzcmd ~ age + Malc + Mcmd, data = dat, family = binomial(link = "logit"))), log(3))
b.b.Mzcmd.B <- cbind(b.b.Mzcmd.A, log(3))
b.b.Mzcmd.C <- cbind(b.b.Mzcmd.A, log(3), log(2))


coef <- list(b.pard, b.asb, b.cmd, b.alc, b.thc, b.I, b.II, b.III, b.IV, b.V, b.VI, b.VII, b.a.Mzcmd.A, b.a.Mzcmd.B, b.a.Mzcmd.C, b.Malc, b.Mcmd.A, b.Mcmd.B, b.Mcmd.C, b.b.Mzcmd.A, b.b.Mzcmd.B, b.b.Mzcmd.C)
coef <- as.data.frame(sapply(coef, '[', seq(max(sapply(coef, length)))))
names(coef)[1:22] <- c("b.pard", "b.asb", "b.cmd", "b.alc", "b.thc", "b.I", "b.II", "b.III", "b.IV", "b.V", "b.VI", "b.VII", "b.a.Mzcmd.A", "b.a.Mzcmd.B", "b.a.Mzcmd.C", "b.Malc", "b.Mcmd.A", "b.Mcmd.B", "b.Mcmd.C", "b.b.Mzcmd.A", "b.b.Mzcmd.B", "b.b.Mzcmd.C")
save(coef, file = "simu.coef.Rda")


#################################################################################
#  thc - Adjust exposure proportion                                             #
#################################################################################
rm(list=ls())
setwd("~/R/Project No.1")
load("simu.coef.Rda")

coef.10 <- coef; coef.30 <- coef; coef.50 <- coef
save(coef.10, file = "simu.coef.10.Rda"); save(coef.30, file = "simu.coef.30.Rda"); save(coef.50, file = "simu.coef.50.Rda")

rm(list=ls())
cl <- makeCluster(32)
parLapply(cl, c(1:32), function(x) {source(file = "source.r")})
system.time(
  thc.prop <- parLapply(cl, seq(0,4,0.001), function(delta){
    adj.coef <- coef
    adj.coef[1,5] <- adj.coef[1,5] + delta
    lapply(seq(1000,30000,1000), function(seed){
      mean(gen.thcMalc(1000000, seed, adj.coef)[,7]) 
    }) %>% ldply(., data.frame) %>% 
      # Three exposure prevalent levels: 10%, 30% and 50%
      apply(., 2, mean) %>% cbind(adj.coef[1,5], abs(. - 0.1), abs(. - 0.3), abs(. - 0.5))
  })  %>% ldply(., data.frame)
)
stopCluster(cl)

coef.10[1,5] <- thc.prop[thc.prop[,3]==min(thc.prop[,3]),2]; coef.30[1,5] <- thc.prop[thc.prop[,4]==min(thc.prop[,4]),2]; coef.50[1,5] <- thc.prop[thc.prop[,5]==min(thc.prop[,5]),2]
save(coef.10, file = "simu.coef.10.Rda"); save(coef.30, file = "simu.coef.30.Rda"); save(coef.50, file = "simu.coef.50.Rda")


#################################################################################
#  coefMiss - Adjust missingness proportion                                     #
#################################################################################
rm(list=ls())
setwd("~/R/Project No.1")
cl <- makeCluster(32)
parLapply(cl, c(1:32), function(x) {source(file = "source.r")})

# All alc, cmd and zcmd are incomplete
# Adjust missingness proprotion for alc
Malc.prop <- parLapply(cl, seq(0, 0.1, 0.001), function(delta){
  n <- 1000000
  adj.coef <- cbind(coef.10, coef.30, coef.50)
  adj.coef[1,c(16,38,60)] <- coef[1,16] + c(-0.5,-0.4,-0.3) + delta
  lapply(seq(1000,30000,1000), function(seed){
    cbind(gen.thcMalc(n, seed, adj.coef[,c(1:22)])[,8], gen.thcMalc(n,seed,adj.coef[,c(23:44)])[,8], gen.thcMalc(n, seed, adj.coef[,c(45:66)])[,8]) %>% 
      # Missingness proportion for alc: 10%
      apply(., 2, function(Malc){abs(mean(Malc)-0.1)}) %>% t(.)
  }) %>% ldply(., data.frame) %>% 
    apply(., 2, mean) %>% t(.) %>% cbind(adj.coef[1,c(16,38,60)], .)
}) %>% ldply(., data.frame) 
stopCluster(cl)

coef.10[1,16] <- Malc.prop[Malc.prop[,4]==min(Malc.prop[,4]),1]; coef.30[1,16] <- Malc.prop[Malc.prop[,5]==min(Malc.prop[,5]),2]; coef.50[1,16] <- Malc.prop[Malc.prop[,6]==min(Malc.prop[,6]),3]
save(coef.10, file = "simu.coef.10.Rda"); save(coef.30, file = "simu.coef.30.Rda"); save(coef.50, file = "simu.coef.50.Rda")

# function: adjust missingness proportion for cmd
source(file = "source.r")
Mcmd.prop <- function(n, coef.thc, C){
  lapply(seq(0, 0.1, 0.001), function(delta){
    adj.coef <- coef.thc
    adj.coef[1, 17:19] <- coef.thc[1, 17:19] + C + delta
    lapply(seq(1000, 30000, 1000), function(seed){
      # Missingness proportion for cmd: 10%
      adj.coef %>% gen.M(n, seed, .) %>% .[, 9:11] %>% lapply(., function(Mcmd){abs(mean(Mcmd) - 0.1)})
    }) %>% ldply(., data.frame) %>% apply(., 2, mean) %>% t(.) %>% cbind(adj.coef[1, 17:19], .)
  }) %>% ldply(., data.frame) %>% {
    lapply(c(1:3), function(i){.[.[,i+3]==min(.[,i+3]),i]})
  }
}

# function: adjust missingness proportion for zcmd in pattern b
Mb.zcmd.prop <- function(n, coef.thc, C){
  lapply(seq(0, 0.1, 0.001), function(delta){
    adj.coef <- coef.thc
    adj.coef[1, 20:22] <- coef.thc[1, 20:22] + C + delta
    lapply(seq(1000, 30000, 1000), function(seed){
      # Missingness proportion for zcmd, pattern b: 20%
      adj.coef %>% gen.M(n, seed, .) %>% .[, 12:14] %>% lapply(., function(Mcmd){abs(mean(Mcmd) - 0.2)}) 
    }) %>% ldply(., data.frame) %>% apply(., 2, mean) %>% t(.) %>% cbind(adj.coef[1, 20:22], .)
  }) %>% ldply(., data.frame) %>% {
    lapply(c(1:3), function(i){.[.[,i+3]==min(.[,i+3]),i]})
  }
}

# function: adjust missingness proportion for zcmd in pattern a
Ma.zcmd.prop <- function(n, coef.thc, C){
  lapply(seq(0, 0.05, 0.001), function(delta){
    adj.coef <- coef.thc
    adj.coef[1, 13:15] <- coef.thc[1, 13:15] + C + delta
    lapply(seq(1000, 30000, 1000), function(seed){
      # Missingness proportion for zcmd, pattern b: 30%
      adj.coef %>% gen.M(n, seed, .) %>% .[, 15:17] %>% lapply(., function(Mcmd){abs(mean(Mcmd) - 0.3)}) 
    }) %>% ldply(., data.frame) %>% apply(., 2, mean) %>% t(.) %>% cbind(adj.coef[1, 13:15], .)
  }) %>% ldply(., data.frame) %>% {
    lapply(c(1:3), function(i){.[.[,i+3]==min(.[,i+3]),i]})
  }
}

# thc = 10%
coef.10[1,17:19] <- Mcmd.prop(1000000, coef.10, c(0.3, -0.4, -0.5)) 
coef.10[1,20:22] <- Mb.zcmd.prop(1000000, coef.10, c(0.75, 0, -0.1)) 
coef.10[1,13:15] <- Ma.zcmd.prop(1000000, coef.10, c(1.2,0.5,0.45)) 

# thc = 30%
coef.30[1,17:19] <- Mcmd.prop(1000000, coef.30, c(0.1, -0.65, -0.9)) 
coef.30[1,20:22] <- Mb.zcmd.prop(1000000, coef.30, c(0.5, -0.25, -0.55)) 
coef.30[1,13:15] <- Ma.zcmd.prop(1000000, coef.30, c(1, 0.25, 0.05)) 

# thc = 50%
coef.50[1,17:19] <- Mcmd.prop(1000000, coef.50, c(-0.1, -0.85, -1.2)) 
coef.50[1,20:22] <- Ma.zcmd.prop(1000000, coef.50, c(0.25, -0.5, -0.9))  
coef.50[1,13:15] <- Mb.zcmd.prop(1000000, coef.50, c(0.75, 0.05, -0.3)) 

save(coef.10, file = "simu.coef.10.Rda"); save(coef.30, file = "simu.coef.30.Rda"); save(coef.50, file = "simu.coef.50.Rda")


#################################################################################
#  coefOut - Adjust coefficients in the outcome regression                      #
#################################################################################

# function: adjust ACE for outcome scenarios
rm(list=ls())
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "source.r")})

ACE <- function(n, coef.thc, C){
  lapply(seq(0, 0.01, 0.001), function(delta){
    adj.coef <- coef.thc
    adj.coef[7, 6:10] <- coef.thc[7, 6:10] + C + delta
    adj.coef[8, 7:10] <- adj.coef[7, 7:10] %>% `*`(c(-1/4, -1/2, 1/4, 1/2))
    lapply(seq(1000,20000,1000), function(seed){
      dat <- adj.coef %>% gen.zcmd(n, seed, .) 
      lapply(dat[,8:12], function(zcmd){
        coeff(lm(zcmd ~ alc + asb + cmd + pard + pari + thc, data = dat))[,7]
      }) %>% ldply(., data.frame) %>% .[,2] %>% t(.) 
    }) %>% ldply(., data.frame) %>% 
      # True ACE: 0.3
      apply(., 2, function(ACE){abs(mean(ACE)-0.3)}) %>% t(.) %>% cbind(adj.coef[7, 6:10], .)
  }) %>% ldply(., data.frame) %>% {
    lapply(c(1:5), function(i){.[.[,i+5]==min(.[,i+5]),i]})
  }
}

# thc = 10%
coef.10[7,6:12] <- ACE(1000000, coef.10, c(-0.05,-0.05,0.05,0.32,-0.2,-0.2,-0.2)) 
coef.10[8,7:12] <- coef.10[7, 7:12] %>% `*`(c(-1/4, -1/2, -3/4, 1/4, 1/2, 3/4))

# thc = 30%
coef.30[7,6:12] <- ACE(1000000, coef.30, c(-0.05,-0.1,0,0.2,-0.11,-0.2,-0.2)) 
coef.30[8,7:12] <- coef.30[7, 7:12] %>% `*`(c(-1/4, -1/2, -3/4, 1/4, 1/2, 3/4))

# thc = 50%
coef.50[7,6:12] <- ACE(1000000, coef.50, c(-0.05,-0.05,0,0.1,-0.15,-0.2,-0.2)) 
coef.50[8,7:12] <- coef.50[7, 7:12] %>% `*`(c(-1/4, -1/2, -3/4, 1/4, 1/2, 3/4))

save(coef.10, file = "simu.coef.10.Rda"); save(coef.30, file = "simu.coef.30.Rda"); save(coef.50, file = "simu.coef.50.Rda")

# Check the true ACEs
rm(list=ls())
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "source.r")})
test.ACE <- lapply(seq(0, 200000, 100), function(seed){
  dat <- gen.zcmd(1000000, seed, coef.50)
  lapply(dat[,8:14], function(zcmd){
    summary(lm(zcmd ~ pari + pard + asb + cmd + alc + thc, data = dat))[["coefficients"]][7,1]
  }) %>% ldply(., data.frame) %>% .[,2] %>% t(.)
}) %>% ldply(., data.frame) %>% apply(., 2, mean)
stopCluster(cl)


#################################################################################
#  SampleSize - Find out the apporate sample sizes                              #
#################################################################################
rm(list=ls())
source(file = "source.r")

sample.size <- lapply(seq(1200, 1400, 100), function(n){
  lapply(seq(0, 200000, 100), function(seed){
    dat <- gen.zcmd(n, seed, coef.10)
    lapply(dat[,8:14], function(zcmd){
      ifelse(summary(lm(zcmd ~ pari + pard + asb + cmd + alc + thc, data = dat))[["coefficients"]][7,4] < 0.05, 1, 0)
    }) %>% ldply(., data.frame) %>% .[,2] %>% t(.)
  }) %>% ldply(., data.frame) %>% apply(., 2, mean) %>% t(.) %>% cbind(n, .)
}) %>% ldply(., data.frame) 


#################################################################################
#  By Exposure Complete Case Proportion - estimate overall complete case        #
#  proprotions in each exposure grpup                                           #
#################################################################################

lapply(seq(0,1000,100), function(seed){
  dat <- gen.M(1000000, seed, coef.thc)
  dat0 <- subset(dat, thc == 0)
  dat1 <- subset(dat, thc == 1)
  cbind(dat0[,c(8,9,12)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat0[,c(8,10,13)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat0[,c(8,11,14)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat1[,c(8,9,12)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat1[,c(8,10,13)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat1[,c(8,11,14)] %>% apply(., 1, sum) %>% '>='(1) %>% mean(.),
        dat0[,15] %>% mean(.),
        dat0[,16] %>% mean(.),
        dat0[,17] %>% mean(.),
        dat1[,15] %>% mean(.),
        dat1[,16] %>% mean(.),
        dat1[,17] %>% mean(.)
  )
}) %>% ldply(., data.frame) %>% apply(., 2, mean)

#################################################################################
#  Generate incomplete data for missingness methods - 50% exposed               #
#################################################################################

n.sim <- 2000
gen.list <- vector("list", n.sim)
gen.list <- lapply(c(1:n.sim), function(i){
  seed <- 100*i
  gen.list[[i]] <- gen(550, seed, coef.50)
})

data.list.50 <- lapply(c(1:42), function(i){
  resort(i, gen.list) 
})

save(data.list.50, file = "incomplete.data.50.Rda")

# Check exposure and complete-case proportions
prop <- function(dat){
  thc.prop <- mean(as.numeric(as.character(dat$thc)))
  comp.prop <- sum(complete.cases(dat))/nrow(dat)
  cbind(thc.prop, comp.prop)
}
load("incomplete.data.50.Rda")
lapply(c(1:42), function(i){
  map(data.list.50[[i]], prop) %>% ldply(., data.frame) %>% apply(., 2, mean)
}) 


#################################################################################
#  methods.r                                                                    # 
#  functions and parameters for simulation.r                                    #
#################################################################################
library(ggrepel, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(boot, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(nnet, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/") 
library('MASS', lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library("mice", lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(car, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/") 
library(dplyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(ggrepel, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(tidyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
library(purrr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/")
load("/home/jiaxin.zhang/R/Project No.1/incomplete.data.50.Rda")
seed <- seq(1,200000,100)

# Methods for handling missing values
# function: Complete-case analysis
CCA <- function(dat, seed){
  CCA <- lm(zcmd ~ pari + pard + asb + cmd + alc + thc, data = dat) 
  as.data.frame(t(round(cbind(summary(CCA)[["coefficients"]],confint(CCA))[7,-3],3)))
}
# function: Whole-cohort multiple imputation
WC.MI <- function(dat, seed){
  WC.MI <- mice(data = dat, m = 30, method = c("","","","","logreg","logreg","","norm"), seed = seed, print = F)
  WC.MI.fit <- with(WC.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  WC.MI.est <- round(summary(pool(WC.MI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)
}
# function: Exposure-outcome multiple imputation
EO.MI <- function(dat, seed){
  dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
  predmat <- make.predictorMatrix(dat)
  predmat["zcmd", "thc.zcmd"] <- 0
  EO.MI <- mice(dat, m = 30, method = c("","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*zcmd)"), predictorMatrix = predmat, seed = seed, print = F)
  EO.MI.fit <- with(EO.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  EO.MI.est <- round(summary(pool(EO.MI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)
}
# function: Exposure-confounder multiple imputation
EC.MI <- function(dat, seed){
  dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
  predmat <- make.predictorMatrix(dat)
  predmat["cmd", "thc.cmd"] <- 0
  EC.MI <- mice(dat, m = 30, method = c("","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))"), predictorMatrix = predmat, seed = seed, print = F)
  EC.MI.fit <- with(EC.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  EC.MI.est <- round(summary(pool(EC.MI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)
}
# function: Exposure-outcome and exposure-confounder multiple imputation
EOC.MI <- function(dat, seed){
  dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
  dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
  predmat <- make.predictorMatrix(dat)
  predmat["cmd", "thc.cmd"] <- 0
  predmat["zcmd", "thc.zcmd"] <- 0
  EOC.MI <- mice(dat, m = 30, method = c("","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))","~I(as.numeric(as.character(thc))*zcmd)"), predictorMatrix = predmat, seed = seed, print = F)
  EOC.MI.fit <- with(EOC.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  EOC.MI.est <- round(summary(pool(EOC.MI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)
}
# function: Exposure-incomplete multiple imputation
EI.MI <- function(dat, seed){
  dat$thc.cmd <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$cmd)))
  dat$thc.alc <- (as.numeric(as.character(dat$thc)))*(as.numeric(as.character(dat$alc)))
  dat$thc.zcmd <- (as.numeric(as.character(dat$thc)))*dat$zcmd
  predmat <- make.predictorMatrix(dat)
  predmat["cmd", "thc.cmd"] <- 0
  predmat["alc", "thc.alc"] <- 0
  predmat["zcmd", "thc.zcmd"] <- 0
  EI.MI <- mice(dat, m = 30, method = c("","","","","logreg","logreg","","norm","~I(as.numeric(as.character(thc))*as.numeric(as.character(cmd)))","~I(as.numeric(as.character(thc))*as.numeric(as.character(alc)))","~I(as.numeric(as.character(thc))*zcmd)"), predictorMatrix = predmat, seed = seed, print = F)
  EI.MI.fit <- with(EI.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  EI.MI.est <- round(summary(pool(EI.MI.fit), conf.int = TRUE, conf.level = 0.95)[7, -c(1,4:5)], 3)
}
# function: Exposure-grouping multiple imputation
EG.MI <- function(dat, seed){
  dat0 <- subset(dat, thc == 0)
  dat1 <- subset(dat, thc == 1)
  predmat <- mice(data = dat, maxit = 0)[["predictorMatrix"]]  
  predmat[,"thc"] <- 0
  # For 10% thc: m=65; For 30% thc: m=50; For 50% thc: m=40
  EG.MI.0 <- mice(data = dat0, m = 40, method = c("","","","","logreg","logreg","","norm"), predictorMatrix = predmat, seed = seed, print = F)
  EG.MI.1 <- mice(data = dat1, m = 40, method = c("","","","","logreg","logreg","","norm"), predictorMatrix = predmat, seed = seed, print = F)
  EG.MI <- rbind(EG.MI.0, EG.MI.1)
  EG.MI.fit <- with(EG.MI, expr = lm(zcmd ~ pari + pard + asb + cmd + alc + thc))
  EG.MI.est <- round(summary(pool(EG.MI.fit), conf.int = TRUE, conf.level = 0.95)[7,-c(1,4:5)], 3)
}

# funtion: performance indicator
perform.ind <- function(inc) {
  Estimate = mean(inc$est, na.rm = T)
  Bias = mean(inc$est, na.rm = T) - 0.3
  PrBias = abs(Bias/0.3*100) # unit: %
  EmpSE = sd(inc$est, na.rm = T)
  ModSE = sqrt(mean(inc$se^2, na.rm = T))
  MCSE = EmpSE/sqrt(2000)
  Coverage = mean(apply(cbind(inc$low, inc$up), 1, function(x){ifelse(x[1]>0.3 | x[2]<0.3, 0, 1)}), na.rm = T) * 100# unit:%
  StBias <- Bias/EmpSE
  Power = mean(inc$p < 0.05, na.rm = T) * 100# unit: %
  return(cbind(Estimate, Bias, PrBias, StBias, MCSE, EmpSE, ModSE, Coverage, Power))
}


#################################################################################
#  simulation.50.r                                                              # 
#  simulation for missingness methods in different exposure prevalents          #
#################################################################################
rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
CCA.50 <- parLapply(cl, c(1:42), function(i){
  map2(data.list.50[[i]], seed, CCA) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "CCA")
stopCluster(cl)
save(CCA.50, file = "/home/jiaxin.zhang/R/Project No.1/CCA.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
WC.MI.50 <- parLapply(cl, c(1:42), function(i){
  map2(data.list.50[[i]], seed, WC.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "WC.MI")
stopCluster(cl)
save(WC.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/WC.MI.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
EO.MI.50 <- parLapply(cl, c(22:42), function(i){
  map2(data.list.50[[i]], seed, EO.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "EO.MI")
stopCluster(cl)
save(EO.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/EO.MI.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
EC.MI.50 <- parLapply(cl, c(1:42), function(i){
  map2(data.list.50[[i]], seed, EC.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "EC.MI")
stopCluster(cl)
save(EC.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/EC.MI.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
EOC.MI.50 <- parLapply(cl, c(22:42), function(i){
  map2(data.list.50[[i]], seed, EOC.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "EOC.MI")
stopCluster(cl)
save(EOC.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/EOC.MI.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
EI.MI.50 <- parLapply(cl, c(22:42), function(i){
  map2(data.list.50[[i]], seed, EI.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "EI.MI")
stopCluster(cl)
save(EI.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/EI.MI.50.Rda")


rm(list=ls())
library(parallel, lib.loc = "/hpc/software/installed/R/3.6.1/lib64/R/library")
library(magrittr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE) 
library(plyr, lib.loc = "/home/jiaxin.zhang/R/x86_64-pc-linux-gnu-library/3.6/", warn.conflicts = FALSE)
cl <- makeCluster(36)
parLapply(cl, c(1:36), function(x) {source(file = "/home/jiaxin.zhang/R/Project No.1/methods.r")})
EG.MI.50 <- parLapply(cl, c(1:42), function(i){
  map2(data.list.50[[i]], seed, EG.MI) %>% ldply(., data.frame) %>% `names<-`(c("est","se","p","low","up")) %>% perform.ind(.)
}) %>% ldply(., data.frame) %>% cbind(., Exposure = "50%", Method = "EG.MI")
stopCluster(cl)
save(EG.MI.50, file = "/home/jiaxin.zhang/R/Project No.1/EG.MI.50.Rda")



#################################################################################
#  simulation.pbs                                                               # 
#  pbs file for simulation using HPC                                            #
#################################################################################
#!/bin/bash

##########################
#                        #
#   The PBS directives   #
#                        #
##########################

# Define the shell in which your jobs should run. Shouldn't really be changed
# unless you have a very specific reason for doing so
#PBS -S /bin/bash

# Define the name for the job
#PBS -N simulation

# Define a new name for the .o file
#PBS -o simulation.o
#
# If you don't want a .o file
# #PBS -o /dev/null

# Define a new name for the .e file
# #PBS -e <new name>
#
# If you don't want a .e file
# #PBS -e /dev/null

# If you want to merge the standard error file into the standard output file
#PBS -j oe

# Defining the wall time for the job
#PBS -l walltime=48:00:00

# Selecting which queue to send the job to
#PBS -q batch

# Defining the amount of memory you require
#PBS -l mem=200GB

# Defining email notifications
## a = notify when job aborts (default)
## b = notify when job begins
## e = notify when job terminates
## n = no mail at all (If specified takes precedence)
# #PBS -m n

# Define the email address to be used in correspondence
#PBS -M jiaxin.zhang@mcri.edu.au


# Define the number of nodes and cores you require
#PBS -l nodes=1:ppn=36

# Define which project i am with
#PBS -A cebu

##########################################
#                                        #
#   Output some useful job information.  #
#                                        #
##########################################

echo ------------------------------------------------------
  echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
  echo PBS: qsub was run on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: temporary directory on node is $TMPDIR
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
  
  runJob()
{
  module load R/3.6.1
  
  Rscript /home/jiaxin.zhang/R/Project\ No.1/simulation.r
  }


# ----- Notifies you that your job terminated early ----- #
earlyTermination()
{
  echo ' '
  echo ' ############ WARNING:  EARLY TERMINATION #############'
  echo ' '
}
trap 'earlyTermination;' 2 9 15
# ------------------------------------------------------- #

# Run the steps required to complete this job
runJob

exit



#################################################################################
#  plots                                                                        #
#################################################################################
library(magrittr) 
library(plyr)
library(ggpubr)
library(ggthemes)
library("ggsci")
library("ggplot2")
library(purrr)

rm(list=ls())
load("CCA.10.Rda");load("CCA.30.Rda");load("CCA.50.Rda")
load("WC.MI.10.Rda");load("WC.MI.30.Rda");load("WC.MI.50.Rda")
load("EO.MI.10.Rda");load("EO.MI.30.Rda");load("EO.MI.50.Rda")
load("EC.MI.10.Rda");load("EC.MI.30.Rda");load("EC.MI.50.Rda")
load("EOC.MI.10.Rda");load("EOC.MI.30.Rda");load("EOC.MI.50.Rda")
load("EI.MI.10.Rda");load("EI.MI.30.Rda");load("EI.MI.50.Rda")
load("EG.MI.10.Rda");load("EG.MI.30.Rda");load("EG.MI.50.Rda")
# function: convert the data form for plotting
trandat <- function(method.thc) {
  if(is.data.frame(method.thc)) { 
  nrow <- nrow(method.thc)
  Outcome <- rep(rep(c("I","II","III","IV","V","VI","VII"), each = 3), 2)[1:nrow]
  ifelse(nrow == 42, {Missing = rep(c("a","b"), each = 21)}, {Missing = rep("b", 21)})
  Missingness <- rep(c("A","B","C"), 14)[1:nrow]
  method.thc <- method.thc %>% 
    cbind(., Outcome, Missing, Missingness) 
  return(method.thc)}
}
result <- lapply(ls(), get) %>% map(., trandat) %>% do.call("rbind", .)

#result$Method <- factor(result$Method, levels = c("CCA","WC.MI","EO.MI","EC.MI","EOC.MI","EI.MI", "EG.MI"))
#result$Methods <- factor(result$Method, levels = c("CCA","WC.MI","EC.MI","EG.MI","EO.MI","EOC.MI","EI.MI"))
levels(result$Method) <- c("CCA","MI-E×C","MI-EG","MI-E×I","MI-E×O","MI-E×OC","MI-NI")
result$Method <- factor(result$Method, levels = c("CCA","MI-NI","MI-E×O","MI-E×C","MI-E×OC","MI-E×I","MI-EG"))
result$Methods <- factor(result$Method, levels = c("CCA","MI-NI","MI-E×C","MI-EG","MI-E×O","MI-E×OC","MI-E×I"))

levels(result$Outcome) <- c("No Interaction","Weak Negative","Moderate Negative","Strong Negative", "Weak Positive", "Moderate Positive", "Strong Positive")
levels(result$Missingness) <- c("A\nexposure","B\nexposure + confounder","C\nexposure × confounder")
result$Outcome <- factor(result$Outcome, levels = c("Strong Negative","Moderate Negative","Weak Negative","No Interaction","Weak Positive", "Moderate Positive", "Strong Positive"))
save(result, file = "simu.result.Rda")
rm(list=setdiff(ls(), "result"))

# function: mean ACE estimate error bar for each exposure proportion
Estimates.MCSE <- function(result){
  ggplot(result, aes(x = Outcome, y = Estimate, group = Method, color = Method)) + 
    geom_hline(data = result, aes(yintercept = 0.3), linetype = "solid", size = 0.2) +
    geom_hline(data = result, aes(yintercept = 0.33), linetype = "dotted", size = 0.2) +
    geom_hline(data = result, aes(yintercept = 0.27), linetype = "dotted", size = 0.2) +
    geom_point(position = position_dodge(width = 0.5), size = 1.5) + 
    geom_errorbar(aes(ymin = Estimate - MCSE, ymax = Estimate + MCSE), width = 0.3, position = position_dodge(width = 0.5), size = 0.5) +
    annotate("rect", xmin = 0.4, xmax = 7.6, ymin = 0.27, ymax = 0.33, alpha = 0.1) + 
    scale_color_nejm(drop = F) + theme_classic() + 
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 20), legend.position = "bottom") +
    facet_grid(Missingness ~ .) +
    guides(color = guide_legend(nrow = 1)) +
    labs(title = expression(paste(sep = "Mean ACE Estimate" %+-% "Monte Carlo Starndard Error")), 
         x = "Outcome Scenarios", y = "Mean Average Causal Effect Estimate", 
         caption = expression(paste("*The gray shade represents 10% bias from true ACE"))) 
}

Estimates <- function(result){
  ggplot(result, aes(x = Outcome, y = Estimate, group = Method, color = Method)) + 
    geom_hline(data = result, aes(yintercept = 0.3), linetype = "solid", size = 0.2) +
    geom_hline(data = result, aes(yintercept = 0.33), linetype = "dotted", size = 0.2) +
    geom_hline(data = result, aes(yintercept = 0.27), linetype = "dotted", size = 0.2) +
    geom_line(size = .5) + 
    annotate("rect", xmin = 0.4, xmax = 7.6, ymin = 0.27, ymax = 0.33, alpha = 0.1) + 
    scale_color_nejm(drop = F) + theme_classic() + 
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 20), legend.position = "bottom") +
    facet_grid(Missingness ~ .) +
    guides(color = guide_legend(nrow = 1)) +
    labs(title = expression(paste(sep = "Mean ACE Estimate")), 
         x = "Outcome Scenarios", y = "Mean Average Causal Effect Estimate", 
         caption = expression(paste("*The gray shade represents 10% bias from true ACE"))) 
}

result[result[, "Exposure"] == "10%" & result[, "Missing"] == "a", ] %>% 
  Estimates.MCSE(.) %>% 
  `+`(labs(subtitle = "Exposure Proportion: 10%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "10%" & result[, "Missing"] == "b", ] %>% 
  Estimates.MCSE(.) %>% 
  `+`(labs(subtitle = "Exposure Proportion: 10%; Incomplete variables: alc, cmd, zcmd"))

result[result[, "Exposure"] == "30%" & result[, "Missing"] == "a", ] %>% 
  Estimates.MCSE(.) %>% 
  `+`(labs(subtitle = "Exposure Proportion: 30%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "30%" & result[, "Missing"] == "b", ] %>% 
  Estimates.MCSE(.) %>% 
 `+`(labs(subtitle = "Exposure Proportion: 30%; Incomplete variables: alc, cmd, zcmd"))

result[result[, "Exposure"] == "50%" & result[, "Missing"] == "a", ] %>%
  Estimates.MCSE(.) %>% 
  `+`(labs(subtitle = "Exposure Proportion: 50%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "50%" & result[, "Missing"] == "b", ] %>%
  Estimates.MCSE(.) %>% 
  `+`(labs(subtitle = "Exposure Proportion: 50%; Incomplete variables: alc, cmd, zcmd"))


result[result[, "Exposure"] == "10%" & result[, "Missing"] == "a", ] %>% 
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 10%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "10%" & result[, "Missing"] == "b", ] %>% 
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 10%; Incomplete variables: alc, cmd, zcmd"))

result[result[, "Exposure"] == "30%" & result[, "Missing"] == "a", ] %>% 
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 30%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "30%" & result[, "Missing"] == "b", ] %>% 
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 30%; Incomplete variables: alc, cmd, zcmd"))

result[result[, "Exposure"] == "50%" & result[, "Missing"] == "a", ] %>%
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 50%; Incomplete variable: zcmd"))
result[result[, "Exposure"] == "50%" & result[, "Missing"] == "b", ] %>%
  Estimates(.) %>%
  `+`(labs(subtitle = "Exposure Proportion: 50%; Incomplete variables: alc, cmd, zcmd"))


## Performance comparison across 27 scenarios
my.labels <- c("CCA","MI-\nWC","MI-\nE×O","MI-\nE×C","MI-\nE×OC","MI-\nE×I","MI-\nEG")

# function: relative bias
tapply(result$PrBias, result$Missing, summary)
ReBias <- function(data){ggplot(data, aes(x = Method, y = Exposure)) + 
    geom_raster(aes(fill = PrBias), hjust = 0.5, vjust = 0.5, interpolate=F) +
    scale_fill_gradient2(name="Bias (%)", low="#005180", mid="#FFFFFFFF", high="#991400", limits = c(-10.01,34), midpoint = 0, breaks=c(-10,0,10,20,30)) +
    scale_y_discrete(name ="Exposure prevalence") +
    scale_x_discrete(labels= my.labels) +
    scale_color_nejm() + theme_classic() + 
    #theme(legend.key.size = unit(1, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.key.size = unit(0.5, "cm"), legend.key.width = unit(5, "cm")) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 20), legend.position = "bottom") +
    facet_grid(Outcome ~ Missingness, labeller = labeller(Outcome = label_wrap_gen(10))) +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(title ="Relative bias")
}
ReBias(result[result[, "Missing"] == "a",]) + labs(subtitle = "Incomplete variable: zcmd", caption = expression(paste("*The range for relative bias is [-7.9, 33.5]"))) 
ReBias(result[result[, "Missing"] == "b",]) + labs(subtitle = "Incomplete variable: alc, cmd, zcmd", caption = expression(paste("*The range for relative bias is [-7.1, 27.1]"))) 


# function: coverage
tapply(result$Coverage, result$Missing, summary)
Coverage <- function(data){ggplot(data, aes(x = Method, y = Exposure)) + 
    geom_raster(aes(fill = Coverage), hjust = 0.5, vjust = 0.5, interpolate=F) +
    #scale_fill_gradientn(name="Bias (%)", colours = c("#005180","#FFFFFFFF","#991400"), midpoint = 95, breaks=c(90, 95, 100), limits = c(88,97)) +
    scale_fill_gradient2(low="#005180", mid="#FFFFFFFF", high="#991400", limits = c(89.99,100), midpoint = 95, breaks=c(90,92.5,95,97.5,100)) +
    scale_y_discrete(name ="Exposure Proportion") +
    scale_color_nejm() + theme_classic() + 
    theme(legend.key.size = unit(0.5, "cm"), legend.key.width = unit(5, "cm")) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 20), legend.position = "bottom") +
    facet_grid(Outcome ~ Missingness, labeller = labeller(Outcome = label_wrap_gen(10))) +
    theme(strip.text.y = element_text(angle = 0)) +
    labs(title ="Coverage")}

Coverage(result[result[, "Missing"] == "a",]) + labs(subtitle = "Incomplete variable: zcmd", caption = expression(paste("*The range for coverage is [90.8, 96.8]"))) 
Coverage(result[result[, "Missing"] == "b",]) + labs(subtitle = "Incomplete variable: alc, cmd, zcmd", caption = expression(paste("*The range for coverage is [91.6, 96.2]"))) 


# function: empirical and model-based standard error
Error <- function(data){ggplot(data, aes(x = Exposure, y = EmpSE, group = Method, color = Method)) + 
    geom_line() +
    geom_errorbar(data = data[abs(data[, "ModSE"]-data[, "EmpSE"]) > 0.005, ], aes(ymin = EmpSE, ymax = ModSE), width = 0.1, size = 0.8) +
    scale_color_nejm(drop = F) + theme_classic() + 
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5), plot.caption = element_text(hjust = 1, size = 20), legend.position = "bottom") +
    facet_grid(Outcome ~ Missingness, labeller = labeller(Outcome = label_wrap_gen(10))) +
    theme(strip.text.y = element_text(angle = 0)) +
    guides(color = guide_legend(nrow = 1)) +
    labs(title = expression(paste(sep ="Empirical SE trends and difference with Model-based SE")), 
         x = "Exposure Proportion", y = "Empirical Standard Error",
         caption = expression(paste("*Error bars represent the difference between EmpSE and ModSE greater than 0.005"))) 
}

Error(result[result[, "Missing"] == "a", ]) +
  labs(subtitle = "Incomplete variable: zcmd")

Error(result[result[, "Missing"] == "b", ]) +
  labs(subtitle = "Incomplete variable: alc, cmd, zcmd")




