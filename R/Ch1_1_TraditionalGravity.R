# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 1: Traditional Gravity Estimates
# Code to replicate the results from Table 1 (page 42)
# Date: 2021/01/28

# Load packages -----------------------------------------------------------
# packages <- c(
#        "haven", "multiwayvcov", "lmtest",
#        "speedglm", "Matrix", "lfe",
#        "dplyr"
# )

# load_packages <- sapply(packages, check_and_install)

# # load speedglm_cluster function
# source("R/aux_functions.R")

library(haven)
library(dplyr)
library(estimatr)
library(lmtest)
library(fixest)

# Read and Prepare Data --------------------------------------------------

data <- read_dta("Data/Chapter1Application1.dta")

# Filter years 1986, 1990, 1994, 1998, 2002 and 2006
# Create new variables
data <- data %>%
  filter(year %in% seq(1986, 2006, 4)) %>%
  mutate(
    ln_trade = log(trade),
    ln_DIST = log(DIST)
  ) %>%
  group_by(exporter, year) %>%
  mutate(
    Y = sum(trade),
    ln_Y = log(Y)
  ) %>%
  group_by(importer, year) %>%
  mutate(
    E = sum(trade),
    ln_E = log(E)
  )

# Filter: only international trade and trade > 0
filter1 <- data$exporter != data$importer & data$trade > 0

# Column 1 (OLS) ----------------------------------------------------------
# OLS estimation ignoring multilateral resistance terms

fit <- lm_robust(
  ln_trade ~ ln_DIST + CNTG + LANG + CLNY + ln_Y + ln_E,
  data = data %>% filter(exporter != importer, trade > 0),
  clusters = pair_id,
  se_type = "stata"
)

summary(fit)

# Reset test
resettest(fit, power = 2)

# Column 2 (OLS - Remotness indexes) ----------------------------------------

# Create remotness indexes variables
data <- data %>%
  group_by(exporter, year) %>%
  mutate(TotEj = sum(E)) %>%
  group_by(year) %>%
  mutate(TotE = max(TotEj)) %>%
  group_by(exporter, year) %>%
  mutate(
    REM_EXP = sum(DIST * TotE / E),
    ln_REM_EXP = log(REM_EXP)
  ) %>%
  group_by(importer, year) %>%
  mutate(TotYi = sum(Y)) %>%
  group_by(year) %>%
  mutate(TotY = max(TotYi)) %>%
  group_by(importer, year) %>%
  mutate(
    REM_IMP = sum(DIST / (Y / TotY)),
    ln_REM_IMP = log(REM_IMP)
  ) %>% 
  ungroup()

# OLS estimation controlling for multilateral resistance terms with remoteness indexes
fit2 <- lm_robust(
  ln_trade ~ ln_DIST + CNTG + LANG + CLNY + ln_Y + ln_E + ln_REM_EXP + ln_REM_IMP,
  data = data %>% filter(exporter != importer, trade > 0),
  cluster = pair_id,
  se_type = "stata"
)

summary(fit2)

# Reset test
resettest(fit2, power = 2)

# Colum 3 (OLS with Fixed Effects) ----------------------------------------

# Create exporter/year e importer/year variables
data <- data %>%
  mutate(
    exp_year = paste0(exporter, year),
    imp_year = paste0(importer, year)
  )

# OLS with fixed effects to control for multilateral resistance terms
# use the fixest package: Fast and user-friendly estimation of econometric
# models with multiple fixed-effects.

fit3 <- feols(
  ln_trade ~ ln_DIST + CNTG + LANG + CLNY | exp_year + imp_year,
  data = data %>% filter(exporter != importer, trade > 0)
)

summary(fit3, cluster = ~pair_id)

# Colum 4 (PPML with Fixed Effects) ----------------------------------------
# First solution: glm()
# It is slower, but we can use the cluster.vcov function

fit4 <- feglm(
  trade ~ ln_DIST + CNTG + LANG + CLNY | exp_year + imp_year,
  family = "poisson",
  data = data %>% filter(exporter != importer)
)

summary(fit4, cluster = ~ pair_id)

# Reset test
pred2 <- predict(fit4, type = "link")^2
#data[data$exporter != data$importer, "pred_2"] <- pred2
data <- data %>% 
  filter(exporter != importer) %>% 
  mutate(pred_2 = pred2)

fit.reset <- feglm(
  trade ~ pred_2 + ln_DIST + CNTG + LANG + CLNY | exp_year + imp_year,
  family = "poisson",
  data =  data %>% 
    filter(exporter != importer) %>% 
    ungroup()
)

summary(fit.reset, cluster = ~ pair_id)
