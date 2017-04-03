# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 1: Traditional Gravity Estimates
# Code to replicate the results from Table 1 (page 42)
# Date: 2017/03/25


# Function to load (and install packages if necessary) --------------------

check_and_install <- function(package_name){
  if(!require(package_name, character.only = TRUE)){
    install.packages(package_name, character.only = TRUE)
  }
}

# Load packages -----------------------------------------------------------
packages <- c('haven', 'dplyr', 'multiwayvcov', 'lmtest',
              'speedglm', 'Matrix', 'lfe')

load_packages <- sapply(packages, check_and_install)

# load speedglm_cluster function
source('R/speedglm_cluster.R')

# Read and Prepare Data --------------------------------------------------

data <- read_dta('Data/Chapter1Application1.dta')

# Filter years 1986, 1990, 1994, 1998, 2002 and 2006
# Create new variables
data <- data %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  mutate(ln_trade = log(trade),
         ln_DIST = log(DIST)) %>% 
  group_by(exporter, year) %>% 
  mutate(Y = sum(trade),
         ln_Y = log(Y)) %>% 
  group_by(importer, year) %>% 
  mutate(E = sum(trade),
         ln_E = log(E)) 

# Filter: only international trade and trade > 0
filter1 <- data$exporter != data$importer & data$trade > 0

# Column 1 (OLS) ----------------------------------------------------------
# OLS estimation ignoring multilateral resistance terms

fit <- lm(ln_trade ~ ln_DIST + CNTG + LANG + CLNY + ln_Y + ln_E,
          subset = filter1,
          data = data)

# Summary with default standard errors
summary(fit)

# Compute clustered standard errors
vcov_cluster <- cluster.vcov(fit, cluster = data[filter1, "pair_id"])

# Tests on coefficients with clustered standard errors
coeftest(fit, vcov_cluster)

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
  mutate(REM_EXP = sum(DIST *  TotE / E),
         ln_REM_EXP = log(REM_EXP)) %>% 
  group_by(importer, year) %>% 
  mutate(TotYi = sum(Y)) %>% 
  group_by(year) %>% 
  mutate(TotY = max(TotYi)) %>% 
  group_by(importer, year) %>% 
  mutate(REM_IMP = sum(DIST / (Y / TotY)),
         ln_REM_IMP = log(REM_IMP))

# OLS estimation controlling for multilateral resistance terms with remoteness indexes
fit2 <- lm(ln_trade ~ ln_DIST + CNTG + LANG + CLNY + ln_Y + ln_E +
             ln_REM_EXP + ln_REM_IMP,
           subset = filter1,
           data = data)

# Compute clustered standard-errors
vcov_cluster2 <- cluster.vcov(fit2, cluster = data[filter1, "pair_id"])

# Tests on coefficients with clustered standard errors
coeftest(fit2, vcov_cluster2)

# Reset test
resettest(fit2, power = 2)


# Colum 3 (OLS with Fixed Effects) ----------------------------------------

# Create exporter/year e importer/year variables
data <- data %>% 
  mutate(exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year))

# OLS with fixed effects to control for multilateral resistance terms
# The felm function fits a linear model with multiple group fixed effects
# It uses the method of alternating projections
# It is more efficient in time
# see ?felm 
# the last block pair_id defines the cluster specification for standard-errors.

fit3 <- felm(ln_trade ~ ln_DIST + CNTG + LANG + CLNY | (exp_year) + (imp_year) | 0 | pair_id,
     subset = filter1,
     data = data)

summary(fit3)

# Colum 4 (PPML with Fixed Effects) ----------------------------------------
# First solution: glm()
# It is slower, but we can use the cluster.vcov function

fit5 <- glm(trade ~ 1 + ln_DIST + CNTG + LANG + CLNY + exp_year + imp_year,
            family = quasipoisson(),
            subset = exporter != importer,
            data = data)

# Compute clustered standard errors
vcov_cluster5 <- cluster.vcov(fit5,
                              cluster = data[data$exporter != data$importer, "pair_id"],
                              df_correction = FALSE)

# Tests on coefficients with clustered standard errors
# Difference in constant due the reference level
c <- coeftest(fit5, vcov_cluster5)[1:5, ]
round(c, 3)

# Reset test

pred <- predict(fit5)
pred <- pred^2
data[data$exporter != data$importer, "pred"] <- pred
fit.reset <- glm(trade ~ 1 + pred + ln_DIST + CNTG + LANG + CLNY + exp_year + imp_year,
            family = quasipoisson(),
            subset = exporter != importer,
            data = data)

vcov_cluster_reset <- cluster.vcov(fit.reset,
                              cluster = data[data$exporter != data$importer, "pair_id"],
                              df_correction = FALSE)

c.reset <- coeftest(fit.reset, vcov_cluster_reset)[2, ]
round(c.reset, 3)

# Second: speedglm
# The speedglm package gives a faster implementation of glm.
# One disadvantage is that it does provide a function to compute robust/clustered
# standard errors.
# Use cluster.vcov.speedglm

X <- sparse.model.matrix(~ +1 + ln_DIST + CNTG + LANG + CLNY + exp_year + imp_year,
                         data = data[data$exporter != data$importer, ])

fit6 <- speedglm.wfit(X = X, y = data[data$exporter != data$importer, "trade"][[1]],
                      sparse = TRUE,
                      family = quasipoisson(link = log),
                      fitted = TRUE)

fit6.cluster <- speedglm.cluster(fit6, X = X,
                                 y = data[data$exporter != data$importer, "trade"][[1]],
                                 cluster = data[data$exporter != data$importer, "pair_id"][[1]])

round(fit6.cluster, 3)[1:5, ]
