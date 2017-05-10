# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 2: The "distance puzzle" resolved
# Code to replicate the results from Table 2 (page 47)
# Date: 2017/04/05


# Function to load (and install packages if necessary) --------------------

check_and_install <- function(package_name){
  if(!require(package_name, character.only = TRUE)){
    install.packages(package_name, character.only = TRUE)
  }
}

# Load packages -----------------------------------------------------------
packages <- c('haven', 'multiwayvcov', 'lmtest',
              'speedglm', 'Matrix', 'lfe', 'dplyr')

load_packages <- sapply(packages, check_and_install)

# load speedglm_cluster function
source('R/aux_functions.R')

# Read and Prepare Data --------------------------------------------------

data <- read_dta('Data/Chapter1Application2.dta')

# Create new variables
data <- data %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  mutate(exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year),
         year_f = as.factor(year))

dist_year <- model.matrix(~ -1 + log(DIST):year_f,
                          data = data)

dist_year <- data.frame(dist_year)

colnames(dist_year) <- c("logDIST1986", "logDIST1990",
                         "logDIST1994", "logDIST1998",
                         "logDIST2002", "logDIST2006")

data <- bind_cols(data, dist_year)
# OLS with exporter/time and importer/time fixed effects ---------------------

fit <- felm(log(trade) ~ logDIST1986 + logDIST1990 + logDIST1994 +
              logDIST1998 + logDIST2002 + logDIST2006 +
              CNTG + LANG + CLNY |
              exp_year + imp_year | 0 | pair_id,
            data = data,
            subset = trade > 0 & importer != exporter)

summary(fit)

# PPML with exporter/time and importer/time fixed effects --------------------
xvars <- c("logDIST1986", "logDIST1990",
           "logDIST1994", "logDIST1998",
           "logDIST2002", "logDIST2006",
           "CNTG", "LANG", "CLNY")

fit2 <- gravity_ppml3(y = "trade", x = xvars,
                     data = data,
                     fixed_effects = c("exp_year", "imp_year"),
                     robust = TRUE,
                     cluster = "pair_id",
                     subset = data$importer != data$exporter)

summary(fit2)

# PPML with intra trade ------------------------------------------------------

# create intra trade distance variable
data <- data %>% 
  mutate(SMCTRY = ifelse(importer == exporter, 1, 0),
         logDIST_INTRA = log(DIST)*SMCTRY) %>% 
  mutate_each(funs(.*(1 - SMCTRY)), contains("logDIST1"), contains("logDIST2"))

xvars <- c(xvars, 'logDIST_INTRA')

fit3 <- gravity_ppml3(y = "trade", x = xvars,
                     data = data,
                     fixed_effects = c("exp_year", "imp_year"),
                     cluster = "pair_id")

summary(fit3)

# PPML with intra trade and Home Bias   --------------------------------------

# create intra trade distance variable
xvars <- c(xvars, 'SMCTRY')

fit4 <- gravity_ppml3(y = "trade", x = xvars,
                     data = data,
                     fixed_effects = c("exp_year", "imp_year"),
                     cluster = "pair_id")

summary(fit4)

# PPML with intra-national fixed effects

data <- data %>% 
  mutate(intra_pair = ifelse(exporter == importer, exporter, "inter"))

xvars <- xvars[1:9]

fit5 <- gravity_ppml3(y = "trade", x = xvars,
                     data = data,
                     fixed_effects = c("exp_year", "imp_year", "intra_pair"),
                     cluster = "pair_id")

summary(fit5)
