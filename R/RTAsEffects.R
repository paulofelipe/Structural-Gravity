# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 3: Reginal trade agreements effects
# Code to replicate the results from Table 3 (page 51)
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

data <- read_dta('Data/Chapter1Application3.dta')

# Create new variables
data <- data %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  mutate(ln_trade = log(trade),
         ln_DIST = log(DIST),
         exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year),
         pair_id2 = ifelse(exporter == importer, "intra", pair_id),
         INTL_BRDR = ifelse(exporter != importer, "inter", pair_id),
         INTL_BRDR2 = ifelse(exporter != importer, 1, 0),
         year = as.factor(year))

# Total trade by pair
data <- data %>% 
  group_by(pair_id) %>% 
  mutate(sum_trade = sum(trade)) %>% 
  ungroup()

# Interaction Border and Year

INTL_BRDR_YEAR <- model.matrix(~ -1 + INTL_BRDR2:year, data = data) %>% 
  data.frame()
data <- bind_cols(data, INTL_BRDR_YEAR)

# Column 1 (OLS) ----------------------------------------------------------
# Traditional Estimates

fit1 <- felm(ln_trade ~ ln_DIST + CNTG + LANG + CLNY + RTA| exp_year + imp_year |
               0 | pair_id,
             data = data,
             subset = trade > 0 & importer != exporter)

# Column 2 (PPML) ----------------------------------------------------------
# Traditional Estimates

fit2 <- gravity_ppml3(y = "trade",
                      x = c("ln_DIST", "CNTG", "LANG", "CLNY", "RTA"),
                      data = data,
                      fixed_effects = c("exp_year", "imp_year"),
                      robust = TRUE,
                      subset = data$exporter != data$importer,
                      cluster = 'pair_id')

summary(fit2)

# Column 3 (PPML) ----------------------------------------------------------
# Traditional Estimates with Intra-national Trade

fit3 <- gravity_ppml3(y = "trade",
                      x = c("ln_DIST", "CNTG", "LANG", "CLNY", "RTA"),
                      data = data,
                      fixed_effects = c("exp_year", "imp_year", "INTL_BRDR"),
                      robust = TRUE,
                      cluster = 'pair_id')

summary(fit3)

# Column 4 (PPML) ----------------------------------------------------------
# Addressing Endogeneity of RTAs

fit4 <- gravity_ppml3(y = "trade",
                     x = c("RTA"),
                     data = data[data$sum_trade > 0, ],
                     fixed_effects = c("exp_year", "imp_year", "pair_id2"),
                     robust = TRUE,
                     cluster = 'pair_id')

summary(fit4)

# Column 5 (PPML) ----------------------------------------------------------
# Testing potential reverse causality between Trade and RTA

fit5 <- gravity_ppml3(y = "trade",
                      x = c("RTA", "RTA_LEAD4"),
                      data = data,
                      fixed_effects = c("exp_year", "imp_year", "pair_id2"),
                      robust = TRUE,
                      cluster = 'pair_id')

summary(fit5)

# Column 6 (PPML) ----------------------------------------------------------
# Allowing for potential non-linear and phasing-in effects of RTAs

fit6 <- gravity_ppml3(y = "trade",
                      x = c("RTA", "RTA_LAG4", "RTA_LAG8", "RTA_LAG12"),
                      data = data[data$sum_trade > 0, ],
                      fixed_effects = c("exp_year", "imp_year", "pair_id2"),
                      robust = TRUE,
                      cluster = 'pair_id')

summary(fit6)

# Column 7 (PPML) ----------------------------------------------------------
# Addressing Globalization Effects

fit7 <- gravity_ppml3(y = "trade",
                      x = c("RTA", "RTA_LAG4", "RTA_LAG8", "RTA_LAG12",
                            "INTL_BRDR2.year1986", "INTL_BRDR2.year1990",
                            "INTL_BRDR2.year1994", "INTL_BRDR2.year1998",
                            "INTL_BRDR2.year2002"),
                      data = data[data$sum_trade > 0, ],
                      fixed_effects = c("exp_year", "imp_year", "pair_id2"),
                      robust = TRUE,
                      cluster = 'pair_id')

summary(fit7)


