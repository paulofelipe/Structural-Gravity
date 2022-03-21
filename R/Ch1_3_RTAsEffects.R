# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 3: Reginal trade agreements effects
# Code to replicate the results from Table 3 (page 51)
# Last update: 2021/02/07

# Load packages -----------------------------------------------------------
library(haven)
library(dplyr)
library(estimatr)
library(lmtest)
library(fixest)

# Read and Prepare Data --------------------------------------------------

data <- read_dta("Data/Chapter1Application3.dta")

# Create new variables
data <- data %>%
  filter(year %in% seq(1986, 2006, 4)) %>%
  mutate(
    ln_trade = log(trade),
    ln_DIST = log(DIST),
    exp_year = paste0(exporter, year),
    imp_year = paste0(importer, year),
    pair_id2 = ifelse(exporter == importer, "intra", pair_id),
    INTL_BRDR = ifelse(exporter != importer, "inter", pair_id),
    INTL_BRDR2 = ifelse(exporter != importer, 1, 0),
    year = as.factor(year)
  )

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

fit1 <- feols(
  ln_trade ~ ln_DIST + CNTG + LANG + CLNY + RTA | exp_year + imp_year,
  data = data %>% 
    filter(trade > 0, importer != exporter)
)

summary(fit1, cluster = ~ pair_id)

# Column 2 (PPML) ----------------------------------------------------------
# Traditional Estimates

fit2 <- fepois(
  trade ~ ln_DIST + CNTG + LANG + CLNY + RTA | exp_year + imp_year,
  data = data %>% 
    filter(importer != exporter)
)

summary(fit2, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# Column 3 (PPML) ----------------------------------------------------------
# Traditional Estimates with Intra-national Trade

fit3 <- fepois(
  trade ~ ln_DIST + CNTG + LANG + CLNY + RTA | exp_year + imp_year + INTL_BRDR,
  data = data
)

summary(fit3, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# Column 4 (PPML) ----------------------------------------------------------
# Addressing Endogeneity of RTAs

fit4 <- fepois(
  trade ~ RTA | exp_year + imp_year + pair_id2,
  data = data %>%
    filter(sum_trade > 0)
)

summary(fit4, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# Column 5 (PPML) ----------------------------------------------------------
# Testing potential reverse causality between Trade and RTA

fit5 <- fepois(
  trade ~ RTA + RTA_LEAD4 | exp_year + imp_year + pair_id2,
  data = data %>%
    filter(sum_trade > 0)
)

summary(fit5, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# Column 6 (PPML) ----------------------------------------------------------
# Allowing for potential non-linear and phasing-in effects of RTAs

fit6 <- fepois(
  trade ~ RTA + RTA_LAG4 + RTA_LAG8 + RTA_LAG12 | exp_year + imp_year + pair_id2,
  data = data %>%
    filter(sum_trade > 0)
)

summary(fit6, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# Column 7 (PPML) ----------------------------------------------------------
# Addressing Globalization Effects

fit7 <- fepois(
  trade ~ RTA + RTA_LAG4 + RTA_LAG8 + RTA_LAG12 + 
  INTL_BRDR2.year1986 + INTL_BRDR2.year1990 + INTL_BRDR2.year1994 + 
  INTL_BRDR2.year1998 + INTL_BRDR2.year2002 | exp_year + imp_year + pair_id2,
  data = data %>%
    filter(sum_trade > 0)
)

summary(fit7, se = "cluster", cluster = ~pair_id, dof(fixef.K = "none"))
