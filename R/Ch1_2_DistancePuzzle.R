# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 1: Partial Equilibrium Trade Policy Analysis with Structural Gravity
# Application 2: The "distance puzzle" resolved
# Code to replicate the results from Table 2 (page 47)
# Last update: 2021/02/07

# Load packages -----------------------------------------------------------
library(haven)
library(dplyr)
library(estimatr)
library(lmtest)
library(fixest)

# Read and Prepare Data --------------------------------------------------

data <- read_dta("Data/Chapter1Application2.dta")

# Create new variables
data <- data %>%
  filter(year %in% seq(1986, 2006, 4)) %>%
  mutate(
    exp_year = paste0(exporter, year),
    imp_year = paste0(importer, year),
    year_f = as.factor(year)
  )

dist_year <- model.matrix(
  ~ -1 + log(DIST):year_f,
  data = data
)

dist_year <- data.frame(dist_year)

colnames(dist_year) <- c(
  "logDIST1986", "logDIST1990",
  "logDIST1994", "logDIST1998",
  "logDIST2002", "logDIST2006"
)

data <- bind_cols(data, dist_year)
# OLS with exporter/time and importer/time fixed effects ---------------------

fit <- feols(
  log(trade) ~ logDIST1986 + logDIST1990 + logDIST1994 +
    logDIST1998 + logDIST2002 + logDIST2006 +
    CNTG + LANG + CLNY |
    exp_year + imp_year,
  data = data %>%
    filter(trade > 0, importer != exporter)
)

summary(fit, cluster = ~pair_id)

# PPML with exporter/time and importer/time fixed effects --------------------

fit2 <- fepois(
  trade ~ logDIST1986 + logDIST1990 + logDIST1994 +
    logDIST1998 + logDIST2002 + logDIST2006 +
    CNTG + LANG + CLNY |
    exp_year + imp_year,
  data = data %>%
    filter(importer != exporter)
)

summary(fit2, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# PPML with intra trade ------------------------------------------------------

# create intra trade distance variable
data <- data %>%
  mutate(
    SMCTRY = ifelse(importer == exporter, 1, 0),
    logDIST_INTRA = log(DIST) * SMCTRY
  ) %>%
  mutate(
    across(
      c(contains("logDIST1"), contains("logDIST2")),
      ~ . * (1 - SMCTRY)
    )
  )

fit3 <- fepois(
  trade ~ logDIST1986 + logDIST1990 + logDIST1994 +
    logDIST1998 + logDIST2002 + logDIST2006 +
    CNTG + LANG + CLNY + logDIST_INTRA |
    exp_year + imp_year,
  data = data
)

summary(fit3, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# PPML with intra trade and Home Bias   --------------------------------------

fit4 <- fepois(
  trade ~ logDIST1986 + logDIST1990 + logDIST1994 +
    logDIST1998 + logDIST2002 + logDIST2006 +
    CNTG + LANG + CLNY + logDIST_INTRA + SMCTRY |
    exp_year + imp_year,
  data = data
)

summary(fit4, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))

# PPML with intra-national fixed effects

data <- data %>%
  mutate(intra_pair = ifelse(exporter == importer, exporter, "inter"))

fit5 <- fepois(
  trade ~ logDIST1986 + logDIST1990 + logDIST1994 +
    logDIST1998 + logDIST2002 + logDIST2006 +
    CNTG + LANG + CLNY |
    exp_year + imp_year + intra_pair,
  data = data
)

summary(fit5, se = "cluster", cluster = ~ pair_id, dof(fixef.K = "none"))