# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 2: General Equilibrium Trade Policy Analysis with Structural Gravity
# Application 2: IMPACT OF REGIONAL TRADE AGREEMENTS
# Date: 2017/06/07

# Function to load (and install packages if necessary) --------------------

rm(list = ls())
check_and_install <- function(package_name){
  if(!require(package_name, character.only = TRUE)){
    install.packages(package_name)
    library(package_name, character.only = TRUE)
  }
}

# load gravity_ppml function
source('R/aux_functions.R')

# Load packages -----------------------------------------------------------
packages <- c('haven', 'lfe', 'dplyr')

load_packages <- sapply(packages, check_and_install)

# Read data
data_nafta <- read_dta('Data/Chapter2Application2.dta')

data_nafta <- data_nafta %>% 
  mutate(exporter = ifelse(exporter == "DEU", "AAA", exporter),
         importer = ifelse(importer == "DEU", "AAA", importer),
         ln_DIST = log(DIST),
         INTL = ifelse(exporter != importer, 1, 0)) %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  group_by(importer) %>% 
  mutate(E = sum(trade)) %>% 
  group_by(exporter) %>% 
  mutate(Y = sum(trade)) %>% 
  ungroup() %>% 
  mutate(E_R_BLN = ifelse(importer == "AAA", E, 0),
         E_R = max(E_R_BLN),
         exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year),
         pair_id2 = ifelse(exporter == importer, "intra", pair_id))

# Sum of trade by pair - all years
data_nafta <- data_nafta %>% 
  group_by(pair_id) %>% 
  mutate(sum_trade = sum(trade)) %>% 
  ungroup()

# Step 1: Solve the baseline gravity model --------------------------------
fit <- gravity_ppml3(y = "trade", x = "RTA",
                     data = data_nafta %>% filter(sum_trade > 0),
                     fixed_effects = c("imp_year", "exp_year", "pair_id2"),
                     robust = TRUE,
                     cluster = "pair_id")

summary(fit)

# Get the fixed effects
fe <- getfe(fit)
fe_imp <- getfe(fit) %>% filter(fe == "imp_year") %>% 
  select(imp_year = idx, fe_imp_year = effect)
fe_exp <- getfe(fit) %>% filter(fe == "exp_year") %>% 
  select(exp_year = idx, fe_exp_year = effect)
fe_pair <- getfe(fit) %>% filter(fe == "pair_id2") %>% 
  select(pair_id2 = idx, fe_pair_id2 = effect)

# Join data 
data_nafta <- data_nafta %>% 
  left_join(fe_imp) %>% 
  left_join(fe_exp) %>% 
  left_join(fe_pair)

# Generate trade costs
data_nafta <- data_nafta %>% 
  mutate(tij_bar = exp(fe_pair_id2),
         tij_bln = exp(fe_pair_id2 + fit$coefficients["RTA",1]*RTA))

# Stage 2: predict trade costs for the pairs with zero trade flows
data_stage2 <- data_nafta %>% 
  filter(year == 1994, exporter != importer) %>% 
  mutate(tij = exp(fe_pair_id2))

fit_stage2 <- glm(tij ~ ln_DIST + CNTG + LANG + CLNY + importer + exporter - 1,
                  data = data_stage2, family = quasipoisson())

data_stage2$tij_noRTA <- predict(fit_stage2, newdata = data_stage2,
                                 type = "response")

data_stage2 <- data_stage2 %>% 
  select(exporter, importer, tij_noRTA)

data_nafta <- data_nafta %>% 
  left_join(data_stage2) %>% 
  mutate(tij_bar = ifelse(is.na(tij_bar), tij_noRTA, tij_bar),
         tij_bln = ifelse(is.na(tij_bln), tij_bar*exp(fit$coefficients["RTA",1]*RTA),
                          tij_bln)) %>% 
  select(-tij_noRTA) %>% 
  mutate(ln_tij_bln = log(tij_bln))


  