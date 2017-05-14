# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 2: General Equilibrium Trade Policy Analysis with Structural Gravity
# Application 1: Trade without borders
# Date: 2017/05/10


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
data_borders <- read_dta('Data/Chapter2Application1.dta')

# Define Germany as reference country
data_borders <- data_borders %>% 
  mutate(exporter = ifelse(exporter == "DEU", "AAA", exporter),
         importer = ifelse(importer == "DEU", "AAA", importer),
         ln_DIST = log(DIST),
         INTL = ifelse(exporter != importer, 1, 0)) %>% 
  filter(year == 2006) %>% 
  group_by(importer) %>% 
  mutate(E = sum(trade)) %>% 
  group_by(exporter) %>% 
  mutate(Y = sum(trade)) %>% 
  ungroup() %>% 
  mutate(E_R_BLN = ifelse(importer == "AAA", E, 0),
         E_R = max(E_R_BLN))

# World output

data_borders <- data_borders %>% 
  mutate(Y_WLD = sum(Y * (1 - INTL)))


#############################################################################
########################## STEP 1 ##########################################
################# Solve the baseline gravity model ########################
###########################################################################

# Estimate Gravity Equation -----------------------------------------------
fit <- gravity_ppml3(y = "trade", x = c("ln_DIST", "CNTG", "INTL"),
                     fixed_effects = c("importer", "exporter"),
                     data = data_borders,
                     robust = TRUE)
summary(fit)

# Predicted trade in baseline model
#data_borders$tradehat_BLN <- data_borders$trade - fit$residuals

# Fixed effects
fe_exp <- getfe(fit) %>% filter(fe == "exporter") %>% 
  select(idx, fe_exp = effect)
fe_imp <- getfe(fit) %>% filter(fe == "importer") %>% 
  select(idx, fe_imp = effect)

# Exporter Fixed Effect
data_borders <- data_borders %>% 
  left_join(fe_exp, by = c("exporter" = "idx")) %>% 
  left_join(fe_imp, by = c("importer" = "idx")) %>% 
  mutate(fe_exp = exp(fe_exp),
         fe_imp = exp(fe_imp))

# Trade costs - Baseline
data_borders <- data_borders %>% 
  mutate(tij_BLN = exp(fit$coefficients["ln_DIST", 1] * ln_DIST +
                         fit$coefficients["CNTG", 1] * CNTG +
                         fit$coefficients["INTL", 1] * INTL))

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>% 
  mutate(OMR_BLN = Y * E_R/ (fe_exp),
         IMR_BLN = E / (fe_imp * E_R))

# Output and Expenditure - Baseline
data_borders <- data_borders %>% 
  mutate(Y_BLN = Y,
         E_BLN = E,
         Y_WLD_BLN = Y_WLD)

# Predicted trade - Baseline
data_borders <- data_borders %>% 
  mutate(tradehat_BLN = (Y * E * tij_BLN)/( OMR_BLN * IMR_BLN))

##############################################################################
########################## STEP 2 ###########################################
################ Define the counterfactual scenario ########################
###########################################################################

data_borders <- data_borders %>% 
  mutate(INTL_BRDR_CFL = 0,
         tij_CFL = exp(fit$coefficients["ln_DIST", 1] * ln_DIST +
                         fit$coefficients["CNTG", 1] * CNTG +
                         fit$coefficients["INTL", 1] * INTL_BRDR_CFL))

##############################################################################
########################## STEP 3 ###########################################
################# Solve the countefactual gravity model #####$##############
###########################################################################

fit_cdl <- gravity_ppml3(y = "trade", x = NULL,
                         offset = log(data_borders$tij_CFL),
                         fixed_effects = c("importer", "exporter"),
                         data = data_borders,
                         robust = TRUE)

head(getfe(fit_cdl) %>% filter(fe == "exporter"))
head(getfe(fit) %>% filter(fe == "exporter"))

#tradehat_CDL <- data_borders$trade - fit_cdl$residuals

fe_exp_cdl <- getfe(fit_cdl) %>% filter(fe == "exporter") %>% 
  select(idx, fe_exp_cdl = effect)
fe_imp_cdl <- getfe(fit_cdl) %>% filter(fe == "importer") %>% 
  select(idx, fe_imp_cdl = effect)

# Exporter and Importer Fixed Effect - Conditional
data_borders <- data_borders %>% 
  left_join(fe_exp_cdl, by = c("exporter" = "idx")) %>% 
  left_join(fe_imp_cdl, by = c("importer" = "idx")) %>% 
  mutate(fe_exp_cdl = exp(fe_exp_cdl),
         fe_imp_cdl = exp(fe_imp_cdl))

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>% 
  mutate(OMR_CDL = Y * E_R/ (fe_exp_cdl),
         IMR_CDL = E / (fe_imp_cdl * E_R))

# Predicted trade - Conditional
data_borders <- data_borders %>% 
  mutate(tradehat_CDL = (Y * E * tij_CFL)/( OMR_CDL * IMR_CDL))

# Estimated level of international trade conditional on given Y and E
# data_borders <- data_borders %>% 
#   mutate(tempXi_CDL = trade - fit_cdl$residuals) %>% 
#   group_by(exporter) %>% 
#   mutate(Xi_CDL = sum(tempXi_CDL)) %>% 
#   ungroup()

# Full endownment general equilibrium effects

# The constant of elasticity of substitution
sigma = 7

# Change in bilateral trade costs
data_borders <- data_borders %>% 
  mutate(change_tij = tij_CFL / tij_BLN)

# Deficit/Surplus parameter
data_borders <- data_borders %>% 
  mutate(phi = ifelse(importer == exporter, E/Y, 0)) %>% 
  group_by(exporter) %>% 
  mutate(phi = max(phi))

# First-order change in prices 

# Changes in prices for exporters
data_borders <- data_borders %>% 
  mutate(change_p_i = ((fe_exp_cdl / E_R) / (fe_exp / E_R))^(1 /(1 - sigma)))

# Changes in prices for importers
data_borders <- data_borders %>% 
  group_by(importer) %>% 
  mutate(change_p_j = ifelse(importer == exporter, change_p_i, 0),
         change_p_j = max(change_p_j)) %>% 
  ungroup()

# Compute change in output and expenditure
data_borders <- data_borders %>% 
  mutate(Y_CFL = Y,
         E_CFL = E,
         Y_WLD_CFL = sum(Y_CFL * (1 - INTL)),
         OMR_CFL = OMR_CDL,
         IMR_CFL = IMR_CDL,
         E_R_CFL = E_R,
         fe_exp_cfl = fe_exp_cdl,
         fe_imp_cfl = fe_imp_cdl)

# Update trade flows - Counterfactual
data_borders <- data_borders %>% 
  mutate(trade_CFL = tradehat_CDL * change_p_i * change_p_j)

# Start loop
data_borders$dif <- 1
data_borders$change_p_i_0 <- data_borders$change_p_i

while(max(abs(data_borders$dif)) > 1e-6){
  print(max(abs(data_borders$dif)))
  fit_cfl <- gravity_ppml3(y = "trade_CFL", x = NULL,
                           offset = log(data_borders$tij_CFL),
                           fixed_effects = c("importer", "exporter"),
                           data = data_borders,
                           robust = TRUE)
  
  fe_exp_cfl <- getfe(fit_cfl) %>% filter(fe == "exporter") %>% 
    select(idx, fe_exp_cfl_tmp = effect)
  fe_imp_cfl <- getfe(fit_cfl) %>% filter(fe == "importer") %>% 
    select(idx, fe_imp_cfl_tmp = effect)
  
  # Exporter and Importer Fixed Effect - Conditional
  data_borders <- data_borders %>% 
    left_join(fe_exp_cfl, by = c("exporter" = "idx")) %>% 
    left_join(fe_imp_cfl, by = c("importer" = "idx")) %>% 
    mutate(fe_exp_cfl_tmp = exp(fe_exp_cfl_tmp),
           fe_imp_cfl_tmp = exp(fe_imp_cfl_tmp))
  
  data_borders <- data_borders %>% 
    mutate(tradehat_CFL = tij_CFL * fe_exp_cfl_tmp * fe_imp_cfl_tmp) %>% 
    group_by(exporter) %>% 
    mutate(Y_CFL_tmp = sum(tradehat_CFL)) %>% 
    ungroup() %>% 
    mutate(E_CFL_tmp = ifelse(importer == exporter, phi * Y_CFL_tmp, 0)) %>% 
    group_by(importer) %>% 
    mutate(E_CFL_tmp = max(E_CFL_tmp)) %>% 
    ungroup() %>% 
    mutate(E_R_CFL_tmp = ifelse(importer == "AAA", E_CFL_tmp, 0),
           E_R_CFL_tmp = max(E_R_CFL_tmp))
  
  # Outward and Inward multilateral resistance terms
  data_borders <- data_borders %>% 
    mutate(OMR_CFL_tmp = Y_CFL_tmp * E_R_CFL_tmp/ (fe_exp_cfl_tmp),
           IMR_CFL_tmp = E_CFL_tmp / (fe_imp_cfl_tmp * E_R_CFL_tmp))
  
  # Predicted trade - Conditional
  data_borders <- data_borders %>% 
    mutate(tradehat_CFL2 = (Y_CFL_tmp * E_CFL_tmp * tij_CFL)/(OMR_CFL_tmp * IMR_CFL_tmp))
  
  # Changes in prices for exporters
  data_borders <- data_borders %>% 
    mutate(change_p_i = ((fe_exp_cfl_tmp / E_R_CFL_tmp) / (fe_exp_cfl / E_R_CFL))^(1 /(1 - sigma)))
  
  # Changes in prices for importers
  data_borders <- data_borders %>% 
    group_by(importer) %>% 
    mutate(change_p_j = ifelse(importer == exporter, change_p_i, 0),
           change_p_j = max(change_p_j)) %>% 
    ungroup()
  
  # Changes in OMR and IMR
  data_borders <- data_borders %>% 
    mutate(change_IMR_FULL = IMR_CFL_tmp/IMR_CFL,
           change_OMR_FULL = OMR_CFL_tmp/OMR_CFL)
  
    # Update trade flows - Counterfactual
  data_borders <- data_borders %>% 
    mutate(trade_CFL = tradehat_CFL2 *  change_p_i * change_p_j /
             (change_IMR_FULL * change_OMR_FULL)) 
  
  # Update dif
  data_borders$dif <- data_borders$change_p_i - data_borders$change_p_i_0
  data_borders$change_p_i_0 <- data_borders$change_p_i
  
  # update counterfactuals
  data_borders <- data_borders %>% 
    mutate(Y_CFL = Y_CFL_tmp,
           E_CFL = E_CFL_tmp,
           fe_exp_cfl = fe_exp_cfl_tmp,
           fe_imp_cfl = fe_imp_cfl_tmp,
           E_R_CFL = E_R_CFL_tmp,
           OMR_CFL = OMR_CFL_tmp,
           IMR_CFL = IMR_CFL_tmp) %>% 
    select(-fe_exp_cfl_tmp, -fe_imp_cfl_tmp)
  
}

data_borders <- data_borders %>% 
  mutate(change_p_i_full = ((fe_exp_cfl / E_R_CFL) / (fe_exp / E_R))^(1 /(1 - sigma)))

