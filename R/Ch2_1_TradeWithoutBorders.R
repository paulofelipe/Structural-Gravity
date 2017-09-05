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

# Plot theme --------------------------------------------------------------

theme_gravity <- function(base_size = 12){
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(color = "grey85"),
        panel.grid.major = element_line(color = "grey85"),
        title = element_text(face = 'bold'),
        legend.position = "bottom",
        legend.key = element_rect(fill = "white"),
        axis.text = element_text(color = "black"),
        axis.ticks = element_blank())
}

# Load packages -----------------------------------------------------------
packages <- c('haven', 'lfe', 'dplyr', 'ggplot2', 'hrbrthemes')

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
data_borders <- data_borders %>% 
  mutate(tradehat_BLN = as.vector(fit$fitted.values)) %>% 
  group_by(exporter) %>% 
  mutate(Xi_BLN = sum(tradehat_BLN * (exporter != importer))) %>% 
  ungroup()

# Fixed effects
data_borders <- data_borders %>% 
  left_join(fit$fixed.effects)

# Trade costs - Baseline
data_borders <- data_borders %>% 
  mutate(tij_BLN = exp(fit$coefficients["ln_DIST", 1] * ln_DIST +
                         fit$coefficients["CNTG", 1] * CNTG +
                         fit$coefficients["INTL", 1] * INTL))

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>% 
  mutate(OMR_BLN = Y * E_R/ exp(fe_exporter),
         IMR_BLN = E / (exp(fe_importer) * E_R))

# Output and Expenditure - Baseline
data_borders <- data_borders %>% 
  mutate(Y_BLN = Y,
         E_BLN = E,
         Y_WLD_BLN = Y_WLD) %>% 
  rename(fe_exp_bln = fe_exporter,
         fe_imp_bln = fe_importer)

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

# Fixed effects
data_borders <- data_borders %>% 
  left_join(fit_cdl$fixed.effects)

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>% 
  mutate(OMR_CDL = Y * E_R/ (exp(fe_exporter)),
         IMR_CDL = E / (exp(fe_importer) * E_R))

# Predicted trade - Conditional
data_borders <- data_borders %>% 
  mutate(tradehat_CDL = as.vector(fit_cdl$fitted.values)) %>% 
  group_by(exporter) %>% 
  mutate(Xi_CDL = sum(tradehat_CDL * (exporter != importer))) %>%
  ungroup() %>% 
  rename(fe_exp_cdl = fe_exporter,
         fe_imp_cdl = fe_importer)

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
  mutate(change_p_i = ((exp(fe_exp_cdl) / E_R) / (exp(fe_exp_bln) / E_R))^(1 /(1 - sigma)))

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
         OMR_CFL_0 = OMR_CDL,
         IMR_CFL_0 = IMR_CDL,
         E_R_CFL_0 = E_R,
         fe_exp_cfl_0 = fe_exp_cdl,
         fe_imp_cfl_0 = fe_imp_cdl)

# Update trade flows - Counterfactual
data_borders <- data_borders %>% 
  mutate(trade_CFL = tradehat_CDL * change_p_i * change_p_j)

# Start loop
data_borders <- data_borders %>% 
  mutate(change_IMR_FULL_0 = 1,
         change_OMR_FULL_0 = 1,
         change_p_i_0 = change_p_i,
         change_p_j_0 = change_p_j,
         tradehat_0 = tradehat_CDL)
max_dif <- 1
sd_dif <- 1
change_price_i_old <- 0
while(sd_dif > 1e-10 | max_dif > 1e-10){
  
  data_borders <- data_borders %>% 
    mutate(trade_1 = tradehat_0 * change_p_i_0 * change_p_j_0 /(change_OMR_FULL_0 * change_IMR_FULL_0))
  
  fit_cfl <- gravity_ppml3(y = "trade_1", x = NULL,
                           offset = log(data_borders$tij_CFL),
                           fixed_effects = c("importer", "exporter"),
                           data = data_borders,
                           robust = TRUE)
  
  # Fixed effects
  data_borders <- data_borders %>% 
    left_join(fit_cfl$fixed.effects)
  
  data_borders <- data_borders %>% 
    mutate(tradehat_1 = fit_cfl$fitted.values) %>% 
    group_by(exporter) %>% 
    mutate(Y_CFL_1 = sum(tradehat_1)) %>% 
    ungroup() %>% 
    mutate(E_CFL_1 = ifelse(importer == exporter, phi * Y_CFL_1, 0)) %>% 
    group_by(importer) %>% 
    mutate(E_CFL_1 = max(E_CFL_1)) %>% 
    ungroup() %>% 
    mutate(E_R_CFL_1 = ifelse(importer == "AAA", E_CFL_1, 0),
           E_R_CFL_1 = max(E_R_CFL_1))
  
  # Changes in prices for exporters
  data_borders <- data_borders %>% 
    mutate(change_p_i_1 = ((exp(fe_exporter) / E_R_CFL_1) / (exp(fe_exp_cfl_0) / E_R_CFL_0))^(1 /(1 - sigma)))
  
  # Changes in prices for importers
  data_borders <- data_borders %>% 
    group_by(importer) %>% 
    mutate(change_p_j_1 = ifelse(importer == exporter, change_p_i_1, 0),
           change_p_j_1 = max(change_p_j_1)) %>% 
    ungroup()
  
  # Outward and Inward multilateral resistance terms
  data_borders <- data_borders %>% 
    mutate(OMR_CFL_1 = Y_CFL_1 * E_R_CFL_1/ (exp(fe_exporter)),
           IMR_CFL_1 = E_CFL_1 / (exp(fe_importer) * E_R_CFL_1))
  
  
  
  # Update dif
  print(summary(data_borders$change_p_i_0 - change_price_i_old))
  max_dif <- abs(max(data_borders$change_p_i_0 - change_price_i_old))
  sd_dif <- sd(data_borders$change_p_i_0 - change_price_i_old)
  change_price_i_old <- data_borders$change_p_i_0
  
  # Changes in OMR and IMR
  data_borders <- data_borders %>% 
    mutate(change_IMR_FULL_1 = IMR_CFL_1/IMR_CFL_0,
           change_OMR_FULL_1 = OMR_CFL_1/OMR_CFL_0,
           IMR_CFL_0 = IMR_CFL_1,
           OMR_CFL_0 = OMR_CFL_1,
           change_IMR_FULL_0 = change_IMR_FULL_1,
           change_OMR_FULL_0 = change_OMR_FULL_1,
           change_p_i_0 = change_p_i_1,
           change_p_j_0 = change_p_j_1,
           tradehat_0 = tradehat_1,
           E_R_CFL_0 = E_R_CFL_1,
           fe_exp_cfl_0 = fe_exporter,
           fe_imp_cfl_0 = fe_importer) %>% 
    select(-fe_exporter, -fe_importer)

}

### CONSERTAR CHANGE E_FULL
data_borders <- data_borders %>% 
  mutate(change_p_i_full = ((exp(fe_exp_cfl_0) / E_R_CFL_0) / (exp(fe_exp_bln) / E_R))^(1 /(1 - sigma)),
         change_p_j_full = change_p_i_full * (exporter == importer)) %>% 
  group_by(importer) %>% 
  mutate(change_p_j_full = max(change_p_j_full)) %>% 
  ungroup() %>% 
  mutate(Y_FULL = change_p_i_full * Y_BLN,
         E_FULL = change_p_j_full * E_BLN * (exporter == importer)) %>% 
  group_by(importer) %>% 
  mutate(E_FULL = max(E_FULL, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(E_FULL_R = E_FULL * (importer == "AAA"),
         E_FULL_R = max(E_FULL_R),
         OMR_FULL = Y_FULL * E_R_CFL_0 / exp(fe_exp_cfl_0),
         IMR_FULL = E_FULL / (exp(fe_imp_cfl_0) * E_FULL_R),
         X_FULL = (Y_FULL * E_FULL * tij_CFL)/(IMR_FULL * OMR_FULL)) %>% 
  group_by(exporter) %>% 
  mutate(Xi_FULL = sum(X_FULL * (importer != exporter))) %>% 
  ungroup()

# Collect indexes of interest ------------------------------------------------
exporter_indexes <- data_borders %>%
  filter(importer == exporter) %>% 
  select(exporter, starts_with("OMR"), change_p_i_full, starts_with("Xi_"),
         Y_BLN, Y_FULL) %>% 
  mutate(change_price_FULL = (change_p_i_full - 1) * 100,
         change_OMR_CDL = (OMR_CDL^(1/(1 - sigma)) / OMR_BLN^(1/(1 - sigma)) - 1) * 100,
         change_OMR_FULL = (OMR_FULL^(1/(1 - sigma)) / OMR_BLN^(1/(1 - sigma)) - 1) * 100,
         change_Xi_CDL = (Xi_CDL / Xi_BLN  - 1) * 100,
         change_Xi_Full = (Xi_FULL / Xi_BLN - 1) * 100)

importer_indexes <- data_borders %>% 
  filter(importer == exporter) %>% 
  select(importer, starts_with("IMR"), change_p_i_full, starts_with("Y_")) %>% 
  mutate(change_price_FULL = (change_p_i_full - 1) * 100,
         change_IMR_FULL = -(IMR_FULL^(1/(1 - sigma)) / IMR_BLN^(1/(1 - sigma)) - 1) * 100,
         rGDP = ((Y_FULL/IMR_FULL^(1/(1 - sigma)))/(Y_BLN/IMR_BLN^(1/(1 - sigma))) - 1) * 100)


# Visualizations --------------------------------------------------------------

# Figure 6 - page 100
ggplot(exporter_indexes %>% filter(exporter != "HKG"),
       aes(x = log(Y_BLN))) +
  geom_point(aes(y = change_Xi_Full, color = "a"), size = 3.5, alpha = 0.7) + 
  geom_point(aes(y = change_Xi_CDL, color = "b"), size = 3.5, alpha = 0.7) + 
  labs(title = "Effects of abolishing international borders on exports",
       y = "Percent change of exports",
       x = "Log value of output") +
  scale_color_manual("",values = c("#280052", "#4A9BE6"),
                      labels = c("Full Endowment General Equilibrium",
                                 "Conditional General Equilibrium")) +
  theme_gravity()

# Figure 7 - page 100
ggplot(importer_indexes %>% filter(importer != "HKG"),
       aes(x = log(Y_BLN))) +
  geom_point(aes(y = rGDP, color = "a", shape = "a"), size = 3.5, alpha = 0.7) +
  geom_point(aes(y = change_price_FULL, color = "b", shape = "b"), size = 3.5, alpha = 0.7) +
  geom_point(aes(y = change_IMR_FULL, color = "c", shape = "c"), size = 3.5, alpha = 0.7) +
  labs(title = "Effects of abolishing international borders on real GDP",
       y = "Percent change",
       x = "Log value of output") +
  scale_color_manual("",values = c("#280052", "#4A9BE6", "#000878"),
                     labels = c("Real GDP",
                                "Factory-gate price",
                                "-(inward multilateral resistances)")) +
  scale_shape_manual("", values = 15:17,
                     labels = c("Real GDP",
                                "Factory-gate price",
                                "-(inward multilateral resistances)")) +
  theme_gravity()


teste <- data.frame(letra = letters[1:8], y = rnorm(8))
ggplot(teste, aes(x = letra, y = y, fill = letra)) +
  geom_col() +
  scale_fill_deaex() +
  theme_gravity()
