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

data$const <- 1
# Create new variables
data <- data %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  mutate(exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year),
         pair_id2 = ifelse(exporter == importer, "intra", pair_id))

fit <- gravity_ppml2(y = "trade",
                     x = c("RTA"),
                     data = data,
                     reference = "imp_yearDEU|intra",
                     fixed_effects = c("exp_year", "imp_year", "pair_id2"),
                     robust = TRUE,
                     cluster = 'pair_id', 
                     trace = TRUE)

summary.ppml(fit)
