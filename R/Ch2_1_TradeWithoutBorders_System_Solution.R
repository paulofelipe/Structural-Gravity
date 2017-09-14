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
packages <- c('haven', 'lfe', 'dplyr', 'nleqslv', 'ggplot2')

load_packages <- sapply(packages, check_and_install)

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


# Functions -----------------------------------------------------------------

# CES Parameters
betas <- function(beta){
  N <- length(beta)
  r <- rep(1, N)
  P <- colSums(beta * tij)
  Pi <- colSums(t(tij)/P * E/sum(Y))
  r[1:(N-1)] <- (beta - Y/sum(Y) * 1/Pi)[1:(N-1)]
  r[N] <- sum(beta) - 1
  r
}

# Multilateral Resistance Terms
mrts <- function(m){
  N <- length(m)
  r <- vector(length = N, mode = "numeric")
  
  ## IMR
  P <- colSums(tij/m[(N/2 + 1):N] * Y/sum(Y))
  
  ## OMR
  PI <-  colSums(t(tij)/m[1:(N/2)] * E/sum(Y))
  
  ## Solution
  r[1] <- m[1] - 1
  r[2:(N/2)] <- m[2:(N/2)] - P[2:(N/2)]
  r[(N/2 + 1):N] <- m[(N/2 + 1):N] - PI[1:(N/2)]
  r
  
}

# Equilibrium Prices
prices <- function(p){
  N <- length(p)
  r <- vector(length = N)
  
  Y <- p * Q
  E <- phi * p * Q
  
  ## IMR
  P <- colSums(beta * p^(1-sigma) * tij)
  
  ## OMR
  PI <-  colSums(t(tij)/P * E/sum(Y))
  
  ## Prices
  r[1] <- P[1] - 1  
  r[2:N] <- (p^(1-sigma) - Y/sum(Y) * ((beta * PI)^(-1)))[2:N]
  
  r
}

# Load Data -------------------------------------------------------
dados <- read_dta('Data/Chapter2Application1.dta')

# Filter and generate new variables
dados <- dados %>% 
  filter(year == 2006) %>% 
  mutate(ln_DIST = log(DIST),
         INTL_BRDR = ifelse(exporter != importer, 1, 0)) %>% 
  group_by(exporter) %>% 
  mutate(Y = sum(trade)) %>% 
  ungroup() %>% 
  group_by(importer) %>% 
  mutate(E = sum(trade)) %>% 
  ungroup()

# Define Germany as reference country
dados <- dados %>%
  mutate(exporter = ifelse(exporter == "DEU", "AAA", exporter),
         importer = ifelse(importer == "DEU", "AAA", importer))

# Estimate the gravitational equation
fit <- gravity_ppml3(y = "trade", x = c("ln_DIST", "CNTG", "INTL_BRDR"),
                     fixed_effects = c("importer", "exporter"),
                     data = dados)

summary(fit)

# Compute trade costs - baseline and counterfactual
dados <- dados %>% 
  mutate(t_bsln = exp(fit$coefficients["ln_DIST",] * ln_DIST +
                        fit$coefficients["CNTG",] * CNTG +
                        fit$coefficients["INTL_BRDR",] * INTL_BRDR),
         t_crfl = exp(fit$coefficients["ln_DIST",] * ln_DIST +
                        fit$coefficients["CNTG",] * CNTG +
                        fit$coefficients["INTL_BRDR",] * INTL_BRDR * 0)) %>% 
  select(exporter, importer, trade, t_bsln, t_crfl) %>% 
  arrange(exporter, importer)


# Compute Y, E, Q, phi -------------------------------------------------------

Y <- dados %>% 
  group_by(exporter) %>% 
  summarise(Y = sum(trade)) %>% 
  .[["Y"]]

E <- dados %>% 
  group_by(importer) %>% 
  summarise(E = sum(trade)) %>% 
  .[["E"]]

Q <- Y

phi <- E/Q

# Trade costs matrix - Baseline -----------------------------------------------
tc_bsln <- matrix(dados$t_bsln, nrow = sqrt(nrow(dados)), byrow = TRUE)
tij <- tc_bsln

# Compute MRTs - Baseline -----------------------------------------------------
sigma <- 7
N <- nrow(tij)

start_mrts <- c(1, rep(0.5, N-1),
                rep(0.01, N))

mrts_bln <- nleqslv(start_mrts,
                    mrts, control = list(maxit = 2000))$x

# Indexes - Baseline ----------------------------------------------------------

# IMR - Baseline
P_BSLN <- mrts_bln[1:N]^(1/(1-sigma))
# OMR - Baseline
PI_BSLN <- mrts_bln[(N + 1):(2 * N)]^(1/(1-sigma))

# RGDP
RGDP_BSLN <- Y/P_BSLN

# Betas (beta^(1 -sigma)) -----------------------------------------------------
beta <- (Y/sum(Y)) * 1/(PI_BSLN^(1-sigma))

# Check if beta is the solution
nleqslv(rep(1.1, N),
        prices)
max(abs(1 - nleqslv(rep(1.1, N), prices)$x)) < 1e-8


# Trade costs matrix - Counterfactual -----------------------------------------
tij <- matrix(dados$t_crfl, nrow = sqrt(nrow(dados)), byrow = TRUE)

# Valores Iniciais
mrts_cdl <- nleqslv(rep(0.9, N * 2),
                       mrts, control = list(maxit = 2000))$x

# IMR - Condicional
P_CDL <- mrts_cdl[1:N]^(1/(1-sigma))

# OMR - Condicional
PI_CDL <- mrts_cdl[(N + 1):(2*N)]^(1/(1-sigma))

# Solution: Counterfactual ----------------------------------------------------
p <- nleqslv(rep(1, N),
             prices)$x

# Indexes - Counterfactual ---------------------------------------------------

# IMR - Counterfactual
P_CRFL <- colSums(beta * p^(1-sigma) * tij)^(1/(1-sigma))

# OMR - Counterfactual
PI_CRFL <- colSums(t(tij)/P_CRFL^(1-sigma) * E/sum(Y))^(1/(1-sigma))

# RGDP - Counterfactual
RGDP_CRFL<- p*Q/P_CRFL

# Merge results 
# Indexes by country

resultados <- data.frame(country = sort(unique(dados$exporter)),
                         p_BSLN = 1,
                         p_CRFL = p,
                         P_BSLN = P_BSLN,
                         P_CDL = P_CDL,
                         P_CRFL = P_CRFL, 
                         PI_BSLN = PI_BSLN,
                         PI_CDL = PI_CDL,
                         PI_CRFL = PI_CRFL,
                         Y_BSLN = Q,
                         Y_CRFL = p * Q,
                         E_BSLN = phi * Q,
                         E_CRFL = phi * p * Q,
                         RGDP_BSLN = RGDP_BSLN,
                         RGDP_CRFL = RGDP_CRFL,
                         stringsAsFactors = FALSE)

# Combine indexes for exporter and importer

dados <- dados %>% 
  left_join(resultados, by = c("exporter" = "country")) %>% 
  left_join(resultados, by = c("importer" = "country"),
            suffix = c("_exporter", "_importer"))

# Compute trade - baseline, conditional e counterfactual
dados <- dados %>% 
  mutate(Y_WLD_BSLN = sum(Y_BSLN_exporter * (1- (exporter != importer))),
         Y_WLD_CRFL = sum(Y_CRFL_exporter * (1- (exporter != importer))),
         X_BSLN = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_bsln/(P_BSLN_importer * PI_BSLN_exporter)^(1-sigma),
         X_CDL = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_crfl/(P_CDL_importer * PI_CDL_exporter)^(1-sigma),
         X_CRFL = (Y_CRFL_exporter * E_CRFL_importer/Y_WLD_CRFL)  * 
           t_crfl/(P_CRFL_importer * PI_CRFL_exporter)^(1-sigma))

# Summarise export results
results_exp <- dados %>% 
  group_by(exporter) %>% 
  summarise(Y_BSLN = mean(Y_BSLN_exporter),
            X_BSLN = sum(X_BSLN * (exporter != importer)),
            X_CDL = sum(X_CDL * (exporter != importer)),
            X_CRFL = sum(X_CRFL * (exporter != importer))) %>% 
  mutate(CHNG_X_CDL = (X_CDL/X_BSLN - 1) * 100,
         CHNG_X_FULL = (X_CRFL/X_BSLN - 1) * 100) %>% 
  filter(exporter != "HKG")

# Plots -----------------------------------------------------------------------

# Figure 6 - page 100
ggplot(results_exp,
       aes(x = log(Y_BSLN))) +
  geom_point(aes(y = CHNG_X_FULL, color = "a"), size = 3.5, alpha = 0.7) + 
  geom_point(aes(y = CHNG_X_CDL, color = "b"), size = 3.5, alpha = 0.7) + 
  labs(title = "Effects of abolishing international borders on exports",
       y = "Percent change of exports",
       x = "Log value of output") +
  scale_color_manual("",values = c("#280052", "#4A9BE6"),
                     labels = c("Full Endowment General Equilibrium",
                                "Conditional General Equilibrium")) +
  theme_gravity()

# Figure 7 - page 100
ggplot(resultados %>% filter(country != "HKG"),
       aes(x = log(Y_BSLN))) +
  geom_point(aes(y = (RGDP_CRFL/RGDP_BSLN - 1) * 100, color = "a", shape = "a"), size = 3.5, alpha = 0.7) +
  geom_point(aes(y = (p_CRFL - 1) * 100, color = "b", shape = "b"), size = 3.5, alpha = 0.7) +
  geom_point(aes(y = -(P_CRFL/P_BSLN - 1) * 100, color = "c", shape = "c"), size = 3.5, alpha = 0.7) +
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
