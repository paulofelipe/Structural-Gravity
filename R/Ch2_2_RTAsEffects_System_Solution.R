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
packages <- c('haven', 'lfe', 'dplyr', 'nleqslv')

load_packages <- sapply(packages, check_and_install)


# Funções -----------------------------------------------------------------

# Parâmetros CES
betas <- function(beta){
  N <- length(beta)
  r <- rep(1, N)
  P <- colSums(beta * tij)
  Pi <- colSums(t(tij)/P * E/sum(Y))
  r[1:(N-1)] <- (beta - Y/sum(Y) * 1/Pi)[1:(N-1)]
  r[N] <- sum(beta) - 1
  r
}

# Termos multilaterais de resistência
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

# Preços de equilíbrio
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

prices2 <- function(x){
  N <- length(x)
  r <- vector(length = N)
  
  # Índices
  idx_P <- 1:(N/3)
  idx_PI <- 1:(N/3) + N/3
  idx_p <- 1:(N/3) + 2*N/3
  
  # Variáveis Endógenas
  P <- x[idx_P]
  PI <- x[idx_PI]
  p <- x[idx_p]
  Y <- p * Q
  E <- phi * p * Q
  
  ## IMR
  r[idx_P] <- P - colSums(tij/PI * Y/sum(Y))
  r[1] <- P[1] - 1
  
  ## OMR
  r[idx_PI] <- PI - colSums(t(tij)/P * E/sum(Y))
  
  ## Prices
  r[idx_p] <- p^(1-sigma) - Y/sum(Y) * ((beta * PI)^(-1))
  
  r
}

# Read data
data_nafta <- read_dta('Data/Chapter2Application2.dta')

data_nafta <- data_nafta %>% 
  mutate(exporter = ifelse(exporter == "DEU", "AAA", exporter),
         importer = ifelse(importer == "DEU", "AAA", importer),
         ln_DIST = log(DIST),
         INTL = ifelse(exporter != importer, 1, 0)) %>% 
  filter(year %in% seq(1986, 2006, 4)) %>% 
  group_by(importer, year) %>% 
  mutate(E = sum(trade)) %>% 
  ungroup() %>% 
  group_by(exporter, year) %>% 
  mutate(Y = sum(trade)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  mutate(E_R_BLN = ifelse(importer == "AAA", E, 0),
         E_R = max(E_R_BLN),
         exp_year = paste0(exporter, year),
         imp_year = paste0(importer, year),
         pair_id2 = ifelse(exporter == importer, "intra", pair_id)) %>% 
  ungroup()

# Sum of trade by pair - all years
data_nafta <- data_nafta %>% 
  group_by(pair_id) %>% 
  mutate(sum_trade = sum(trade)) %>% 
  ungroup()

# Step 1: Solve the baseline gravity model --------------------------------

# Stage 1: Obtain the estimates of pair fixed effects and RTAs
fit <- gravity_ppml(y = "trade", x = "RTA",
                     data = data_nafta %>% filter(sum_trade > 0),
                     fixed_effects = c("imp_year", "exp_year", "pair_id2"),
                     robust = TRUE,
                     cluster = "pair_id")

summary(fit)

# Get the fixed effects
data_nafta <- data_nafta %>% 
  left_join(fit$fixed.effects)

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
  filter(year == 1994) %>% 
  left_join(data_stage2) %>% 
  mutate(tij_bar = ifelse(is.na(tij_bar), tij_noRTA, tij_bar),
         tij_bln = ifelse(is.na(tij_bln), tij_bar*exp(fit$coefficients["RTA",1]*RTA),
                          tij_bln)) %>% 
  select(-tij_noRTA) %>% 
  mutate(ln_tij_bln = log(tij_bln))

# Estimate trade on "constrained" trade costs

fit_bln <- gravity_ppml(y = "trade",
                         x = NULL,
                         fixed_effects = c("importer", "exporter"),
                         offset = data_nafta$ln_tij_bln,
                         data = data_nafta,
                         robust = TRUE,
                         cluster = "pair_id")

data_nafta <- data_nafta %>% 
  mutate(tradehat_bln = as.vector(fit_bln$fitted.values))

# Construct baseline indexes
# Fixed effects
data_nafta <- data_nafta %>% 
  left_join(fit_bln$fixed.effects)

# Outward and Inward multilateral resistance terms
data_nafta <- data_nafta %>% 
  mutate(OMR_BLN = Y * E_R/ exp(fe_exporter),
         IMR_BLN = E / (exp(fe_importer) * E_R))

# Output and Expenditure - Baseline
data_nafta <- data_nafta %>% 
  mutate(Y_BLN = Y,
         E_BLN = E) %>% 
  rename(fe_exp_bln = fe_exporter,
         fe_imp_bln = fe_importer)

# Step 2: Define a couterfactual scenario -------------------------------------

# RTA = 0 if pair in NAFTA

NAFTA <- c("CAN", "MEX", "USA")
data_nafta <- data_nafta %>% 
  mutate(RTA_NO_NAFTA = ifelse(exporter %in% NAFTA & importer %in% NAFTA,
                               0,
                               RTA),
         tij_cfl = tij_bar * exp(fit$coefficients["RTA",1]*RTA_NO_NAFTA),
         ln_tij_cfl = log(tij_cfl))

# Select variables for the system
data_nafta <- data_nafta %>% 
  select(exporter, importer, trade, t_bsln = tij_bln, t_crfl = tij_cfl) %>% 
  arrange(exporter, importer)


# Computar Y, E, Q, phi -------------------------------------------------------

Y <- data_nafta %>% 
  group_by(exporter) %>% 
  summarise(Y = sum(trade)) %>% 
  .[["Y"]]

E <- data_nafta %>% 
  group_by(importer) %>% 
  summarise(E = sum(trade)) %>% 
  .[["E"]]

Q <- Y

phi <- E/Q

# Custos de Comércio - Baseline -----------------------------------------------
tc_bsln <- matrix(data_nafta$t_bsln, nrow = sqrt(nrow(data_nafta)), byrow = TRUE)
tij <- tc_bsln
# Calibrar betas --------------------------------------------------------------
sigma <- 7
N <- nrow(tij)


mrts_bln <- nleqslv(rep(0.9, N * 2),
                    mrts, control = list(maxit = 2000))$x


# Índices Baseline ------------------------------------------------------------

# IMR - Baseline
P_BSLN <- mrts_bln[1:N]^(1/(1-sigma))
# OMR - Baseline
PI_BSLN <- mrts_bln[(N + 1):(2 * N)]^(1/(1-sigma))

# Renda Real - RGDP
RGDP_BSLN <- Y/P_BSLN

# Betas (beta^(1 -sigma)) -----------------------------------------------------
beta <- (Y/sum(Y)) * 1/(PI_BSLN^(1-sigma))

# Checar se os betas formam a solução
nleqslv(rep(1.1, N),
        prices)
max(abs(1 - nleqslv(rep(1.1, N), prices)$x)) < 1e-8


# Índices - Condicional -------------------------------------------------------
tij <- matrix(data_nafta$t_crfl, nrow = sqrt(nrow(data_nafta)), byrow = TRUE)

# Valores Iniciais
mrts_cdl <- nleqslv(rep(0.9, N * 2),
                       mrts, control = list(maxit = 2000))$x

# IMR - Condicional
P_CDL <- mrts_cdl[1:N]^(1/(1-sigma))

# OMR - Condicional
PI_CDL <- mrts_cdl[(N + 1):(2*N)]^(1/(1-sigma))

# Solução: Contra Factual ----------------------------------------------------
p <- nleqslv(rep(1, N),
             prices)$x

# Índices Contrafactual ---------------------------------------------------

# IMR - Contrafactual
P_CRFL <- colSums(beta * p^(1-sigma) * tij)^(1/(1-sigma))

# OMR - Contractual
PI_CRFL <- colSums(t(tij)/P_CRFL^(1-sigma) * E/sum(Y))^(1/(1-sigma))

# Renda Real - RGDP - Contrafactual
RGDP_CRFL<- p*Q/P_CRFL

# Combinar os resultados com os data_nafta originais

resultados <- data.frame(country = sort(unique(data_nafta$exporter)),
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

data_nafta <- data_nafta %>% 
  left_join(resultados, by = c("exporter" = "country")) %>% 
  left_join(resultados, by = c("importer" = "country"),
            suffix = c("_exporter", "_importer"))

data_nafta <- data_nafta %>% 
  mutate(Y_WLD_BSLN = sum(Y_BSLN_exporter * (1- (exporter != importer))),
         Y_WLD_CRFL = sum(Y_CRFL_exporter * (1- (exporter != importer))),
         X_BSLN = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_bsln/(P_BSLN_importer * PI_BSLN_exporter)^(1-sigma),
         X_PART = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_crfl/(P_BSLN_importer * PI_BSLN_exporter)^(1-sigma),
         X_CDL = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_crfl/(P_CDL_importer * PI_CDL_exporter)^(1-sigma),
         X_FULL = (Y_CRFL_exporter * E_CRFL_importer/Y_WLD_CRFL)  * 
           t_crfl/(P_CRFL_importer * PI_CRFL_exporter)^(1-sigma))

indexes_exp <- data_nafta %>% 
  group_by(exporter) %>% 
  mutate(Xi_BSLN = sum(X_BSLN * (exporter != importer)),
         Xi_PART = sum(X_PART * (exporter != importer)),
         Xi_CDL = sum(X_CDL * (exporter != importer)),
         Xi_FULL = sum(X_FULL * (exporter != importer))) %>% 
  select(exporter, starts_with("Xi_"), Y_BSLN = Y_BSLN_exporter,
         Y_FULL = Y_CRFL_exporter,
         PI_BSLN = PI_BSLN_exporter, PI_CDL = PI_CDL_exporter,
         PI_FULL = PI_CRFL_exporter,
         p_BSLN = p_BSLN_exporter,
         p_FULL = p_CRFL_exporter,) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(exporter = ifelse(exporter == "AAA", "DEU", exporter)) %>% 
  arrange(exporter) %>% 
  mutate(CHNG_X_PART = (Xi_BSLN/Xi_PART - 1) * 100,
         CHNG_X_CDL = (Xi_BSLN/Xi_CDL - 1) * 100,
         CHNG_X_FULL = (Xi_BSLN/Xi_FULL - 1) * 100,
         CHNG_OMR_FULL = (PI_BSLN/PI_FULL -  1) * 100,
         CHNG_p_FULL = (p_BSLN/p_FULL - 1) * 100) %>% 
  select(exporter, starts_with("CHNG"), starts_with("Y_"))

indexes_imp <- data_nafta %>% 
  select(importer, P_BSLN = P_BSLN_importer,
         P_CDL = P_CDL_importer,
         P_FULL = P_CRFL_importer) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(importer = ifelse(importer == "AAA", "DEU", importer)) %>% 
  arrange(importer) %>%
  mutate(CHNG_IMR_CDL = (P_BSLN/P_CDL - 1) * 100,
         CHNG_IMR_FULL = (P_BSLN/P_FULL - 1) * 100)
  
left_join(indexes_exp, indexes_imp,
          by = c("exporter" = "importer")) %>% 
  mutate(rGDP_BSLN = Y_BSLN/P_BSLN,
         rGDP_FULL = Y_FULL/P_FULL,
         CHNG_rGDP_FULL = (rGDP_BSLN/rGDP_FULL - 1) * 100) %>% 
  select(country = exporter,
         CHNG_X_PART, CHNG_X_CDL, CHNG_X_FULL, CHNG_rGDP_FULL,
         CHNG_p_FULL, CHNG_OMR_FULL, CHNG_IMR_FULL) %>% 
  mutate_if(is.numeric, round, digits = 2) %>% 
  print(n = 69)
