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
  
  ## IMR
  P <- colSums(beta * p^(1-sigma) * tij)
  
  ## OMR
  PI <-  colSums(t(tij)/P * E/sum(Y))
  
  ## Prices
  r[1] <- P[1] - 1  
  r[2:N] <- (p^(1-sigma) - Y/sum(Y) * ((beta * PI)^(-1)))[2:N]
  
  r
}

# Carregar os dados -------------------------------------------------------
dados <- read_dta('Data/Chapter2Application1.dta')

# Filtra os dados e criar novas variáveis
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

# Define o país de referência como AAA
dados <- dados %>%
  mutate(exporter = ifelse(exporter == "DEU", "AAA", exporter),
         importer = ifelse(importer == "DEU", "AAA", importer))

# Estima o modelo gravitacional
fit <- gravity_ppml3(y = "trade", x = c("ln_DIST", "CNTG", "INTL_BRDR"),
                     fixed_effects = c("importer", "exporter"),
                     data = dados)

summary(fit)

# Computa os custos de comércio - baseline e contra factual
dados <- dados %>% 
  mutate(t_bsln = exp(fit$coefficients["ln_DIST",] * ln_DIST +
                        fit$coefficients["CNTG",] * CNTG +
                        fit$coefficients["INTL_BRDR",] * INTL_BRDR),
         t_crfl = exp(fit$coefficients["ln_DIST",] * ln_DIST +
                        fit$coefficients["CNTG",] * CNTG +
                        fit$coefficients["INTL_BRDR",] * INTL_BRDR * 0)) %>% 
  select(exporter, importer, trade, t_bsln, t_crfl) %>% 
  arrange(exporter, importer)


# head(dados)
# dados <- dados %>%
#   filter(exporter %in% c("CAN","BRA", "MEX", "USA") &
#            importer %in% c("CAN","BRA", "MEX", "USA"))

# Computar Y, E, Q, phi -------------------------------------------------------

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

# Custos de Comércio - Baseline -----------------------------------------------
tc_bsln <- matrix(dados$t_bsln, nrow = sqrt(nrow(dados)), byrow = TRUE)
tij <- tc_bsln
# Calibrar betas --------------------------------------------------------------
sigma <- 7
N <- nrow(tij)
# start_values <- c(0.0498231619392814, 0.0068197580875335, 0.0232387261695378, 
#                   0.00723142130490437, 0.0147143603449896, 0.00186137961886111, 
#                   0.000668271991502307, 0.0385672833591634, 0.0267505831909108, 
#                   0.0112209700412317, 0.0110160521713741, 0.151007769702764, 0.000548925527240701, 
#                   0.00465563612595528, 0.00147522715892983, 0.000932499713168157, 
#                   0.00530681746099243, 0.00147405501456217, 0.00265609872058792, 
#                   0.022845803328183, 0.00856718075450557, 0.031946412579548, 0.0226219499805236, 
#                   0.00432721871481442, 0.000470674707004898, 0.00541033778456974, 
#                   0.020589388332436, 0.0353297370328381, 0.00860810524648814, 0.00678740616236631, 
#                   0.000882518536038507, 0.0050770246890381, 0.036210511848402, 
#                   0.00175044214353065, 0.0933102306777749, 0.0011952601227573, 
#                   0.0385871687058756, 0.00102630608686662, 0.00357101350049923, 
#                   0.000153092937648639, 0.00277796720090255, 0.019475113480638, 
#                   0.0003446129913758, 0.00340127867994067, 0.000407589741189976, 
#                   0.000302833330765114, 0.0242532450491626, 5.59683502443296e-05, 
#                   0.00831951525265885, 0.0115770621520492, 0.00621020573145802, 
#                   0.000543938808466601, 0.000502897650423773, 0.00958000670242452, 
#                   0.0106265333117643, 0.00537531845173408, 0.000841903579322441, 
#                   0.00443273923842532, 0.000427263563189535, 0.00763043577964612, 
#                   0.013255883914777, 0.0151939955182756, 0.0014115438930697, 0.00186964100335118, 
#                   0.0149485506391318, 0.000371853000871041, 0.00126744178851435, 
#                   0.109516216674201, 0.0158416630068266)
beta <- nleqslv(rep(0.1, N), betas)$x
# Converte beta para que P_{AAA} = 1
beta <- beta/colSums(beta * tij)[1]

# Checar se os betas formam a solução
nleqslv(rep(1, N),
        prices)
max(abs(1 - nleqslv(rep(1, N), prices)$x)) < 1e-6


# Índices Baseline ------------------------------------------------------------

# IMR - Baseline
P_BSLN <- colSums(beta * tij)^(1/(1-sigma))
# OMR - Baseline
PI_BSLN <- colSums(t(tij)/P_BSLN^(1-sigma) * E/sum(Y))^(1/(1-sigma))

# Renda Real - RGDP
RGDP_BSLN <- Y/P_BSLN

# Índices - Condicional -------------------------------------------------------
tij <- matrix(dados$t_crfl, nrow = sqrt(nrow(dados)), byrow = TRUE)

# Valores Iniciais
termos_mult <- nleqslv(rep(0.9, N * 2),
                       mrts, control = list(maxit = 2000))$x

# IMR - Condicional
P_CDL <- termos_mult[1:N]^(1/(1-sigma))

# OMR - Condicional
PI_CDL <- termos_mult[(N + 1):(2*N)]^(1/(1-sigma))

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

# Combinar os resultados com os dados originais

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

dados <- dados %>% 
  left_join(resultados, by = c("exporter" = "country")) %>% 
  left_join(resultados, by = c("importer" = "country"),
            suffix = c("_exporter", "_importer"))

dados <- dados %>% 
  mutate(Y_WLD_BSLN = sum(Y_BSLN_exporter * (1- (exporter != importer))),
         Y_WLD_CRFL = sum(Y_CRFL_exporter * (1- (exporter != importer))),
         X_BSLN = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_bsln/(P_BSLN_importer * PI_BSLN_exporter)^(1-sigma),
         X_CDL = (Y_BSLN_exporter * E_BSLN_importer/Y_WLD_BSLN)  * 
           t_crfl/(P_CDL_importer * PI_CDL_exporter)^(1-sigma),
         X_CRFL = (Y_CRFL_exporter * E_CRFL_importer/Y_WLD_CRFL)  * 
           t_crfl/(P_CRFL_importer * PI_CRFL_exporter)^(1-sigma))

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
