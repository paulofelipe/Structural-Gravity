library(haven)
library(lfe)
library(dplyr)
library(BB)

source("R/aux_functions.R")


# Função para criar os termos multilaterais de resistência ----------------

# Observação: O primeiro país (em ordem alfabética deve ter sido usado como
# referência na estimação

MRTS <- function(p){
  r <- rep(NA, length(p))
  n <- length(p)/2
  idx_P <- 2:n
  idx_Pi <- (n+1):(2*n)
  
  Pi <- matrix(rep(p[idx_Pi], n), n, byrow = FALSE) 
  P <- matrix(rep(p[1:n], n), n, byrow = TRUE) 
  
  # IMR
  r[1] <- 1 - p[1]
  r[idx_P] <- p[idx_P] - as.vector((t((tc / Pi)) %*% out_share)[-1,])
  # OMR
  r[idx_Pi] <- p[idx_Pi] - as.vector((tc / P) %*% exp_share)
  
  r
}

get_alphas <- function(p){
  n <- length(p)
  r <- rep(NA, n)
  
  n <- length(p)/3
  idx_P <- 2:n
  idx_Pi <- (n+1):(2*n)
  idx_alphas <- (2*n + 1):(3*n)
  
  # Variáveis Endógenas
  Pi <- matrix(rep(p[idx_Pi], n), n, byrow = FALSE) 
  P <- matrix(rep(p[1:n], n), n, byrow = TRUE) 
  
  # IMR
  r[1] <- 1 - p[1]
  r[idx_P] <- p[idx_P] - as.vector((t((tc / Pi)) %*% out_share)[-1,])
  
  # OMR
  r[idx_Pi] <- p[idx_Pi] - as.vector((tc / P) %*% exp_share)
  
  # Alphas
  r[idx_alphas] <- p[idx_alphas] - as.vector(out_share) * 1/as.vector(Pi[,1])
  
  r
}

get_prices <- function(p){
  r <- rep(NA, length(p))
  n <- length(p)/3
  idx_P <- 2:n
  idx_Pi <- (n+1):(2*n)
  idx_price <- (2*n + 1):(3*n)
  
  # Variáveis Endógenas
  Pi <- matrix(rep(p[idx_Pi], n), n, byrow = FALSE) 
  P <- matrix(rep(p[1:n], n), n, byrow = TRUE) 
  price <- p[idx_price]
  out_share <- price * Q / sum(price * Q)
  out_share <- matrix(out_share, nrow = n)
  exp_share <- (phi * price * Q) / sum(phi * price * Q)
  exp_share <- matrix(exp_share, nrow = n)
  
  # IMR
  r[1] <- 1 - p[1]
  r[idx_P] <- p[idx_P] - as.vector((t((tc / Pi)) %*% out_share)[-1,])
  
  # OMR
  r[idx_Pi] <- p[idx_Pi] - as.vector((tc / P) %*% exp_share)
  
  # Price
  r[idx_price] <- p[idx_price]^(1-sigma) - as.vector(out_share) * 1/(alphas * as.vector(Pi[,1]))
  
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


head(dados)

# Dados usados no sistema ----------------------------------------------------

# Matriz de custos de comércio - Baseline
tc_bsln <- matrix(dados$t_bsln, nrow = sqrt(nrow(dados)), byrow = TRUE)

# Termos multilaterais de resistência: contra factual
tc_crfl <- matrix(dados$t_crfl, nrow = sqrt(nrow(dados)), byrow = TRUE)

# Dados de produção 
y_i <- dados %>% 
  group_by(exporter) %>% 
  summarise(y_i = sum(trade))

# Dados de dispêndio
e_i <- dados %>% 
  group_by(importer) %>% 
  summarise(e_i = sum(trade))

# Dados de participação na produção
out_share <- matrix(y_i$y_i/sum(y_i$y_i), nrow = nrow(y_i))

# Dados de participação no dispêndio
exp_share <- matrix(e_i$e_i/sum(e_i$e_i), nrow = nrow(e_i))

# Dados de dotação
Q <- matrix(y_i$y_i, nrow = nrow(y_i))

# Parâmetro phi
phi <- left_join(y_i, e_i, by = c("exporter" = "importer"))
phi <- phi %>% 
  mutate(phi = e_i/y_i)
phi <- phi$phi

# Elasticidade de substituição
sigma <- 7

# Soluções -------------------------------------------------------------------

# Termos multilaterais de resistência: baseline
tc <- tc_bsln
p <- rep(1, nrow(tc)*2)
sol <- dfsane(p, MRTS, control = list(tol = 1e-12, maxit = 10000, M = 0, noimp = 1000))
mrt_bsln <- sol$par

# Recuperar os alphas
p <- c(mrt_bsln, rep(10, nrow(tc)))
sol <- dfsane(p, get_alphas, control = list(tol = 1e-12, maxit = 10000, M = 300, noimp = 1000))
alphas <- sol$par[(nrow(tc)*2 + 1):(nrow(tc)*3)]

# Termos multilaterais de resistência: contra factual
tc <- tc_crfl
p <- mrt_bsln
sol <- dfsane(p, MRTS, control = list(tol = 1e-12, maxit = 10000, M = 300, noimp = 1000))
mrt_crfl <- sol$par

# Checar se os alphas são a solução para custos base
p <- c(mrt_bsln, rep(1, nrow(tc)) + 0.1)
tc <- tc_bsln
sol <- dfsane(p, get_prices, control = list(tol = 1e-12, maxit = 10000, M = 300, noimp = 1000))

# Solução com os custos de comércio contra factual
tc <- tc_crfl
p <- c(mrt_crfl, rep(1, nrow(tc)))
sol <- dfsane(p, get_prices, control = list(tol = 1e-12, maxit = 10000, M = 0, noimp = 1000))
sol

## Comparar resultados

indexes <- read_dta('D:/Users/paulo.alencar/Downloads/GE PPML/all_indexes_geppml.dta')
data <- data.frame(nls = sol$par[83:123], ppml = indexes$p_full[c(41,1:40)])
library(ggplot2)
ggplot(data, aes(x = nls, y = ppml)) + geom_point(size = 3) + 
  geom_abline(col = "red")
