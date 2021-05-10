# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 2: General Equilibrium Trade Policy Analysis with Structural Gravity
# Application 1: Trade without borders
# Date: 2017/05/10


# Load packages -----------------------------------------------------------
library(haven)
library(dplyr)
library(ggplot2)
library(estimatr)
library(lmtest)
library(fixest)
library(nleqslv)

# Function for general equilibrium analysis
source("R/ge_functions.R", encoding = "UTF-8")

# Plot theme --------------------------------------------------------------

theme_gravity <- function(base_size = 12) {
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    title = element_text(face = "bold"),
    axis.line = element_line(color = "#000000"),
    axis.ticks.length=unit(.2, "cm"),
    legend.position = "bottom",
    legend.key = element_rect(fill = "white"),
    axis.text = element_text(color = "black")
  )
}

# Load the data -------------------------------------------------------------
dados <- read_dta("Data/Chapter2Application1.dta")

# Filter and generate new variables
dados <- dados %>%
  filter(year == 2006) %>%
  mutate(
    ln_dist = log(DIST),
    intl_brdr = ifelse(exporter != importer, 1, 0)
  ) %>%
  rename(cntg = CNTG) %>% 
  group_by(exporter) %>%
  mutate(Y = sum(trade)) %>%
  ungroup() %>%
  group_by(importer) %>%
  mutate(E = sum(trade)) %>%
  ungroup()

# Define Germany as reference country
dados <- dados %>%
  mutate(
    exporter = ifelse(exporter == "DEU", "AAA", exporter),
    importer = ifelse(importer == "DEU", "AAA", importer)
  )

# Estimate the gravity equation
fit <- fepois(
  trade ~ ln_dist + cntg + intl_brdr | importer + exporter,
  data = dados %>%
    arrange(exporter, importer)
)

# Compute trade costs - baseline and counterfactual
dados <- dados %>%
  mutate(
    # t_bln is tc^(1-sigma); tc = trade costs
    t_bln = exp(
      fit$coefficients["ln_dist"] * ln_dist +
        fit$coefficients["cntg"] * cntg +
        fit$coefficients["intl_brdr"] * intl_brdr
    ),
    t_cfl = exp(
      fit$coefficients["ln_dist"] * ln_dist +
        fit$coefficients["cntg"] * cntg +
        fit$coefficients["intl_brdr"] * intl_brdr * 0
    )
  ) %>%
  select(exporter, importer, trade, t_bln, t_cfl) %>%
  arrange(exporter, importer)


# Compute Y, E, Q, phi, trade costs (tc) and alpha -----------------------------

Y <- dados %>%
  group_by(exporter) %>%
  summarise(Y = sum(trade), .groups = "drop") %>%
  pull(Y)

E <- dados %>%
  group_by(importer) %>%
  summarise(E = sum(trade), .groups = "drop") %>%
  pull(E)

Q <- Y

phi <- E / Q

N <- sqrt(nrow(dados))

X <- fit$fitted.values
X <- matrix(X, N, byrow = TRUE)

sigma <- 7

tc_bln <- matrix(dados$t_bln, nrow = N, byrow = TRUE)^(1 / (1 - sigma))
tc_cfl <- matrix(dados$t_cfl, nrow = N, byrow = TRUE)^(1 / (1 - sigma))

#alpha <- exp(fixef(fit)$exporter)^(1 / (1 - sigma))

P <- (E / exp(fixef(fit)$importer))^(1 / (1 - sigma))
P <- P / P[1]
PI <- P
for(i in 1:N){
  PI[i] <- sum((tc_bln[i, ] / P)^(1 - sigma) * E / sum(Y))^(1 / (1 - sigma))
}

alpha <- (Y/sum(Y))^(1/(1-sigma)) * 1/PI

# Normalize alpha to make P[AAA] = 1
#alpha <- alpha / P[1]
# Compute  - Baseline ----------------------------------------------------------

vars_bln <- ge_gravity(phi, Q, alpha, tc_bln, N, full_endowment = TRUE)

P_bln <- vars_bln$P
PI_bln <- vars_bln$PI
RGDP_bln <- vars_bln$RGDP
p_bln <- vars_bln$p
exp_bln <- rowSums(vars_bln$X) - diag(vars_bln$X)

# Compute  - Conditional  ------------------------------------------------------
vars_cdl <- ge_gravity(phi, Q, alpha, tc_cfl, N, full_endowment = FALSE)

P_cdl <- vars_cdl$P
PI_cdl <- vars_cdl$PI
exp_cdl <- rowSums(vars_cdl$X) - diag(vars_cdl$X)

# Compute: Counterfactual ----------------------------------------------------
vars_cfl <- ge_gravity(phi, Q, alpha, tc_cfl, N, full_endowment = TRUE)

P_cfl <- vars_cfl$P
PI_cfl <- vars_cfl$PI
RGDP_cfl <- vars_cfl$RGDP
p_cfl <- vars_cfl$p
exp_cfl <- rowSums(vars_cfl$X) - diag(vars_cfl$X)


# Merge results
# Indexes by country

results <- dados %>%
  select(country = exporter) %>%
  distinct() %>%
  mutate(
    country = ifelse(country == "AAA", "DEU", country),
    p_bln = p_bln,
    P_bln = P_bln,
    PI_bln = PI_bln,
    RGDP_bln = RGDP_bln,
    exp_bln = exp_bln,
    exp_cdl = exp_cdl,
    p_cfl = p_cfl,
    P_cfl = P_cfl,
    PI_cfl = PI_cfl,
    RGDP_cfl = RGDP_cfl,
    exp_cfl = exp_cfl
  ) %>%
  arrange(country) %>% 
  mutate(
    chng_exp_cdl = round((exp_cdl / exp_bln - 1) * 100, 2),
    chng_exp_cfl = round((exp_cfl / exp_bln - 1) * 100, 2),
    chng_RGDP = round((RGDP_cfl / RGDP_bln - 1) * 100, 2),
    chng_IMR = round((P_cfl / P_bln - 1) * 100, 2),
    chng_OMR = round((PI_cfl / PI_bln - 1) * 100, 2),
    chng_p = round((p_cfl / p_bln - 1) * 100, 2)
  ) %>% 
  select(country, starts_with("chng"), RGDP_bln) %>%
  as.data.frame()

# Plots -----------------------------------------------------------------------

# Figure 6 - page 100
ggplot(
  results,
  aes(x = log(RGDP_bln))
) +
  geom_point(aes(y = chng_exp_cfl, color = "a", shape = "a"), size = 3.5, alpha = 0.7) +
  geom_point(aes(y = chng_exp_cdl, color = "b", shape = "b"), size = 3.5, alpha = 0.7) +
  labs(
    title = "Effects of abolishing international borders on exports",
    y = "Percent change of exports",
    x = "Log value of output"
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(6, 16, 1),
    limits = c(6, 16)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 450)
  ) +
  scale_color_manual(
    values = c("#280052", "#4A9BE6"),
    labels = c(
      "Full Endowment General Equilibrium",
      "Conditional General Equilibrium"
    )
  ) +
  scale_shape_discrete(
    labels = c(
      "Full Endowment General Equilibrium",
      "Conditional General Equilibrium"
    )
  ) +
  labs(
    color = "",
    shape = ""
  ) +
  theme_gravity()

ggsave(
  "figures/figure6.png",
  width = 8,
  height = 6,
  dpi = "print"
)

# Figure 7 - page 100
ggplot(
  results %>% filter(!country %in% c("HKG", "NER")),
  aes(x = log(RGDP_bln))
) +
  geom_point(
    aes(y = chng_RGDP, color = "a", shape = "a"),
    size = 3.5,
    alpha = 0.7
  ) +
  geom_point(
    aes(y = chng_p, color = "b", shape = "b"),
    size = 3.5,
    alpha = 0.7
  ) +
  geom_point(
    aes(y = -chng_IMR, color = "c", shape = "c"),
    size = 3.5,
    alpha = 0.7
  ) +
  labs(
    title = "Effects of abolishing international borders on real GDP",
    y = "Percent change",
    x = "Log value of output"
  ) +
  scale_color_manual("",
    values = c("#280052", "#4A9BE6", "#000878"),
    labels = c(
      "Real GDP",
      "Factory-gate price",
      "-(inward multilateral resistances)"
    )
  ) +
  scale_shape_manual("",
    values = 15:17,
    labels = c(
      "Real GDP",
      "Factory-gate price",
      "-(inward multilateral resistances)"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0),
    breaks = seq(6, 16, 1),
    limits = c(6, 16)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = seq(-20, 80, 20),
    limits = c(-20, 80)
  ) +
  theme_gravity()

ggsave(
  "figures/figure7.png",
  width = 8,
  height = 6,
  dpi = "print"
)