# An Advanced Guide to Trade Policy Analysis: The Structural Gravity Model
# Chapter 2: General Equilibrium Trade Policy Analysis with Structural Gravity
# Application 1: Trade without borders
# Last update 2021/03/06


# Load packages -----------------------------------------------------------
library(haven)
library(dplyr)
library(estimatr)
library(lmtest)
library(fixest)
library(ggplot2)
library(tidyr)
library(janitor)

# Plot theme --------------------------------------------------------------

theme_gravity <- function(base_size = 12) {
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    title = element_text(face = "bold"),
    axis.line = element_line(color = "#000000"),
    axis.ticks.length = unit(.2, "cm"),
    legend.position = "bottom",
    legend.key = element_rect(fill = "white"),
    axis.text = element_text(color = "black")
  )
}

# Read data ---------------------------------------------------------------
data_borders <- read_dta("Data/Chapter2Application1.dta") %>%
  clean_names() %>%
  mutate(trade = trade / 1e3)

# Define Germany as reference country
data_borders <- data_borders %>%
  mutate(
    exporter = ifelse(exporter == "DEU", "AAA", exporter),
    importer = ifelse(importer == "DEU", "AAA", importer),
    ln_dist = log(dist),
    intl = ifelse(exporter != importer, 1, 0)
  ) %>%
  filter(year == 2006) %>%
  group_by(importer) %>%
  mutate(e = sum(trade)) %>%
  group_by(exporter) %>%
  mutate(y = sum(trade)) %>%
  ungroup() %>%
  mutate(
    e_r_bln = ifelse(importer == "AAA", e, 0),
    e_r = max(e_r_bln),
    # World output
    y_wld = sum(y * (1 - intl))
  ) %>%
  arrange(importer)


#############################################################################
########################## STEP 1 ##########################################
################# Solve the baseline gravity model ########################
###########################################################################

# Estimate Gravity Equation -----------------------------------------------
fit <- fepois(
  trade ~ ln_dist + cntg + intl | exporter + importer,
  data = data_borders,
  se = "hetero",
  ssc = ssc(adj = FALSE)
)

summary(fit)

# Predicted trade in baseline model
data_borders <- data_borders %>%
  mutate(tradehat_bln = predict(fit)) %>%
  group_by(exporter) %>%
  mutate(xi_bln = sum(tradehat_bln * (exporter != importer))) %>%
  ungroup()

# Fixed effects
data_borders <- data_borders %>%
  mutate(
    fe_exporter = fixef(fit)$exporter[exporter],
    fe_importer = fixef(fit)$importer[importer]
  )

# Trade costs - Baseline
data_borders <- data_borders %>%
  mutate(
    tij_bln = exp(
      fit$coefficients["ln_dist"] * ln_dist +
        fit$coefficients["cntg"] * cntg +
        fit$coefficients["intl"] * intl
    )
  )

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>%
  mutate(
    omr_bln = y * e_r / exp(fe_exporter),
    imr_bln = e / (exp(fe_importer) * e_r)
  )

# Output and Expenditure - Baseline
data_borders <- data_borders %>%
  mutate(
    y_bln = y,
    e_bln = e,
    y_wld_bln = y_wld
  ) %>%
  rename(
    fe_exp_bln = fe_exporter,
    fe_imp_bln = fe_importer
  )

##############################################################################
########################## STEP 2 ###########################################
################ Define the counterfactual scenario ########################
###########################################################################

data_borders <- data_borders %>%
  mutate(
    intl_brdr_cfl = 0,
    tij_cfl = exp(
      fit$coefficients["ln_dist"] * ln_dist +
        fit$coefficients["cntg"] * cntg +
        fit$coefficients["intl"] * intl_brdr_cfl
    )
  )

##############################################################################
########################## STEP 3 ###########################################
################# Solve the counterfactual gravity model #####$##############
###########################################################################

fit_cdl <- fepois(
  trade ~ 0 | exporter + importer,
  offset = ~ log(tij_cfl),
  data = data_borders
)

# Fixed effects
data_borders <- data_borders %>%
  mutate(
    fe_exporter = fixef(fit_cdl)$exporter[exporter],
    fe_importer = fixef(fit_cdl)$importer[importer]
  )

# Outward and Inward multilateral resistance terms
data_borders <- data_borders %>%
  mutate(
    omr_cdl = y * e_r / (exp(fe_exporter)),
    imr_cdl = e / (exp(fe_importer) * e_r)
  )

# Predicted trade - Conditional
data_borders <- data_borders %>%
  mutate(tradehat_cdl = predict(fit_cdl)) %>%
  group_by(exporter) %>%
  mutate(xi_cdl = sum(tradehat_cdl * (exporter != importer))) %>%
  ungroup() %>%
  rename(
    fe_exp_cdl = fe_exporter,
    fe_imp_cdl = fe_importer
  )

# full endownment general equilibrium effects

# The constant of elasticity of substitution
sigma <- 7

# Change in bilateral trade costs
data_borders <- data_borders %>%
  mutate(change_tij = tij_cfl / tij_bln)

# Deficit/Surplus parameter
data_borders <- data_borders %>%
  mutate(phi = ifelse(importer == exporter, e / y, 0)) %>%
  group_by(exporter) %>%
  mutate(phi = max(phi))

# First-order change in prices

# Changes in prices for exporters
data_borders <- data_borders %>%
  mutate(
    change_p_i = ((exp(fe_exp_cdl) / e_r) / (exp(fe_exp_bln) / e_r))^(1 / (1 - sigma))
  )

# Changes in prices for importers
data_borders <- data_borders %>%
  group_by(importer) %>%
  mutate(
    change_p_j = ifelse(importer == exporter, change_p_i, 0),
    change_p_j = max(change_p_j)
  ) %>%
  ungroup()

# Compute change in output and expenditure
data_borders <- data_borders %>%
  mutate(
    y_cfl = y,
    e_cfl = e,
    y_WLD_cfl = sum(y_cfl * (1 - intl)),
    omr_cfl_0 = omr_cdl,
    imr_cfl_0 = imr_cdl,
    e_r_cfl_0 = e_r,
    fe_exp_cfl_0 = fe_exp_cdl,
    fe_imp_cfl_0 = fe_imp_cdl
  )

# Update trade flows - Counterfactual
data_borders <- data_borders %>%
  mutate(trade_cfl = tradehat_cdl * change_p_i * change_p_j)

# Start loop
data_borders <- data_borders %>%
  mutate(
    change_imr_full_0 = 1,
    change_omr_full_0 = 1,
    change_p_i_0 = change_p_i,
    change_p_j_0 = change_p_j,
    tradehat_0 = tradehat_cdl
  )

max_dif <- 1
sd_dif <- 1
change_price_i_old <- 0
iter <- 1
while (sd_dif > 1e-8 | max_dif > 1e-8) {
  data_borders <- data_borders %>%
    mutate(
      trade_1 = tradehat_0 * change_p_i_0 * change_p_j_0 /
        (change_omr_full_0 * change_imr_full_0)
    )

  fit_cfl <- fepois(
    trade_1 ~ 0 | exporter + importer,
    offset = ~ log(tij_cfl),
    data = data_borders,
    glm.iter = 100
  )

  # Fixed effects
  data_borders <- data_borders %>%
    mutate(
      fe_exporter = fixef(fit_cfl)$exporter[exporter],
      fe_importer = fixef(fit_cfl)$importer[importer]
    )

  data_borders <- data_borders %>%
    mutate(tradehat_1 = predict(fit_cfl)) %>%
    group_by(exporter) %>%
    mutate(y_cfl_1 = sum(tradehat_1)) %>%
    ungroup() %>%
    mutate(e_cfl_1 = ifelse(importer == exporter, phi * y_cfl_1, 0)) %>%
    group_by(importer) %>%
    mutate(e_cfl_1 = max(e_cfl_1)) %>%
    ungroup() %>%
    mutate(
      e_r_cfl_1 = ifelse(importer == "AAA", e_cfl_1, 0),
      e_r_cfl_1 = max(e_r_cfl_1)
    )

  # Changes in prices for exporters
  data_borders <- data_borders %>%
    mutate(change_p_i_1 = ((exp(fe_exporter) / e_r_cfl_1) / (exp(fe_exp_cfl_0) / e_r_cfl_0))^(1 / (1 - sigma)))

  # Changes in prices for importers
  data_borders <- data_borders %>%
    group_by(importer) %>%
    mutate(
      change_p_j_1 = ifelse(importer == exporter, change_p_i_1, 0),
      change_p_j_1 = max(change_p_j_1)
    ) %>%
    ungroup()

  # Outward and Inward multilateral resistance terms
  data_borders <- data_borders %>%
    mutate(
      omr_cfl_1 = y_cfl_1 * e_r_cfl_1 / (exp(fe_exporter)),
      imr_cfl_1 = e_cfl_1 / (exp(fe_importer) * e_r_cfl_1)
    )



  # Update diff
  # print(summary(data_borders$change_p_i_0 - change_price_i_old))
  max_dif <- abs(max(data_borders$change_p_i_0 - change_price_i_old))
  sd_dif <- sd(data_borders$change_p_i_0 - change_price_i_old)
  change_price_i_old <- data_borders$change_p_i_0
  cat("Iteration: ", iter, "\n")
  iter <- iter + 1
  cat("Summary of changes in factory gate prices\n")
  cat("Maximum: ", max_dif, "\n")
  cat("Standard-deviation: ", sd_dif, "\n")

  # Changes in omr and imr
  data_borders <- data_borders %>%
    mutate(
      change_imr_full_1 = imr_cfl_1 / imr_cfl_0,
      change_omr_full_1 = omr_cfl_1 / omr_cfl_0,
      imr_cfl_0 = imr_cfl_1,
      omr_cfl_0 = omr_cfl_1,
      change_imr_full_0 = change_imr_full_1,
      change_omr_full_0 = change_omr_full_1,
      change_p_i_0 = change_p_i_1,
      change_p_j_0 = change_p_j_1,
      tradehat_0 = tradehat_1,
      e_r_cfl_0 = e_r_cfl_1,
      fe_exp_cfl_0 = fe_exporter,
      fe_imp_cfl_0 = fe_importer
    ) %>%
    dplyr::select(-fe_exporter, -fe_importer)
}

### CONSERTAR CHANGE E_full
data_borders <- data_borders %>%
  mutate(
    change_p_i_full = ((exp(fe_exp_cfl_0) / e_r_cfl_0) / (exp(fe_exp_bln) / e_r))^(1 / (1 - sigma)),
    change_p_j_full = change_p_i_full * (exporter == importer)
  ) %>%
  group_by(importer) %>%
  mutate(change_p_j_full = max(change_p_j_full)) %>%
  ungroup() %>%
  mutate(
    y_full = change_p_i_full * y_bln,
    E_full = change_p_j_full * e_bln * (exporter == importer)
  ) %>%
  group_by(importer) %>%
  mutate(E_full = max(E_full, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    E_full_R = E_full * (importer == "AAA"),
    E_full_R = max(E_full_R),
    omr_full = y_full * e_r_cfl_0 / exp(fe_exp_cfl_0),
    imr_full = E_full / (exp(fe_imp_cfl_0) * E_full_R),
    X_full = (y_full * E_full * tij_cfl) / (imr_full * omr_full)
  ) %>%
  group_by(exporter) %>%
  mutate(xi_full = sum(X_full * (importer != exporter))) %>%
  ungroup()

# Collect indexes of interest ------------------------------------------------
exporter_indexes <- data_borders %>%
  filter(importer == exporter) %>%
  select(
    exporter, starts_with("omr"), change_p_i_full, starts_with("xi_"),
    y_bln, y_full
  ) %>%
  mutate(
    change_price_full = (change_p_i_full - 1) * 100,
    change_omr_cdl = (omr_cdl^(1 / (1 - sigma)) / omr_bln^(1 / (1 - sigma)) - 1) * 100,
    change_omr_full = (omr_full^(1 / (1 - sigma)) / omr_bln^(1 / (1 - sigma)) - 1) * 100,
    change_xi_cdl = (xi_cdl / xi_bln - 1) * 100,
    change_xi_full = (xi_full / xi_bln - 1) * 100
  ) %>%
  pivot_longer(
    cols = c(-exporter, -y_bln),
    names_to = "variable",
    values_to = "value"
  )

importer_indexes <- data_borders %>%
  filter(importer == exporter) %>%
  select(importer, starts_with("imr"), change_p_i_full, starts_with("y_")) %>%
  mutate(
    change_price_full = (change_p_i_full - 1) * 100,
    change_imr_full = -(imr_full^(1 / (1 - sigma)) / imr_bln^(1 / (1 - sigma)) - 1) * 100,
    rgdp = ((y_full / imr_full^(1 / (1 - sigma))) / (y_bln / imr_bln^(1 / (1 - sigma))) - 1) * 100
  ) %>%
  pivot_longer(
    cols = c(-importer, -y_bln),
    names_to = "variable",
    values_to = "value"
  )


# Visualizations --------------------------------------------------------------

# Figure 6 - page 100
ggplot(
  exporter_indexes %>%
    filter(
      exporter != "HKG",
      variable %in% c("change_xi_full", "change_xi_cdl")
    ),
  aes(x = log(y_bln * 1e3), y = value, color = variable, shape = variable)
) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_color_manual(
    values = c("#280052", "#4A9BE6"),
    labels = c(
      "full Endowment General Equilibrium",
      "Conditional General Equilibrium"
    )
  ) +
  scale_shape_discrete(
    labels = c(
      "full Endowment General Equilibrium",
      "Conditional General Equilibrium"
    )
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
  labs(
    title = "Effects of abolishing international borders on exports",
    y = "Percent change of exports",
    x = "Log value of output",
    color = "",
    shape = ""
  ) +
  theme_gravity()


# Figure 7 - page 100
ggplot(
  importer_indexes %>%
    filter(
      importer != "HKG",
      variable %in% c("rgdp", "change_price_full", "change_imr_full")
    ),
  aes(x = log(y_bln * 1e3), y = value, color = variable, shape = variable)
) +
  geom_point(size = 3.5, alpha = 0.7) +
  scale_color_manual(
    values = c("#280052", "#4A9BE6", "#000878"),
    labels = c(
      "Real GDP",
      "Factory-gate price",
      "-(inward multilateral resistances)"
    )
  ) +
  scale_shape_manual(
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
  labs(
    title = "Effects of abolishing international borders on real GDP",
    y = "Percent change",
    x = "Log value of output",
    color = "",
    shape = ""
  ) +
  theme_gravity()