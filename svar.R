# ==========================
# 📦 Load Libraries
# ==========================
library(fredr)
library(dplyr)
library(zoo)
library(purrr)
library(ggplot2)
library(corrr)
library(vars)
library(svars)
library(tseries)      # For ADF test
library(urca)         # For cointegration test

# ==========================
# 🔑 Set Your FRED API Key
# ==========================
fredr_set_key("7167653a5bb571e358c87d7a5d466bf8")

# ==========================
# 🗕️ Set Date Range
# ==========================
start_date <- as.Date("2000-01-01")
end_date <- as.Date("2019-12-31")

# ==========================
# 🧾 Download Monthly Series First (Need Averaging)
# ==========================
unrate <- fredr(series_id = "UNRATE", observation_start = start_date, observation_end = end_date)
unrate_q <- unrate %>% mutate(quarter = as.yearqtr(date)) %>% group_by(quarter) %>% summarise(unrate = mean(value))

employment <- fredr(series_id = "CE16OV", observation_start = start_date, observation_end = end_date)
employment_q <- employment %>% mutate(quarter = as.yearqtr(date)) %>% group_by(quarter) %>% summarise(employment = mean(value))

population <- fredr(series_id = "CNP16OV", observation_start = start_date, observation_end = end_date)
population_q <- population %>% mutate(quarter = as.yearqtr(date)) %>% group_by(quarter) %>% summarise(population = mean(value))

# ==========================
# 🧾 Download Quarterly Series
# ==========================
gdp <- fredr(series_id = "GDPC1", observation_start = start_date, observation_end = end_date) %>% mutate(quarter = as.yearqtr(date)) %>% select(quarter, gdp = value)

spending <- fredr(series_id = "FGEXPND", observation_start = start_date, observation_end = end_date) %>% mutate(quarter = as.yearqtr(date)) %>% select(quarter, gov_spending_nominal = value)

taxes <- fredr(series_id = "FGRECPT", observation_start = start_date, observation_end = end_date) %>% mutate(quarter = as.yearqtr(date)) %>% select(quarter, tax_revenue_nominal = value)

deflator <- fredr(series_id = "GDPDEF", observation_start = start_date, observation_end = end_date) %>% mutate(quarter = as.yearqtr(date)) %>% select(quarter, gdp_deflator = value)

# ==========================
# 🔗 Merge All Data
# ==========================
data_all <- reduce(list(gdp, spending, taxes, deflator, unrate_q, employment_q, population_q), full_join, by = "quarter") %>% arrange(quarter)

# ==========================
# ⚙️ Real Terms & Log Transforms
# ==========================
data_all <- data_all %>% mutate(
  gov_spending_real = gov_spending_nominal / (gdp_deflator / 100),
  tax_revenue_real = tax_revenue_nominal / (gdp_deflator / 100),
  log_gdp = log(gdp),
  log_spending = log(gov_spending_real),
  log_taxes = log(tax_revenue_real)
)

# ==========================
# 🦯 Create First Differences
# ==========================
data_diff <- data_all %>%
  mutate(
    d_log_gdp      = log_gdp      - lag(log_gdp, 1),
    d_log_spending = log_spending - lag(log_spending, 1),
    d_log_taxes    = log_taxes    - lag(log_taxes, 1),
    arra_dummy     = if_else(quarter >= as.yearqtr("2009 Q1") & quarter <= as.yearqtr("2009 Q3"), 1, 0)
  ) %>%
  filter(!is.na(d_log_gdp))

# ==========================
# 🎯 Lag Selection for Differenced Data
# ==========================
var_data_diff <- data_diff %>% select(d_log_spending, d_log_taxes, d_log_gdp)
VARselect(var_data_diff, lag.max = 8, type = "const")

# ==========================
# 🤜 Adjust Taxes for Automatic Stabilizers (Blanchard-Perotti)
# ==========================
data_diff <- data_diff %>% mutate(d_log_taxes_adjusted = d_log_taxes - 2.08 * d_log_gdp)

var_data_adj <- data_diff %>% select(d_log_spending, d_log_taxes_adjusted, d_log_gdp)
VAR_model_adj <- VAR(var_data_adj, p = 2, type = "const")
SVAR_adj <- id.chol(VAR_model_adj)
irf_spend_adj <- irf(SVAR_adj, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12)
plot(irf_spend_adj)

# ==========================
# 📊 VAR with ARRA Dummy (IRFs Only, No SVAR)
# ==========================
aligned_data <- data_diff %>%
  select(d_log_spending, d_log_taxes, d_log_gdp, arra_dummy) %>%
  na.omit()

var_data_diff <- aligned_data %>% select(d_log_spending, d_log_taxes, d_log_gdp)
arra_dummy_df <- aligned_data %>% select(arra_dummy)

VAR_model_diff <- VAR(
  var_data_diff,
  p = 1,
  type = "const",
  exogen = arra_dummy_df
)

irf_spend_diff <- irf(VAR_model_diff, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12, boot = TRUE)
plot(irf_spend_diff)

irf_tax_diff <- irf(VAR_model_diff, impulse = "d_log_taxes", response = "d_log_gdp", n.ahead = 12, boot = TRUE)
plot(irf_tax_diff)

# ==========================
# 📈 GDP Growth Visualization
# ==========================
ggplot(data_diff, aes(x = quarter, y = d_log_gdp)) +
  geom_line(color = "steelblue") +
  labs(title = "GDP Growth Rate (∆log(GDP))", x = "Quarter", y = "∆log(GDP)") +
  theme_minimal()

ggplot(data_diff, aes(x = quarter, y = d_log_gdp)) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = as.yearqtr("2009 Q1"), linetype = "dashed", color = "red") +
  geom_vline(xintercept = as.yearqtr("2009 Q3"), linetype = "dashed", color = "red") +
  labs(title = "GDP Growth Rate with ARRA Window", y = "∆log(GDP)", x = "Quarter") +
  theme_minimal()

#pangay

# ==========================
# 📊 ROBUSTNESS + EXTENSIONS (Post-Estimation)
# ==========================
library(tseries)
library(tsDyn)

# ==========================
# 🌐 Alternative Identification: Reverse Cholesky
# ==========================
var_data_rev <- data_diff %>% select(d_log_taxes, d_log_spending, d_log_gdp)
VAR_rev <- VAR(var_data_rev, p = 1, type = "const")
SVAR_rev <- id.chol(VAR_rev)
irf_rev <- irf(SVAR_rev, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12)
plot(irf_rev)

# ==========================
# ⏳ Lag Sensitivity: p = 1 and p = 3
# ==========================
VAR_lag1 <- VAR(var_data_diff, p = 1, type = "const")
VAR_lag3 <- VAR(var_data_diff, p = 3, type = "const")
irf_lag1 <- irf(VAR_lag1, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_lag1)
irf_lag3 <- irf(VAR_lag3, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_lag3)

# ==========================
# ⚖️ Subsample Stability
# ==========================
data_pre <- data_diff %>% filter(quarter < as.yearqtr("2008 Q1"))
data_post <- data_diff %>% filter(quarter >= as.yearqtr("2010 Q1"))

VAR_pre <- VAR(data_pre %>% select(d_log_spending, d_log_taxes, d_log_gdp), p = 1)
VAR_post <- VAR(data_post %>% select(d_log_spending, d_log_taxes, d_log_gdp), p = 1)

irf_pre <- irf(VAR_pre, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_pre)

irf_post <- irf(VAR_post, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_post)

# ==========================
# 🔍 Stationarity Confirmations
# ==========================
PP.test(data_diff$d_log_gdp)
PP.test(data_diff$d_log_spending)
PP.test(data_diff$d_log_taxes)

# ==========================
# ⚡ Add Control Variables
# ==========================
fedfunds <- fredr(series_id = "FEDFUNDS", observation_start = start_date, observation_end = end_date) %>%
  mutate(quarter = as.yearqtr(date)) %>%
  group_by(quarter) %>%
  summarise(fedfunds = mean(value))

sentiment <- fredr(series_id = "UMCSENT", observation_start = start_date, observation_end = end_date) %>%
  mutate(quarter = as.yearqtr(date)) %>%
  group_by(quarter) %>%
  summarise(sentiment = mean(value))

macro_controls <- reduce(list(fedfunds, sentiment), full_join, by = "quarter")
data_diff <- left_join(data_diff, macro_controls, by = "quarter")

aligned_macro <- data_diff %>%
  select(d_log_spending, d_log_taxes, d_log_gdp, fedfunds, sentiment) %>%
  na.omit()

VAR_macro <- VAR(
  dplyr::select(aligned_macro, d_log_spending, d_log_taxes, d_log_gdp),
  p = 1,
  type = "const",
  exogen = dplyr::select(aligned_macro, fedfunds, sentiment)
)
irf_macro <- irf(VAR_macro, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_macro)

# ==========================
# 🤓 Spending Elasticity Adjustment (BP 2002 style)
# ==========================
data_diff <- data_diff %>%
  mutate(d_log_spending_adjusted = d_log_spending - 0.1 * d_log_gdp)

VAR_spend_adj <- VAR(data_diff %>% select(d_log_spending_adjusted, d_log_taxes, d_log_gdp), p = 1)
SVAR_spend_adj <- id.chol(VAR_spend_adj)
irf_spend_adj2 <- irf(SVAR_spend_adj, impulse = "d_log_spending_adjusted", response = "d_log_gdp")
plot(irf_spend_adj2)

# ==========================
# ⏰ Anticipation Effect: ARRA Lead
# ==========================
data_diff <- data_diff %>% mutate(arra_lead = lead(arra_dummy, 1))

anticipation_data <- data_diff %>%
  select(d_log_spending, d_log_taxes, d_log_gdp, arra_lead) %>%
  na.omit()

VAR_anticipate <- VAR(
  dplyr::select(anticipation_data, d_log_spending, d_log_taxes, d_log_gdp),
  p = 1,
  type = "const",
  exogen = dplyr::select(anticipation_data, arra_lead)
)
irf_anticipate <- irf(VAR_anticipate, impulse = "d_log_spending", response = "d_log_gdp")
plot(irf_anticipate)

# ==========================
# 🚀 Nonlinear Effects: Threshold VAR
# ==========================
tvar_model <- TVAR(var_data_diff, lag = 2, nthresh = 1, thDelay = 1, trim = 0.1)
summary(tvar_model)

# ==========================
# 🔁 Residual Autocorrelation Check (Portmanteau Test)
# ==========================
serial_test <- serial.test(VAR_model_diff, lags.pt = 12, type = "PT.asymptotic")
print(serial_test)

# ==========================
# 📊 Residual Normality Test (Jarque-Bera)
# ==========================
norm_test <- normality.test(VAR_model_diff)
print(norm_test)

# ==========================
# ⚡ ARCH Effects (Heteroskedasticity Test)
# ==========================
arch_test <- arch.test(VAR_model_diff, lags.multi = 5)
print(arch_test)

# ==========================
# 🔒 Stability Check (Roots of Characteristic Polynomial)
# ==========================
stability <- roots(VAR_model_diff)
print(stability)       # Display eigenvalues
modulus_check <- all(Mod(stability) < 1)
cat("VAR Stability Condition Satisfied:", modulus_check, "\n")

# ==========================
# 🔀 Granger Causality: Fiscal Policy → GDP
# ==========================
granger_spending <- causality(VAR_model_diff, cause = "d_log_spending")
granger_taxes <- causality(VAR_model_diff, cause = "d_log_taxes")

print(granger_spending)
print(granger_taxes)

# ==========================
# 📊 Forecast Error Variance Decomposition (FEVD)
# ==========================
fevd_result <- fevd(VAR_model_diff, n.ahead = 12)
plot(fevd_result)

# ==========================
# 📈 Cumulative IRFs
# ==========================
c_irf_spend <- irf(VAR_model_diff, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_spend)

c_irf_tax <- irf(VAR_model_diff, impulse = "d_log_taxes", response = "d_log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_tax)

# ==========================
# 📈 Comparative IRFs: GDP Response to Spending Shocks
# ==========================

library(tidyr)  # Needed for drop_na()

# Helper to extract IRFs into a tidy format
extract_irf <- function(irf_obj, label) {
  tibble(
    quarter = 1:length(irf_obj$irf$d_log_spending),
    response = as.vector(irf_obj$irf$d_log_spending),
    model = label
  )
}

# Collect IRFs from all your earlier objects
irf_base      <- extract_irf(irf(VAR_model_diff, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "ARRA Dummy (Baseline)")
irf_svar      <- extract_irf(irf(SVAR_adj, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "SVAR (Adj. Taxes)")
irf_rev       <- extract_irf(irf(SVAR_rev, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "Reverse Ordering")
irf_pre       <- extract_irf(irf(VAR_pre, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "Pre-GFC")
irf_post      <- extract_irf(irf(VAR_post, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "Post-GFC")
irf_controls  <- extract_irf(irf(VAR_macro, impulse = "d_log_spending", response = "d_log_gdp", n.ahead = 12), "With Controls")

# Combine and clean
irf_combined <- bind_rows(irf_base, irf_svar, irf_rev, irf_pre, irf_post, irf_controls) %>%
  drop_na(response)

# Plot
ggplot(irf_combined, aes(x = quarter, y = response, color = model)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Impulse Responses: Government Spending Shock → GDP",
    x = "Quarter", y = "∆log(GDP)",
    color = "Specification"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# ==========================
# 📈 Residual ACF Plots
# ==========================
acf(resid(VAR_model_diff)[, "d_log_gdp"], main = "ACF of GDP Residuals")
acf(resid(VAR_model_diff)[, "d_log_spending"], main = "ACF of Spending Residuals")
acf(resid(VAR_model_diff)[, "d_log_taxes"], main = "ACF of Tax Residuals")

#final

# ==========================
# 📊 Clean Summary Statistics Table
# ==========================
summary_stats <- data_diff %>%
  select(d_log_gdp, d_log_spending, d_log_taxes, fedfunds, sentiment) %>%
  pivot_longer(cols = everything(), names_to = "variable") %>%
  group_by(variable) %>%
  reframe(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    min  = min(value, na.rm = TRUE),
    max  = max(value, na.rm = TRUE)
  )

print(summary_stats)

# ==========================
# 🔍 VAR Lag Selection Criteria Table
# ==========================
VARselect(var_data_diff, lag.max = 8, type = "const")$criteria

# ==========================
# 🔍 DF-GLS Test for Unit Roots
# ==========================
library(urca)

summary(ur.ers(data_all$log_gdp, type = "DF-GLS", model = "trend", lag.max = 4))
summary(ur.ers(data_all$log_spending, type = "DF-GLS", model = "trend", lag.max = 4))
summary(ur.ers(data_all$log_taxes, type = "DF-GLS", model = "trend", lag.max = 4))

# ==========================
# 🔗 Johansen Cointegration Test
# ==========================
library(urca)
library(dplyr)  # For %>%, select(), etc.

coint_test <- ca.jo(data_all %>% select(log_gdp, log_spending, log_taxes) %>% na.omit(),
                    type = "trace", ecdet = "const", K = 2)
summary(coint_test)

#==== Moment of realization ====

# ==========================
# 🧠 Convert VECM to VAR Representation
# ==========================
library(vars)

vecm_irf_model <- vec2var(coint_test, r = 1)

# Plot Impulse Response
plot(irf(vecm_irf_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12))

# ==========================
# 📉 VECM Residual Diagnostics
# ==========================
# Portmanteau test for residual autocorrelation
serial.test(vecm_irf_model, lags.pt = 12, type = "PT.asymptotic")

# Jarque-Bera test for residual normality
normality.test(vecm_irf_model)

# ARCH-LM test for heteroskedasticity
arch.test(vecm_irf_model, lags.multi = 5)

# Stability check (roots of characteristic polynomial)
vecm_roots <- roots(vecm_irf_model)
print(vecm_roots)
cat("VECM Stability Condition Satisfied:", all(Mod(vecm_roots) < 1), "\n")

# ==========================
# 📊 Forecast Error Variance Decomposition (FEVD) for VECM
# ==========================
fevd_vecm <- fevd(vecm_irf_model, n.ahead = 12)
plot(fevd_vecm)

# ==========================
# 📈 Cumulative IRFs from VECM
# ==========================
cumulative_irf_spending <- irf(vecm_irf_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(cumulative_irf_spending)

cumulative_irf_taxes <- irf(vecm_irf_model, impulse = "log_taxes", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(cumulative_irf_taxes)

# ==========================
# 💡 Subsample Robustness Using VECM
# ==========================

# Pre-GFC VECM (before 2008 Q1)
vecm_pre <- ca.jo(
  data_all %>% 
    filter(quarter < as.yearqtr("2008 Q1")) %>%
    dplyr::select(log_gdp, log_spending, log_taxes) %>%
    na.omit(),
  type = "trace", ecdet = "const", K = 2
)
summary(vecm_pre)

vecm_pre_model <- vec2var(vecm_pre, r = 1)
plot(irf(vecm_pre_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12))

# Post-GFC VECM (from 2010 Q1 onwards)
vecm_post <- ca.jo(
  data_all %>%
    filter(quarter >= as.yearqtr("2010 Q1")) %>%
    dplyr::select(log_gdp, log_spending, log_taxes) %>%
    na.omit(),
  type = "trace", ecdet = "const", K = 2
)
summary(vecm_post)

vecm_post_model <- vec2var(vecm_post, r = 1)
plot(irf(vecm_post_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12))


# ==========================
# 🧠 Summary of VECM Implications
# ==========================

cat("Key Takeaways:\n")
cat("- Johansen test confirms 1 cointegrating relationship: levels of log(GDP), log(spending), and log(taxes) are not independent.\n")
cat("- VECM replaces differenced VAR to retain long-run equilibrium structure.\n")
cat("- Impulse responses show cumulative effects of fiscal shocks on log(GDP).\n")
cat("- Diagnostics: No autocorrelation (Portmanteau test), but normality and homoskedasticity assumptions are rejected. Robust standard errors advisable.\n")
cat("- FEVD shows spending and tax shocks both explain meaningful variation in GDP.\n")

# Cumulative IRF - Pre-GFC
c_irf_pre <- irf(vecm_pre_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_pre)

# Cumulative IRF - Post-GFC
c_irf_post <- irf(vecm_post_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_post)

# Cumulative IRF - Pre-GFC
c_irf_pre <- irf(vecm_pre_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_pre)

# Cumulative IRF - Post-GFC
c_irf_post <- irf(vecm_post_model, impulse = "log_spending", response = "log_gdp", n.ahead = 12, cumulative = TRUE)
plot(c_irf_post)
 




