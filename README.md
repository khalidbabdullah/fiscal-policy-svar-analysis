# U.S. Fiscal Policy Impact Analysis: VECM and SVAR Models

A rigorous econometric study examining the dynamic effects of U.S. government spending and taxation on real GDP from 2000-2019. This research employs advanced time-series methods including Vector Error Correction Models (VECM) and Structural Vector Autoregression (SVAR) to identify fiscal multipliers and analyze policy effectiveness before and after the 2008 Global Financial Crisis.

## 📊 Research Overview

**Period:** 2000-2019 (quarterly data, 80 observations)  
**Data Source:** Federal Reserve Economic Data (FRED)  
**Methods:** VECM, SVAR, Johansen Cointegration, Chow Structural Break Test  
**Tools:** R (vars, urca, strucchange packages)

## 🎯 Research Questions

1. Do government spending and tax revenue have a long-run equilibrium relationship with GDP?
2. What is the dynamic impact of fiscal shocks on output in the short run?
3. Did fiscal policy effectiveness change after the 2008 financial crisis?
4. Are tax multipliers larger than spending multipliers?

## 🔑 Key Findings

### Long-Run Cointegration
- **Strong evidence** of one cointegrating relationship among GDP, government spending, and tax revenue (Johansen trace test: p < 0.001)
- Variables return to equilibrium after shocks, validating VECM framework

### Fiscal Multiplier Estimates
- **Spending shocks**: Modest positive effect on GDP, peaking in Q2 but dissipating quickly
- **Tax shocks**: Larger and more persistent negative impact on output
- **Statistical significance**: Both effects remain statistically insignificant at conventional levels (95% bootstrap CIs include zero)

### Structural Break: Pre- vs. Post-Crisis
- **Chow test** confirms significant structural break around 2008-2009 (supF = 107.85, p < 0.001)
- **Pre-GFC (2000-2007)**: Spending shocks weak, tax shocks dominant
- **Post-GFC (2010-2019)**: Spending multipliers larger and more persistent, suggesting higher fiscal effectiveness during economic slack

### Variance Decomposition
- GDP forecast errors primarily explained by own shocks (~90%)
- Fiscal variables contribute modestly (~10% combined by Q12)
- Suggests limited short-run fiscal transmission, consistent with weak IRF estimates

## 📈 Methodology

### Unit Root Testing
- Elliott-Rothenberg-Stock DF-GLS test with trend
- All variables non-stationary in levels → supports cointegration analysis

### Cointegration Analysis
- Johansen trace test with 2 lags
- Finds r = 1 cointegrating vector

### VECM Specification
ΔYₜ = ΠYₜ₋₁ + Σ Γᵢ ΔYₜ₋ᵢ + μ + εₜ
where Π = αβ′ captures long-run equilibrium

### SVAR Identification
- Recursive (Cholesky) ordering: Spending → Taxes → GDP
- Assumes spending does not respond contemporaneously to GDP or taxes
- Follows Blanchard-Perotti (2002) framework

### Robustness Checks
- Bootstrap confidence intervals (100 runs) for non-normal residuals
- Portmanteau test: No autocorrelation (p = 0.564)
- ARCH-LM test: Heteroskedasticity detected (p < 0.001) → bootstrap inference used
- Subsample analysis: Pre-GFC vs. Post-GFC regimes

## 📁 Repository Contents

- `svar.R` - Complete R script with data import, model estimation, diagnostics, and visualization
- `U_S__Fiscal_Policy___Output.pdf` - Full research paper (21 pages) with literature review, methodology, results, and appendices

## 💡 Technical Skills Demonstrated

- **Econometric Methods**: VECM, SVAR, cointegration testing, impulse response functions, forecast error variance decomposition
- **Statistical Testing**: Unit root tests (DF-GLS), Johansen cointegration, Chow structural break, residual diagnostics (Portmanteau, ARCH-LM, Jarque-Bera)
- **Programming in R**: Time-series analysis with vars, urca, strucchange packages; data manipulation with fredr and dplyr; visualization with ggplot2
- **Research Design**: Subsample analysis, bootstrap inference, model comparison

## 📚 Key References

- Blanchard & Perotti (2002) - SVAR identification strategy
- Johansen (1991) - Cointegration methodology
- Auerbach & Gorodnichenko (2012) - State-dependent fiscal multipliers
- Romer & Romer (2010) - Narrative approach to tax shocks
- Ramey (2011) - Timing of fiscal policy shocks

## 🎓 Academic Context

This paper was completed as part of an advanced econometrics course at Knox College (June 2025). It demonstrates mastery of time-series econometrics, structural identification strategies, and empirical macroeconomic analysis. The research contributes to ongoing debates about fiscal policy effectiveness and highlights how multipliers vary across economic regimes.

## ⚠️ Limitations

- Main impulse responses statistically insignificant (confidence intervals include zero)
- SVAR recursive identification may miss contemporaneous feedback effects
- Omitted variables: interest rates, global shocks, private sector debt dynamics
- Small subsample sizes (~32 quarters each) limit precision in regime-specific estimates

## 🔬 Policy Implications

1. **Context matters**: Fiscal policy appears more effective during economic downturns (post-2008)
2. **Tax vs. spending**: Tax shocks show larger effects than spending shocks, consistent with Romer & Romer (2010)
3. **Timing is critical**: Multiplier estimates vary significantly across economic regimes

---

**Author:** Khalid Bin Abdullah  
**Institution:** Knox College, Department of Economics  
**Completion Date:** June 2025
