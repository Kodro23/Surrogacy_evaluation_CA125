# Evaluation of the Surrogacy of CA-125 Dynamics in Women with Newly Diagnosed Advanced Ovarian Cancer: A Meta-analysis

The main objective of this project is to assess whether the trajectory of CA-125 can serve as a surrogate endpoint for progression-free survival (PFS) in ovarian cancer.

## Data

The analysis is based on a subset of the **MAOV meta-analysis**, consisting of 13 randomized clinical trials including at least 60 patients with previously untreated ovarian cancer.

- **Study period:** January 1, 2001 to September 25, 2016
- **Number of randomized clinical trials:** 13
- **Number of patients:** 8,219
- **Number of repeated CA-125 measurements:** 102,809

## Joint Modeling of CA-125 Trajectory and PFS

### Objective

Estimate the treatment effect on both:

- progression-free survival (PFS), and
- the longitudinal trajectory of CA-125.

### Longitudinal Model

Let \(y(t_{ijk})\) denote the \(k^{th}\) CA-125 measurement for patient \(j\) in trial \(i\).

The longitudinal model for log-transformed CA-125 is:

$$
\log\left(y(t_{ijk})\right)
=
(\beta_0 + b_{0ij})
+
\sum_{k=1}^{3}
(\beta_k + b_{kij})\,ns_k(t_{ijk})
+
\sum_{k=1}^{3}
\gamma_k\,ns_k(t_{ijk})T_{ij}
+
\varepsilon_{ijk}.
$$

The treatment indicator is defined as:

$$
T_{ij}
=
\begin{cases}
1, & \text{if patient } j \text{ in trial } i \text{ received the investigational treatment}, \\
0, & \text{otherwise}.
\end{cases}
$$

The patient-specific random effects and residual errors are assumed to follow:

$$
b_{ij} \sim \mathcal{N}(0,D),
\qquad
\varepsilon_{ijk} \sim \mathcal{N}(0,\sigma^2 I).
$$

Here, \(ns_k(t)\) denotes the \(k^{th}\) basis function of a natural cubic spline representation of time.

The coefficients \(\gamma_k\) characterize the treatment effect on the CA-125 trajectory.

### Survival Model

Progression-free survival is modeled using a proportional hazards model:

$$
h_i(t)
=
h_0(t)\exp(\alpha T_i),
$$

where \(\alpha\) represents the treatment effect on PFS.

The baseline hazard is represented using cubic M-spline basis functions:

$$
h_0(t)
=
\sum_{m=1}^{M}\lambda_m B_m(t),
$$

where \(B_m(t)\) denotes the \(m^{th}\) cubic M-spline basis function.

## Surrogacy Evaluation

### Objective

Quantify the association between:

- the treatment effect on PFS, and
- the treatment effect on the CA-125 trajectory.

The surrogacy analysis follows a correlation-based approach.

First, summary measures of the predicted average log(CA-125) trajectories are computed separately for the investigational and control treatment groups using the joint model.

Treatment effects on these trajectory summaries are then derived for each trial.

A trial-level linear regression is fitted with the estimated treatment effect on PFS, expressed as the log hazard ratio:

$$
\log(HR_i),
$$

as the response variable and the corresponding treatment effect on a CA-125 trajectory summary as the explanatory variable.

For a surrogate measure \(S_i\), the trial-level model can be written as:

$$
\log(HR_i)
=
\beta_0
+
\beta_1 S_i
+
\varepsilon_i.
$$

The coefficient of determination,

$$
R^2,
$$

is used to quantify the strength of the association between the treatment effects on the surrogate endpoint and PFS.

Higher values of \(R^2\) indicate a stronger trial-level association and therefore greater potential for the CA-125 trajectory measure to act as a surrogate endpoint for PFS.

## References

1. **Paoletti X, Lewsley LA, Daniele G, et al.**  
   Assessment of progression-free survival as a surrogate end point of overall survival in first-line treatment of ovarian cancer: a systematic review and meta-analysis.  
   *JAMA Network Open*. 2020;3(1):e1918939.  
   https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2758473

2. **Karamouza E, Glasspool RM, Kelly C, et al.**  
   CA-125 early dynamics to predict overall survival in women with newly diagnosed advanced ovarian cancer based on meta-analysis data.  
   *Cancers*. 2023;15(6):1823.  
   https://doi.org/10.3390/cancers15061823

3. **De Witte D, Molenberghs G, Neyens T, et al.**  
   Evaluation of surrogate endpoints for survival outcomes using the surrogate package in R.  
   *Computer Methods and Programs in Biomedicine*. 2026:109530.  
   https://doi.org/10.1016/j.cmpb.2026.109530

4. **Corbaux P, You B, Glasspool RM, et al.**  
   Survival and modelled cancer antigen-125 ELIMination rate constant K score in ovarian cancer patients in first-line before poly (ADP-ribose) polymerase inhibitor era: A Gynaecologic Cancer Intergroup meta-analysis.  
   *European Journal of Cancer*. 2023;191:112966.  
   https://doi.org/10.1016/j.ejca.2023.112966

5. **Burzykowski T, Coart E, Saad ED, et al.**  
   Evaluation of continuous tumor-size-based end points as surrogates for overall survival in randomized clinical trials in metastatic colorectal cancer.  
   *JAMA Network Open*. 2019;2(9):e1911750.  
   https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2751555

6. **Rustand D, van Niekerk J, Krainski ET, Rue H.**  
   *Bayesian Survival, Longitudinal, and Joint Models with INLA*.  
   Chapman & Hall/CRC; 2026.  
   https://doi.org/10.1201/9781003646822