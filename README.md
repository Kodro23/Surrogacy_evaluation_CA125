# Evaluation of the Surrogacy of CA-125 Dynamics in Women with Newly Diagnosed Advanced Ovarian Cancer: A Meta-analysis

The main objective of this project is to assess whether the trajectory of CA-125 can serve as a surrogate endpoint for progression-free survival (PFS) in ovarian cancer.

## Data

Data is a subset of the **MAOV meta-analysis**, consisting of 13 randomized clinical trials including 8,219 patients with previously untreated ovarian cancer (47,463 observations), included from January 1, 2001 to September 25, 2016.

## Joint Modeling of CA-125 Trajectory and PFS

To estimate the treatment effect on both PFS, and the longitudinal trajectory of CA-125.

**Longitudinal Model**

For $y(t_{ijk})$, the $k^{th}$ value of biomarker CA-125 for patient $j$, in trial $i$:

$$
\log(y(t_{ijk})) = (\beta_0 + b_{0ij}) + \sum_{k=1}^{3} (\beta_k + b_{k_{ij}}) \, ns(t_{ijk}) + \sum_{k=1}^{3} \gamma_k \ ns(t_{ijk}) T_{ij}+ \varepsilon_{ijk}
$$

where

$$
T_{ij} =
\begin{cases}
1 & \text{if patient } j \text{ in trial } i \text{ received the investigational treatment} \\
0 & \text{otherwise}
\end{cases}
$$

and

$$
b_{ij} \sim \mathcal{N}(0, D),
\qquad
\varepsilon_{ijk} \sim \mathcal{N}(0, \sigma^2 I)
$$

**Survival model**

$$
h_i(t) = h_0(t)\exp(\alpha \, T_i)
$$

with

$$
h_0(t) = \sum_{m=1}^{M} \lambda_m B_m(t)
$$

where $B_m(t)$ is a cubic M-spline basis function.

**Joint association structure**

Henderson's shared frailty model where the survival and longitudinal models share common random effects:

$$
h_i(t \mid b_i) =
h_0(t)\exp\left(\alpha \, T_i + \eta^\top b_i\right)
$$

## Surrogacy Evaluation
To quantify the association between the treatment effect on PFS, and the treatment effect on the CA-125 trajectory.

The surrogacy analysis follows a correlation-based approach: first, summary measures of the predicted average log(CA-125) trajectories are computed, accounting fo treatment group; thenaA trial-level linear regression is fitted with the estimated treatment effect on PFS, expressed as the log hazard ratio. For a surrogate measure $S_i$, the trial-level model can be written as:

$\log(HR_i)=\beta_0+\beta_1 S_i+ \varepsilon_i.
$

The coefficient of determination, $R^2$ is used to quantify the strength of the association between the treatment effects on the surrogate endpoint and PFS.

## References

- Paoletti X, Lewsley LA, Daniele G, et al. Assessment of progression-free survival as a surrogate end point of overall survival in first-line treatment of ovarian cancer: a systematic review and meta-analysis. *JAMA Network Open*. 2020;3(1):e1918939. https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2758473

- Karamouza E, Glasspool RM, Kelly C, et al. CA-125 early dynamics to predict overall survival in women with newly diagnosed advanced ovarian cancer based on meta-analysis data. *Cancers*. 2023;15(6):1823.https://doi.org/10.3390/cancers15061823

- De Witte D, Molenberghs G, Neyens T, et al. Evaluation of surrogate endpoints for survival outcomes using the surrogate package in R. *Computer Methods and Programs in Biomedicine*. 2026:109530. https://doi.org/10.1016/j.cmpb.2026.109530

- Corbaux P, You B, Glasspool RM, et al. Survival and modelled cancer antigen-125 ELIMination rate constant K score in ovarian cancer patients in first-line before poly (ADP-ribose) polymerase inhibitor era: A Gynaecologic Cancer Intergroup meta-analysis. *European Journal of Cancer*. 2023;191:112966. https://doi.org/10.1016/j.ejca.2023.112966

- Burzykowski T, Coart E, Saad ED, et al. Evaluation of continuous tumor-size-based end points as surrogates for overall survival in randomized clinical trials in metastatic colorectal cancer. *JAMA Network Open*. 2019;2(9):e1911750. https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2751555

- Rustand D, van Niekerk J, Krainski ET, Rue H. *Bayesian Survival, Longitudinal, and Joint Models with INLA*. Chapman & Hall/CRC; 2026. https://doi.org/10.1201/9781003646822