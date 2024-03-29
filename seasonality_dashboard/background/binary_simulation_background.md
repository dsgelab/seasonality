---
title: "Genetics of seasonality - Disease endpoints"
author: "Sölvi Rögnvaldsson"
date: "2023-03-17"
output: html_document
---




## Introduction
In order to get a better understanding of the results from the seasonal GWAS, we aim to simulate the process generating our data. This allows us to carry out experiments using different assumptions to evaluate statistical power in a GWAS setting for different seasonal phenotypes (binary vs. quantitative) as well as assessing the robustness of associations towards the lag between disease onset and disease diagnosis. 

## Data generating process
Let $T$ be a random variable representing the time of year (in months) of disease onset for a particular individual, $T\in [0,12)$. Further let $g$ be a genetic variant and the genotype of $g$ given by $x_g$. We propose a periodic density function with the parametric form
$$
f(t) =  \frac{1}{P} + \frac{1}{P}(\beta_b+\beta_gx_g) sin\left(\frac{2\pi}{P}t + \phi\right)
$$
where $P=12$ is the period, $\beta_b$ is the baseline seasonal component, $\beta_g$ is the interaction effect between the the genotype and the season and $\phi$ is the phase angle of the seasonal fluctuation. Notice that if $\beta_b=0$ and $\beta_g=0$ then we have a uniform distribution, reflecting no seasonality. $f(t)$ is a density if it integrates to 1 and $f(t) \geq 0$ for all $t$. The first condition can be quickly shown to hold true

$$
\int_0^P f(t) = \int_0^P \left ( \frac{1}{P} + \frac{1}{P}(\beta_b+\beta_gx_g) sin\left(\frac{2\pi}{P}t + \phi\right) \right) dt 
= \left [ \frac{t}{P} - \frac{1}{2 \pi}(\beta_b+\beta_gx_g) cos\left(\frac{2\pi}{P}t + \phi\right)\right]_0^P 
= 1 - \frac{1}{2 \pi}(\beta_b+\beta_gx_g) cos(\phi) + \frac{1}{2 \pi}(\beta_b+\beta_gx_g) cos(\phi) 
= 1
$$


However, the second condition of being non-negative is only true under certain conditions on $\beta_b$ and $\beta_g$, namely
$$
-1 \leq \beta_b + x\beta_g \leq 1
$$
for all $x \in [0,2]$. Note that because of linearity, it is sufficient that the condition holds for $x=0$ and $x=2$

## Modeling diagnosis lag
The random variable $T$ is the true disease onset time but in reality we observe a diagnosis time $T'$ which might differ from T substantially, due to many factors s.a. onset of symptoms, holidays, personality of the patient, data anonymization etc. In this simulation study, we assume
$$
T'=T+T_\ell+\frac{T_a}{28}
$$
where $T_\ell$ is the lag time and $T_a$ is the data anonymization perturbation (in days) as performed in the FinnGen data. We let
$$
T_\ell \sim \text{gamma}(\mu,\sigma=\sqrt{\mu})
$$
and the FinnGen anonymization is
$$
T_a \sim Unf\{-15,-14,...,-1,1,...,15\}
$$
