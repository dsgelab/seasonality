---
title: "Modeling the underlying process of seasonality with genetic seasonal effect"
author: "Sölvi Rögnvaldsson"
date: "May 16th"
output: html_document
---




## Introduction
In order to get a better understanding of how useful the seasonality phenotypes are to capture true genetic seasonal effects when used in GWAS, a process that can generate realistic data sets containing such genetic seasonal effects is needed. Simulating from such a process can help assessing various statistical properties of the seasonality GWAS framework. This dashboad presents an interactive user interface to simulated data from an individual-level time-of-onset density introduced below.

## Data generating process
Let $T$ be a random variable representing the time of year (in months) of disease onset for a particular individual, $T\in [0,12)$. Further let $g\in \{0,1,2\}$ be the genotype of a particular genetic variant. A simple periodic density function is proposed with the parametric form
$$
    f(t) =  \frac{1}{P} + \frac{1}{P}(\beta_b+\beta_{g}g) \sin\left(\frac{2\pi}{P}t + \phi\right)
$$
where $P=12$ is the period, $\beta_b$ is the baseline seasonal effect, $\beta_{g}$ is the genetic seasonal effect and $\phi$ is the phase angle of the seasonal fluctuation. Notice that if $\beta_b=0$ and $\beta_{g}=0$ then $f(t)$ is a the density of the uniform distribution, reflecting no seasonality. To validate that $f(t)$ is truly a density, it needs to be checked whether it integrates to 1 and $f(t) \geq 0$ for all $t$. The first condition can quickly be shown to hold true

$$
\int_0^P f(t) = \int_0^P \left ( \frac{1}{P} + \frac{1}{P}(\beta_b+\beta_gx_g) sin\left(\frac{2\pi}{P}t + \phi\right) \right) dt 
= \left [ \frac{t}{P} - \frac{1}{2 \pi}(\beta_b+\beta_gx_g) cos\left(\frac{2\pi}{P}t + \phi\right)\right]_0^P 
= 1
$$

using the fact that $cos(\phi+2\pi)=cos(\phi)$. However, the second condition of being non-negative is only true under certain conditions on $\beta_b$ and $\beta_{g}$, namely
$$
-1 \leq \beta_b + x\beta_g \leq 1
$$
for all $x \in [0,2]$. Note that because of linearity, it is sufficient that the condition holds for $x=0$ and $x=2$

## Modeling diagnosis lag
The density in Equation \ref{eq:seasonal_dens} models the true disease onset time $T$ but in reality a diagnosis time $T'$ is observed which might differ from $T$ substantially, due to many factors such as slow progression of symptoms, holidays and data anonymization. The diagnosis lag therefore highly depends on the disease endpoint in question. When analyzing seasonality in chronic diseases, the diagnosis lag can be quite substantial due to waiting times for an appointment with a specialist or personal decisions to seek help after the holiday \citep{wang_delays_2004}. In the case of more acute endpoints, such as stroke, myocardial infarction or infectious diseases, the diagnosis times are more reflective of the true disease onset times. To model the diagnosis lag in a simulated setting, we assumed
$$
T'=T+T_\ell+\frac{T_a}{28}
$$
where $T_\ell$ is the lag time (in months) and $T_a$ is a data anonymization perturbation (in days) as performed in the FinnGen data. We further assume a gamma distribution for the diagnosis lag
$$
T_\ell \sim \text{gamma}(\mu,\sigma=\sqrt{\mu}+1)
$$
The FinnGen anonymization is a perturbation by a number drawn from the distribution
$$
T_a \sim \text{uniform}\{-15,-14,...,-1,1,...,15\}
$$
