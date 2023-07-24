---
title: "Genetics of seasonality - Quantitative traits"
author: "Sölvi Rögnvaldsson"
date: "2023-04-01"
output: html_document
---



## Introduction

## Data generating process 
The assumed data generating process for a seasonally affect QT phenotype is the following
$$
E(y_i) = \beta_g x_g + (\beta_a + \beta_{ag}x_g)\sin\left(\frac{2 \pi}{P}t+\phi+\beta_{pg}x_g\right)
$$
where
* $\beta_g$ - Additive genetic effect
* $\beta_a$ - Baseline seasonal amplitude effect
* $\beta_{ag}$ - Genetic-seasonal amplitude effect
* $\phi$ - Baseline seasonal phase angle
* $\beta_{pg}$ - genetic-seasonal phase angle


