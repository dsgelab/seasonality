---
title: "Seasonality of clinical endpoints in FinRegistry"
author: "Sölvi Rögnvaldsson"
date: "May 16th 2023"
output:
  html_document: default
---


## Introduction and data
This dashboard presents the results from  a seasonality study analyzing the seasonal patterns in 1,759 clinical endpoints in the Finnish population. Longitudinal data from the FinRegistry project was utilized to collect monthly counts of the first diagnosis of each endpoint between 1998 and 2019 ([https://www.finregistry.fi/](https://www.finregistry.fi/)). Some filtering was done on the data. 

  * Endpoint was required to be recorded in both FinRegistry and FinnGen study
  * At least one diagnosis in each of the years in FinRegistry
  * Endpoint related to external causes of morbidity and uncategorized endpoints were excluded (as determined by the FinnGen endpoint classification at ([https://www.finngen.fi/en/researchers/clinical-endpoints](https://www.finngen.fi/en/researchers/clinical-endpoints))

## Model definition
The modeling approach is inspired by Madaniyazi et al. (2022). For each disease we have a time series 
$$(y_{k,m})_{k\in\{1998,...,2019\},\ m \in \{January,...,December\}}$$
where $y_{k,m}$ is the number of first disease endpoints in a particular month in a particular year. The aim is to extract the annual trend and seasonal trend from the time series. Since we are dealing with counts, we assume they come from a overdispersed Poisson model family
$$
x_{y,m} \sim \text{quasi-poisson}(\mu,\theta)
$$
Quasi-Poission is used as it can account for over- or underdispersion in the count data. In particular, the variance is assumed to be proportional to $\mu$
$$
\text{Var}(x_{y,m}) = \theta \mu
$$

and thus $\theta >1$ implies overdispersion while $\theta<1$ implies underdispersion. The mean $\mu$ is the daily mean count, which we parametrize in the following way
$$
\log(\mu_{k,m}) = log(N_m)+\alpha + s_a(k) + s_s(m)
$$
where $N_m$ is the total number of days in month $m$, $\alpha$ is the intercept term,  $s_a$ is a smooth function corresponding to the annual trend and $s_s$ is a smooth function corresponding to the seasonal trend. The smooth terms are a linear combination of B-spline basis functions with discrete penalization terms on the spline coefficients (called p-splines). The seasonal spline has the additional constraint that it is cyclical, i.e. splines for end of December beginning of January should continuously coincide. This type of model is called a Generalized additive model (GAM).   

The determine whether the seasonal smooth term is significantly different from the zero function, a joint hypothesis test on the spline coefficients presented by Wood et al. (2014) was used. The peak-to-trough ratio (PTR) is a useful metric to characterize the extent of seasonality. It is the ratio between the peak and trough value of the seasonal smooth term, namely

$$
\text{PTR}=\frac{\max_{x \in [1,13)} \exp(s_s(m))}{\min_{x \in [1,13)} \exp(s_s(m))}
$$
To obtain bounds on PTR as well as the the peak time and trough time, a Monte Carlo simulation was performed. In each iteration of the simulation, spline coefficients were sampled from their posterior distribution, and the PTR and peak and trough times calculated from the induced seasonal smooth term. The 2.5\% and 97.5\% quantiles from the simulation were then taken to be 95\% posterior intervals for the different metrics 

## Adjusting for seasonal behaviour
When extracting the seasonal patterns for non-acute endpoints, it cannot be avoided that certain seasonal behavioral pattern will be captured in the process. Namely, that people in Finland tend to take vacations in the summer and spend time with family during the Christmas period, which can impact hospital activity and hence leading to fewer disease diagnosis of non-acute endpoints. Since this type of seasonality is not necessarily of interest, an attempt was made to adjust for such seasonal behavior when fitting the GAM. In the same spirit as Madaniyazi et al. (2022), three approaches were developed and will briefly be explained below.

### Binary covariates for July and December
The first approach simply adds two binary covariates to the mean structure, indicating whether month $m$ is July or December, respectively. Namely, functions $f_j$ and $f_d$ were defined such that
$$
    f_j(m) = \begin{cases} 1 & \text{if m=July} \\
                           0 & \text{otherwise} 
             \end{cases}
$$

$$
    f_d(m) = \begin{cases} 1 & \text{if m=December} \\
                           0 & \text{otherwise} 
             \end{cases}
$$
The two covariates were then added to the mean structure of the GAM
$$
        \log(\mu_{k,m}) = \log(N_m) + \alpha + s_a(k) + s_s(m) + \gamma_j f_j(m) + \gamma_d f_d(m)
$$
where $\gamma_j$ and $\gamma_d$ are parameters specifically accounting for counts in July and December, respectively.  

### Normalizing by mean seasonal curve
The second approach was to fit an unadjusted GAM on the whole set of endpoints and then assume that the mean seasonal curve reflects seasonal behavior. Formally, the mean seasonal smooth  term $\bar{s}_{s}$ was calculated in the following way
$$
   \bar{s}_s(x) = \frac{1}{K}\sum_{i=1}^K s_{s,i}(x) \label{eq:avg_seasonal_comp}
$$
where $K$ is the total number of endpoints considered and $s_{s,i}(x)$ is the seasonal smooth term for endpoint $i$. The mean seasonal smooth term at each month was then added as an offset to the mean structure in the GAM to act as a normalizing factor
$$
        \log(\mu_{k,m}) = \log(N_m) + \bar{s}_s(m) + \alpha + s_a(k) + s_s(m) 
$$
   
### Mean seasonal curve as a covariate
The third approach involves adjusting for the mean seasonal curve as a linear covariate in the mean structure of the GAM, allowing the normalization to be vary between different endpoints
$$
        \log(\mu_{k,m}) = \log(N_m) + \alpha + s_a(k) + s_s(m) + \gamma\bar{s}_s(m)
$$
where $\gamma$ is parameter reflecting the extent to which the mean seasonal component is accounted for. 

## References

Madaniyazi L and Tobias A and Kim Y, Chung Y, Armstrong B, Hashizume  M. (2022). "Assessing seasonality and the role of its potential drivers in environmental epidemiology: a tutorial". *International Journal of Epidemiology* 51.5, pp. 1677-1686. doi: 10.1093/ije/dyac115.  

Wood, S. N. (2013). "On p-values for smooth components of an extended generalized additive model". *Biometrika* 100.1, pp. 221–228. doi: 10.1093/biomet/ass048.




