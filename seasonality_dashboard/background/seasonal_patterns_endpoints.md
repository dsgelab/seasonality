---
title: "Seasonal pattern of endpoints in FinRegistry"
author: "Sölvi Rögnvaldsson"
date: "2023-01-16"
output:
  html_document: default
---


### Background
The incidence of several diseases follows a seasonal pattern. Most famously, influenza hits in the autumn but winter peaks are also observed, among others, in cardiovascular disease and mood disorders (Altizer et al. 2006; Stewart et al. 2017; Garbazza and Benedetti 2018). Many factors could explain this pattern, availability of healthcare due to holidays, environmental conditions, genetic factors, demographic variables and more. Unwrapping the seasonal dependency in its full complexity and estimating the influence of each factor is therefore challenging. Despite its apparent importance, the genetic influence on the seasonal variation of disease has received little attention to the best of our knowledge. The reason for this could be that it requires high quality data, linking health registry and sequencing data, which until recently has not been available. The purpose of this project is to gain insight into the genetic influence on seasonal variation in disease and improve our understanding of the phenomenon. To achieve this, we leverage data from the Finnish health registry and Finnish genotype data to explore the effect of genetic mutations on the diagnosis time of different diseases. 
  
### Data used and filtering
In this dashboard we aim to describe the seasonality of disease endpoints in the FinRegistry data set. To determine which endpoints to analyze we used the following filtering criteria:

  * Endpoint has at least 1 significant GWAS loci in the FinnGen dataset.
  * The total number of cases diagnosed in summer and winter is sufficiently large in the FinnGen dataset. The cutoff value used here is 200.   
  
Monthly endpoint counts were gathered in FinRegsitry and the extent of seasonality evaluated using a statistical model. The modeling approach is inspired by Madaniyazi et al.

### Model definition
For each disease we have a time series 
$$(x_{y,m})_{y\in\{1998,...,2019\},m \in \{January,...,December\}}$$
where $x_{y,t}$ is the number of first disease endpoints in a particular month in a particular year. The aim is to extract the annual trend and seasonal trend from the time series. Since we are dealing with counts, we assume they come from a modified Poisson distribution
$$
x_{y,m} \sim \text{quasi-poisson}(30\mu,\theta)
$$
Quasi-Poission is used as it can account for over- or underdispersion, $\theta$, in the count data. In particular, the variance is assumed to be proportional to $\theta$
$$
Var(x_{y,m}) = \theta \mu
$$

and thus $\theta >0$ implies overdispersion while $\theta<0$ implies underdispersion. The mean $\mu$ is the daily mean count, which we parametrize in the following way
$$
\log(\mu) = \alpha + s_t(y) + s_s(m)
$$
where $\alpha$ is the intercept term,  $s_t$ is a smooth function corresponding to the annual trend and $s_s$ is a smooth function corresponding to the seasonal trend. The smooths are a linear combination of B-spline basis functions with discrete penalization terms on the spline coefficients (called p-splines). The seasonal spline has the additional constraint that it is cyclical, i.e. splines for end of December beginning of January should continuously coincide. This type of model is called a Generalized additive model (GAM).
### Seasonality metrics
The metrics used to assess the extent of seasonality are the following   
#### Extreme values at peak and trough
One approach is to evaluate the effect of seasonality at the peak and trough of the cyclical spline function $s_s(m)$. Exponentiation of these values gives the multiplicative effect of seasonality on the mean $\mu$. We therefore calculate the values 
$$
x_p = \max_{m \in \{1,...,12\}} exp(s_s(m))
$$
and   
$$
x_t = \min_{m \in \{1,...,12\}} exp(s_s(m))
$$

\end{align*}
where $p$ and $t$ stand for peak and trough, respectively.   

#### Peak-to-Trough ratio
As the name suggests the peak-to-trough ratio is the ratio between $x_p$ and $x_t$, i.e.
$$
PTR=\frac{x_p}{x_t} =\frac{\max_{m \in \{1,...,12\}} exp(s_s(m))}{\min_{m \in \{1,...,12\}} exp(s_s(m))}
$$

#### Attributable fraction
Attributable fraction due to seasonality is interpreted as the fractional reduction of endpoint counts had the true seasonal process been constant at the estimated seasonal trough value. It is calculated as follows
$$
AF=\sum_{m=1}^{12} p_m(1-exp(s_s(m) - log(x_t)))
$$
where $p_m$ stands for the fraction of all individuals that were diagnosed in month $m$. The term $exp(s_s(m) - log(x_t))$ is called relative risk.     

#### Fraction of variance explained

Fraction of variance explained is the fractional reduction of variance of the residuals obtained by adding the seasonal component to the model. Formally:
$$
VE=max\left(0,1-\frac{s^2_{t,s}}{s^2_{t}}\right)
$$
where $s^2_{t,s}$ is the estimated residual variance under a model with both a trend and seasonal component while $s^2_{t,s}$ is the estimated residual variance under a trend only model.   

### Adjusting for median seasonal pattern

When observing the endpoint count data, we noticed that there is a substantial dip in July for many of the endpoints, one that could be explained simply by holidays. This affects the resulting seasonal components in the models, reflecting in that a great majority of endpoint counts showed substantial seasonality:  

<img src="figure/unnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" style="display: block; margin: auto;" />
   
Since we are not interested in capturing the seasonality of healthcare availability, we attempt to adjust for this by using the median seasonal component, $M_s(x)$ from all endpoints as a linear covariate in our GAM. This component looks like this:

<img src="figure/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" style="display: block; margin: auto;" />
   
Note the dip in the summer and the slight dip during christmas. Our adjusted model then models the mean in the Quasi-Poisson model as
$$
log(\mu) = \alpha + \beta M_s(m)  +s_t(y) + s_s(m)
$$
Note that $M_s(m)$ is fixed and not to be estimated.

## References
Altizer S, Dobson A, Hosseini P, Hudson P, Pascual M, and Rohani P. 2006. 'Seasonality and the Dynamics of Infectious Diseases'. *Ecology Letters* 9 (4): 467–84. https://doi.org/10.1111/j.1461-0248.2005.00879.x.   

Madaniyazi L and Tobias A and Kim Y, Chung Y, Armstrong B, Hashizume  M. 2022. 'Assessing seasonality and the role of its potential drivers in environmental epidemiology: a tutorial'. *International Journal of Epidemiology* 51 (5):1677-1686. https://doi.org/10.1093/ije/dyac115.   

Garbazza C, Benedetti F. 2018. 'Genetic Factors Affecting Seasonality, Mood, and the Circadian Clock'. *Frontiers in Endocrinology* 9 (August): 481. https://doi.org/10.3389/fendo.2018.00481.   

Stewart S, Keates AK, Redfern A,McMurray JJV. 2017. 'Seasonal Variations in Cardiovascular Disease'. *Nature Reviews Cardiology* 14 (11): 654–64. https://doi.org/10.1038/nrcardio.2017.76.



