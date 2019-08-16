<h3 align="center"> EFAST </h3>
<p align="center">
<img src="./img/efa_uncorr_met.png" width="500px"></img><br/><br/>
<a href="https://travis-ci.org/vankesteren/efast"><img src="https://travis-ci.org/vankesteren/efast.svg?branch=master"></img></a>
</p>
<h4 align="center">Exploratory Factor Analysis with Structured Residuals</h4>

<p align="center">
Code accompanying the manuscript. Work in progress!
</p>

```r
remotes::install_github("vankesteren/efast")
library(efast)
simdat <- simulate_efast()
fit_sim <- efast_hemi(
  data   = simdat, 
  M      = 4, 
  lh_idx = 1:17, 
  rh_idx = 18:34
)
efast_loadings(fit_sim)
```

