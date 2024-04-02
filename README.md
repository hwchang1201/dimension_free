# Simulation code for "Dimension-free Relaxation Times of Informed MCMC Samplers on Discrete Spaces".

## Quick start

### Bayesian variable selection simulation

```{r}
Rscript bvs.R cov method clip num_var
```
* cov: the dependent structure of the design matrix ("high", "moderate")
* method: the MCMC sampler used ("single", "informed")
* clip: the clipping threshold for the informed sampler ("yes", "no")
* num_var: the number of variables in the initial model (integer within 0 to 200)


### Bayesian variable selection simulation

```{r}
Rscript sbm.R method initial
```

* method: the MCMC sampler used ("single", "informed")
* initial: initialization scheme in the main text ("good", "bad")
