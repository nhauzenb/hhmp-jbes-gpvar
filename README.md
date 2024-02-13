### Code package for Hauzenberger, N., F. Huber, M. Marcellino & N. Petz (202x). Gaussian process vector autoregressions and macroeconomic uncertainty, *Journal of Business & Economic Statistics*, just accepted.

**Publication (open access).** link following soon

[**Final working paper version.**](https://www.dropbox.com/scl/fi/07qcpbq3049zg47cqw3ny/HHMP_JBES_GPVAR-finalwp.pdf?rlkey=l5aopxlbeh7kiwhi41ycbwy3m&dl=0)

[**Online Appendix.**](https://www.dropbox.com/scl/fi/w6o64l6k9wpyjjf9y1u5d/HHMP_JBES_GPVAR-appendix.pdf?rlkey=yfrj506lojz11q93yci68k26c&dl=0)

### Data for Monte Carlo exercise.
For the Monte Carlo exercise, we use the dynamic stochastic general equilibrium (DSGE) model proposed in [Basu & Bundick (2017, ECTA)](https://doi.org/10.3982/ECTA13960) to simulate time series of length T=120 and use these to back out the responses of the model economy to unexpected shocks in uncertainty. The uncertainty shock is identified by ordering the stock market volatility (our measure of uncertainty) first, implying an immediate reaction of all real economic quantities in the system following an uncertainty shock. We provide a single realization from the DGP as a .rda file [`BB_realization`](./data/BB_realization.rda).

### Data for the empirical application. 
For the empirical application,  we use the macroeconomic uncertainty measure of [Jurado, Ludvigson, and Ng (2015, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20131193) provided (and regularly updated) on the [web page of Sydney C. Ludvigson](https://www.sydneyludvigson.com/macro-and-financial-uncertainty-indexes}{sydneyludvigson.com/macro-and-financial-uncertainty-indexes). For the different macroeconomic variables, we rely on the popular FRED-QD dataset provided by the *Federal Reserve Bank of St. Louis* and publicly available [here](https://research.stlouisfed.org/econ/mccracken/fred-databases/). Our quarterly sample spans from 1965Q1 to 2019Q4. Sub-section 4.2 in the paper and Table B.1 in the Online Appendix shows the set of variables included for different model sizes. We provide the data as a .rda file [`unc QD`](./data/unc_QD.rda).

### Estimation files: 

1.) [`A simple univariate example with a GP regression`](!uni_GPreg.R): In Sub-section 2.2, we illustrate the GP regression by means of a simple univariate example. We model US GDP growth as a function of the first lag of a macroeconomic uncertainty measure for a sub-sample around the global financial crisis. This stylized example, highlights the role of the kernel and its hyperparameters crucially impacts the posterior estimates of the conditional mean function. The corresponding estimation file replicates Figure 1 of the paper. 

2.) [`Conjugate GP-VAR with SV`](!GPVAR_main.R): Based on the single realization of the Basu & Bundick (2017, ECTA) DSGE model [`BB_realization`](./data/BB_realization.rda), this file allows to estimate a conjugate Gaussian process vector autoregression (GP-VAR) on an equation-by-equation basis. The conjugate GP-VAR can be estimated either with homoskedastic errors (*sv == "homo"*) or with stochastic volatility (*sv == "SV"*). The folder [`gpvar funcs`](./gpvar_funcs/) contains the MCMC samplers for the GP-VAR and a few auxiliary functions:

* [`Direct sampler for conjugate BVAR with classic symmetric Minnesota prior`](./bvar_funcs/conjVARstd_func.R) 
* [`Direct sampler for BVAR with an asymmetric conjugate Minnesota prior`](./bvar_funcs/conjVARasym_func.R)
* [`Direct sampler for conjugate BVAR with a subspace shrinkage prior`](./bvar_funcs/conjVARsub_func.R)
* [`Gibbs sampler for non-conjugate BVAR with global-local shrinkage priors`](./bvar_funcs/nconjVAR_func.R)


Replication codes come without technical support of any kind.
