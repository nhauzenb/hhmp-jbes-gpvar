### Code package: N. Hauzenberger, F. Huber, M. Marcellino & N. Petz (202x). Gaussian process vector autoregressions and macroeconomic uncertainty, *Journal of Business & Economic Statistics*, forthcoming.

[**Publication (open access).**](https://doi.org/10.1080/07350015.2024.2322089)

[**Final working paper version.**](https://www.dropbox.com/scl/fi/j43tbl3pxpd7m40xerwhg/HHMP_JBES_GPVAR-finalwp.pdf?rlkey=vv9nevscnf4dvlvojcwucor6n&dl=0)

[**Online Appendix.**](https://www.dropbox.com/scl/fi/w6o64l6k9wpyjjf9y1u5d/HHMP_JBES_GPVAR-appendix.pdf?rlkey=yfrj506lojz11q93yci68k26c&dl=0)

### Data for Monte Carlo exercise.
For the Monte Carlo exercise, we use the dynamic stochastic general equilibrium (DSGE) model proposed in [Basu & Bundick (2017, ECTA)](https://doi.org/10.3982/ECTA13960) to simulate time series of length T=120 and use these to back out the responses of the model economy to unexpected shocks in uncertainty. The uncertainty shock is identified by ordering the stock market volatility (our measure of uncertainty) first, implying an immediate reaction of all real economic quantities in the system following an uncertainty shock. We provide a single realization from the DGP as a .csv file [`BB_realization`](./data/BB_realization.csv).

### Data for the empirical application. 
For the empirical application,  we use the macroeconomic uncertainty measure of [Jurado, Ludvigson, and Ng (2015, AER)](https://www.aeaweb.org/articles?id=10.1257/aer.20131193) provided (and regularly updated) on the [web page of Sydney C. Ludvigson](https://www.sydneyludvigson.com/macro-and-financial-uncertainty-indexes). For the different macroeconomic variables, we rely on the popular FRED-QD dataset provided by the *Federal Reserve Bank of St. Louis* and publicly available [here](https://research.stlouisfed.org/econ/mccracken/fred-databases/). Our quarterly sample spans from 1965Q1 to 2019Q4. Sub-section 4.2 in the paper and Table B.1 in the Online Appendix shows the set of variables included for different model sizes. We provide the data as a .rda file [`MacroUncQ`](./data/MacroUncQ.rda).

### Estimation files: 

**1.) [`A simple univariate example with a GP regression`](!univariate-GPreg.R):** In Sub-section 2.2, we illustrate the GP regression by means of a simple univariate example. We model US GDP growth as a function of the first lag of a macroeconomic uncertainty measure for a sub-sample around the global financial crisis. This stylized example, highlights the role of the kernel and its hyperparameters crucially impacts the posterior estimates of the conditional mean function. The corresponding estimation file replicates Figure 1 of the paper. 

**2.) Conjugate GP-VAR with SV:** Based on the single realization of the Basu & Bundick (2017, ECTA) DSGE model [`BB_realization`](./data/BB_realization.csv), the file [`GPVAR_eqbyeq`](!GPVAR_eqbyeq.R) allows to estimate of a conjugate Gaussian process vector autoregression (GP-VAR) on an equation-by-equation basis. The [`GPVAR_girfs`](!GPVAR_girfs.R) collects the equation-wise estimates and allows to compute generalized impulse response functions (GIRFs). In addition, the folder [`gpvar funcs`](./gpvar_funcs/) contains the MCMC sampler for a GP regression and a C++ function to create a squared exponential kernel:

* [`MCMC sampler for conjugate GP-VAR with SV`](./gpvar_funcs/gp_eqbyeq_mcmc.R) 
* [`Squared exponential kernel function`](./gpvar_funcs/sqexp_kernel.cpp)


Replication codes come without technical support of any kind.
