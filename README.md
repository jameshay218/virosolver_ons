# virosolver_ons
Scripts and readme demonstrating virosolver methods

## Requisites
- Clone this git repo to your PC and change the `HOME_WD` variable in the `main.R` script.
- You will need to install my `lazymcmc` R package, which is used for the MCMC procedure. This is easy to do with `devtools::install_github("jameshay218/lazymcmc")`.
- You will also need to clone the `virosolver` package from my github. It's currently private, so send me your git user name and I will add you.

## Scripts
Simply run `main.R`, making sure that the file paths at the top of the script are correct. You can change any of the flags to eg. rerun the MCMC chains or not, use only detectable Ct values (change the `use_pos` variable), change the viral kinetics assumptions etc. This script will source a number of smaller scripts, just split up for readibility. Just step through each line to see what is going on.

## Viral kinetics parameters
Assumptions about viral/Ct kinetics are based on our data, so may not be the same for yours. For example, we assume a limit of detection at a Ct of 40 (the "intercept" parameter). You can change these parameters in the appropriate `partab` file in the `pars` directory. 

The `partab` structure is used throughout the analyses and by the MCMC framework. This specifies the parameter bounds, whether the parameter is estimated or not (if fixed=0, then we estimate it. If fixed=1, then it is, unsurprisingly, kept fixed throughout).

## Outputs
The code saves the MCMC chains to the `mcmc_chains` folder. Some plots are also automatically generated in the `plots` folder. Otherwise, some output plots are just generated as R objects in this R session.

The things to look at will be the incidence and prevalence from the simulated line list data, the distribution of simulated Ct values over time, the cross-sectional growth rate estimates from the exp model (matching the preprint), and the full inferred incidence curve from the Gaussian process model.

## MCMC
It is really important to check the chain convergence. Look for the lines that call the function `gelman_diagnostics`. I have run short chains to provide examples, but if you want to actually use the results the chains should be run for longer. This can be changed in the `mcmcPars` variables in the fitting scripts. Increase the `interations` and `adaptive_period` values.
