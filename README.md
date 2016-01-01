This repository contains the data and source code for the following paper:

* J. Urbano, "[Test Collection Reliability: A Study of Bias and Robustness to Statistical Assumptions via Stochastic Simulation](http://julian-urbano.info/files/publications/064-test-collection-reliability-study-bias-robustness-statistical-assumptions-stochastic-simulation.pdf)", *Information Retrieval Journal*, 2015.

A [single ZIP file](https://github.com/julian-urbano/irj2015-reliability/archive/master.zip) can be downloaded as well.


## Project Structure

* `bin/` Shell scripts to run the code in Linux.
* `config/` Machine-dependent configuration files.
* `data/` Input data files.
* `output/` Generated output files.
* `src/` Source code in R.
* `scratch/` Temporary files generated in the process.

## How to replicate the results in the paper 

It takes several days to run all the code, so it is ready to use an SGE cluster, parallelizing across collections, measures and assumptions. You can still run everything on a single machine using just one core. **It is important that you always run from the base directory**.

1. Edit file `config/environment.sh`:
    * If you want to use a cluster, set variable `SGE=true`. Edit file `bin/qsub.sub` as well to change the notification e-mail and make sure that R is properly loaded in the SGE job.
    * If you don't want to use a cluster, set variable `SGE=false` and make the variable `RSCRIPT` point to the correct path in your machine.
2. Run script `bin/01-simulation.sh`. This simulates the new collections and stores them in `scratch/01-simulation/`.
3. Run script `bin/02-diagnosis.sh`. This computes all diagnosis indicators and stores them in `scratch/02-diagnosis/`.
4. Run script `bin/03-diagnosis-compile.sh`. This compiles all diagnosis data, runs the variance components analysis, and stores the results in `output/diagnosis/`.
5. Run script `bin/11-measures.sh`. This computes all accuracy and reliability scores, and stores them in `scratch/11-measures/`.
6. Run script `bin/12-measures-compile.sh`. This compiles all accuracy and reliability scores, runs the variance components analysis, and stores the results in `output/measures/`.
7. Run script `bin/99-paper.sh`. This generates all figures and tables in the paper and stores them in `output/paper/`.

## How to customize and extend the results in the paper

You can easily customize and extend the code to run with your own initial test collections, compute your own diagnosis indicators or your own measures of test collection accuracy. 

If you only want to use certain collections or measures, edit file `config/params.sh` and follow the instructions. If you want to analyze different topic set sizes or use a different number of trials, edit file `config/params.R`.

Note that the script `src/99-paper.R` is only intended to generate the figures and tables in the paper. If you customize something and want a similar analysis, you will need to extend this script yourself.

### Custom test collections

Simply add new CSV files with the topic-by-system effectiveness matrices in directory `data/` (no row names). Take a look for instance at the file from the [Robust 2003 collection](/data/robust2003.csv). After adding your own files, run all the scripts again as detailed [above](#how-to-replicate-the-results-in-the-paper).

### Custom diagnosis indicators

You can modify function `diagnose.simulation` in file `src/simulation.R` to customize the diagnosis indicators. This function is defined as:

    diagnose.simulation <- function(X, Y,
                                    decX = decompose.effects(X),
                                    F_tX = estimate.cdf(decX$NU_t, xmin = -1, xmax = 1),
                                    F_EX = apply(decX$E, 2, function(e) estimate.cdf(e, xmin = min(e), xmax = max(e))),
                                    gX = g.study(X))

Parameters `X` and `Y` are the original and simulated collections, respectively. The rest of parameters are passed for efficiency purposes. They are the effect decomposition of the original collection, its kernel-density estimated marginal distributions, and its G-study. Take a look at `src/simulation.R` for the details, especially functions `decompose.effects` and `configure.simulation`.

This function is called from `src/02-diagnosis.R` for each simulated collection. You can add a new indicator by adding its name and value in the returned object. After that, run script `bin/02-diagnosis.sh` again to compute all diagnosis scores, and `bin/03-diagnosis-compile.sh` to aggregate results and run the variance components analysis.

### Custom accuracy/reliability measures

You can add new measures by creating a file `src/measures/measure.<name>.R` for each of them. This file must contain two functions:

* `measure.<name>.estimate`, which receives a simulated collection `X` and new topic set sizes `n_t`. It must return the *expected accuracy* of new collections of those sizes. This corresponds to *g(X, n_t')* in the paper.
* `measure.<name>.actual`, which receives a simulated collection `X_estimate` and its original collection `X_truth`. It must return the *actual accuracy* of the simulated collection. This corresponds to *A(X, mu)* in the paper.

Take a look for instance at measure [`EtauAP`](/src/measures/measure.EtauAP.R). These functions are called from `src/11-measures.R` for each simulated collection. After adding your own files, run script `bin/11-measures.sh` again to compute all estimates, and `bin/12-measures-compile.sh` to aggregate results and run the variance components analysis.

## How to simulate new arbitrary collections

You can easily simulate new, arbitrary collections following these steps. First, load your original collection `X`. It must be a matrix or data frame where columns are systems and rows are topics:

    X <- read.csv("collection.csv")
    head(X)
        sys1   sys2   sys3 ...
    1 0.3878 0.3685 0.3951 ...
    2 0.2000 0.0667 0.0092 ...
    3 0.3290 0.3455 0.0472 ...
    ...

Load the script `src/simulation.R` (it also loads `src/gt4ireval.R` and `src/measures/tauAP.R` for diagnosis; you can ignore them if you just want to simulate):

    source("src/simulation.R")

Configure the simulation with your original data and desired statistical assumptions:

    cfg <- configure.simulation(X,
                                normal = FALSE, homoscedastic = FALSE,
                                uncorrelated = FALSE, random = TRUE)

Simulate a new collection `Y` based on the configuration and with the number of topics of your choice:
 
    Y <- simulate.collection(cfg, 100)

If you want, you can also diagnose the simulated collection. If you just want to diagnose a few collections, you can directly run:

    diagnosis <- diagnose.simulation(X, Y)

If you plan on diagnosing many collections, for efficiency it is advised to first precompute some things about the original collection and pass them to the diagnosis function:

    decX <- decompose.effects(X)
    gX <- g.study(X)
    F_tX <- estimate.cdf(decX$NU_t, xmin = -1, xmax = 1)
    F_EX <- apply(decX$E, 2, function(e)
      estimate.cdf(e, xmin = min(-1, min(e)), xmax = max(1, max(e))))

    diagnosis <- diagnose.simulation(X, Y, decX, F_tX, F_EX, gX)

## License

* The TREC results in `data/` are anonymized and posted here with permission from the organizers.
* Databases and their contents are distributed under the terms of the [Creative Commons Attribution-ShareAlike 4.0 International License](http://creativecommons.org/licenses/by-sa/4.0/).
* Software is distributed under the terms of the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0-standalone.html).

When using this archive, please [cite](CITE.bib) the above paper:

    @article{Urbano2015:reliability,
      author = {Urbano, Juli\'{a}n},
      journal = {Information Retrieval Journal},
      title = {{Test Collection Reliability: A Study of Bias and Robustness to Statistical Assumptions via Stochastic Simulation}},
      year = {2015}
    }
