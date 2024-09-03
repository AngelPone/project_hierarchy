# Constructing hierarchical time series through clustering: Is there an optimal way for forecasting?

This repository supports the publication "Constructing hierarchical time series through clustering: Is there an optimal way for forecasting?" written by Bohan Zhang, Anastasios Panagiotelis and Han Li. In this paper, xtend existing work that uses time series clustering to construct hierarchies, with the goal of improving forecast accuracy, in threeways. First, we investigate multiple approaches to clustering, including not only different clustering algorithms, but also the way time series are represented and how distance between time series is defined. We find that cluster-based hierarchies lead to improvements in forecast
accuracy relative to two-level hierarchies. Second, we devise an approach based on random permutation of hierarchies, keeping the structure of the hierarchy fixed, while time series are randomly allocated to clusters. In doing so, we find that improvements in forecast accuracy that accrue from using clustering do not arise from grouping together similar series but from the structure of the hierarchy. Third, we propose an approach based on averaging forecasts across hierarchies constructed using different clustering methods, that is shown to outperform any single clustering method. All analysis is carried out on two benchmark datasets and a
simulated dataset. Our findings provide new insights into the role of hierarchy construction in forecast reconciliation and offer valuable guidance on forecasting practice.

This repository contains the mortality dataset, tourism dataset and necessary code to produce tables and figures in the paper.





## Overview

- `data` folder contains the raw data of the U.S. mortality dataset, which is processed using the `data/dataprocessing.R` and generate `mortality/data.rds`. The `mortality/data.rds` contains the summing matrix and bottom-level series of mortality dataset used in our empirical study. Similarly, `tourism/data.rds` contains the summing matrix and bottom-level series of tourism dataset used in our empirical study.
- `R` folder includes all the code to reproduce the results. 
- There are several steps to produce the results in the paper. Note that due to large size, the intermediary results are not included in this repository. You can download them from this link.
    - produce base forecasts for all evaluation windows
    - produce clusters based on base forecasts 
    - evaluate performance of benchmarks, cluster hierarchies, and combination hierarchies, resulting in Table 3, 7 Figure 4, 12
    - produce, forecast and evaluate twin hierarchies, resulting in Figure 6, 7, 13
    - use other code to produce other figures and tables.
- Most scripts enable parallel computing using the `future` and `furrr` R packages. You can set number of cores in line 8 of `R/utils.R`. 

## Empirical studies

- Use the following code to produce base forecasts for the two datasets. For each rolling window, it will generate base forecasts, calculate features of raw time series and in-sample error, and calculate distance matrix used as input of clustering algorithms.  The results are saved in `mortality` and `tourism` folder with names `batch_xx.rds`.

```shell
Rscript R/run_base.R tourism
Rscript R/run_base.R mortality
```

- Use the following code to generate cluster hierarchies.

```shell
Rscript R/run_nl.R tourism
Rscript R/run_nl.R mortality
```


- Use the following code to generate reconciled forecasts for cluster hierarchies and combination hierarchies (combination and grouped).

```shell
Rscript R/run_nlf.R tourism
Rscript R/run_nlf.R mortality
Rscript R/run_comb.R
```


- Use the following code to evaluate cluster hierarchies and produce Table 3 and 8, Figure 4 and 12. The plots are saved in `manuscript/figures/mortality` and `manuscript/figures/tourism`.

```shell
mkdir -p manuscript/figures/mortality
mkdir -p manuscript/figures/tourism
Rscript R/run_eval_cluster.R mortality
Rscript R/run_eval_cluster.R tourism
```

- Use the following code to generate twin hierarchies for natural and best cluster hierarchy, and produce reconciled forecasts for twin hierarchies.

```shell
Rscript R/run_permute.R tourism
Rscript R/run_permute.R mortality

Rscript R/run_nlf.R mortality
Rscript R/run_nlf.R tourism
```

- Use the following code to generate evaluation and figures used in twin hierarchies vs natural/best cluster hierarchy.

```shell
Rscript R/run_summary_permute.R tourism
Rscript R/run_summary_summary.R mortality
```



- Use the following code to generate other tables in the manuscript.

```
Rscript R/features.R
```

## Simulation

Run the following code to simulate time series and generate reconciliation results and figures:

```shell
Rscript simulation/simu.R
Rscript simulation/simu_summary.R
```

The first line of code will create `simulation.rds` file in `simulation` folder. The second line of code will generate figures `manuscript/figures/hierarchy_rmsse/simulation/*.pdf`, which will be used in the simulation section of the paper.




