# Constructing hierarchical time series through clustering: Is there an optimal way for forecasting?

Manuscript and code for the paper "Constructing hierarchical time series through clustering: Is there an optimal way for forecasting?". The code is used to reproduce the results and figures in the paper.


## Simulation

Run the following code to simulate time series and generate reconciliation results and figures:

```shell
Rscript simulation/simu.R
Rscript simulation/simu_summary.R
```

The first line of code will create `simulation.rds` file in `simulation` folder. The second line of code will generate figures `manuscript/figures/hierarchy_rmsse/simulation/*.pdf`, which will be used in the simulation section of the paper.


## Empirical studies

- `data` folder contains the raw data of the U.S. mortality dataset, which is processed using the `data/dataprocessing.R` and generate `mortality/data.rds`. The `mortality/data.rds` contains the summing matrix and bottom-level series of mortality dataset used in our empirical study. Similarly, `tourism/data.rds` contains the summing matrix and bottom-level series of tourism dataset used in our empirical study.
- `R` folder includes all the code to reproduce the results.
- Use the following code to produce base forecasts for the two datasets. For each rolling window, it will generate base forecasts, calculate features of raw time series and in-sample error, and calculate distance matrix used as input of clustering algorithms. It will produce two folders for saving the computation results `mortality` and `tourism`.
The results are saved in `mortality` and `tourism` folder with names `batch_xx.rds`.

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


- Use the following code to evaluate cluster hierarchies and produce Table 3 and 8, Figure 4 and 12. The plots are saved in `manuscript/figures/mortality`.

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


- Use the following code to generate evaluation.

```shell
Rscript R/run_eval.R tourism
Rscript R/run_eval.R mortality
```


- Use the following code to generate figures used in empirical studies Section.

```shell
Rscript R/run_summary_permute.R tourism ets
Rscript R/run_summary_summary.R mortality ets
```