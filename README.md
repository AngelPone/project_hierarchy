# On the performance of hierarchy construction: impacts of structure and clustering

Manuscript and code for the paper "On the performance of hierarchy construction: impacts of structure and clustering". The code is used to reproduce the results and figures in the paper.


## Simulation

Run the following code to simulate time series and generate reconciliation results and figures:

```shell
Rscript simulation/simu.R
Rscript simulation/simu_summary.R
```

The first line of code will create `simulation.rds` file in `simulation` folder. The second line of code will generate figures `manuscript/figures/simu*.pdf`, which will be used in the simulation section of the paper.


## Empirical studies

- `data` folder contains the raw data of the U.S. mortality dataset, which is processed using the `data/dataprocessing.R` and generate `mortality/data.rds`. The `mortality/data.rds` contains the summing matrix and bottom-level series of mortality dataset used in our empirical study. Similarly, `tourism/data.rds` contains the summing matrix and bottom-level series of tourism dataset used in our empirical study.
- `R` folder includes all the code to reproduce the results.
- Use the following code to produce base forecasts for the two datasets. For each rolling window, it will generate base forecasts, calculate features of raw time series and in-sample error, and calculate distance matrix used as input of clustering algorithms. It will produce two folders for saving the computation results `mortality/ets` and `tourism/ets`.

```shell
Rscript R/run_base.R tourism ets
Rscript R/run_base.R mortality ets
```

- Use the following code to generate cluster hierarchies.

```shell
Rscript R/run_nl.R tourism ets
Rscript R/run_nl.R mortality ets
```

- Use the following code to generate random hierarchies.

```shell
Rscript R/run_nlr.R tourism ets
Rscript R/run_nlr.R mortality ets
```

- Use the following code to generate reconciled forecasts for all hierarchies.

```shell
Rscript R/run_nlf.R tourism ets
Rscript R/run_nlf.R mortality ets
```

- Use the following code to generate evaluation.

```shell
Rscript R/run_eval.R tourism ets
Rscript R/run_eval.R mortality ets
```


- Use the following code to generate figures used in empirical studies Section.

```shell
Rscript R/summary.R tourism ets
Rscript R/summary.R mortality ets
```