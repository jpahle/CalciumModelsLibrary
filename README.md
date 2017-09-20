
<!-- README.md is generated from README.Rmd. Please edit that file -->
CalciumModelsLibrary
====================

The goal of CalciumModelsLibrary is to provide the user with functions to study the differential activation of calcium-sensitive proteins, such as calmodulin, CaM Kinase II, glycogen phosphorylase and others.

Installation
------------

You can install CalciumModelsLibrary from github with:

``` r
# install.packages("devtools")
devtools::install_github("jpahle/CalciumModelsLibrary")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(CalciumModelsLibrary)

sim_params <- c(timestep = 0.05,
                endTime = 100)
model_params <- list(vols      = c(vol = 5e-14),
                     init_conc = c(Prot_inact = 5,
                                   Prot_act = 0))
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f
output <- sim_calmodulin(input_df, sim_params, model_params)
output <- as.data.frame(output)
colnames(output) <- c("time", "calcium", "Prot_inact", "Prot_act")
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), ylim=c(0,15), type="l", xlab="time", ylab="concentration")
lines(output$time, output$Prot_act, col="red", type="l")
legend("topright", legend=c("calcium", "Prot_act"), col=c("blue", "red"), lty=c(1,1))
```

![](README-example-1.png)
