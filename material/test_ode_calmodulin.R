# Define simulation parameters
sim_params <- c(timestep = 2,
                endTime = 400)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 5e-14),
                     init_conc = c(Cal_off = 5,
                                   Cal_on = 0))
# Input calcium time series
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert Ca concentration from input table to particle number (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f
# Simulate Model
output <- calmodulin_detSim(input_df,
                            sim_params,
                            model_params)
# Plot output
plot(output$time, output$Ca, col = "blue", type = "l", xlab = "time [s]", ylab = "concentration [nMol/l]")
lines(output$time, output$Cal_on, col = "red", type = "l")
legend("topright", legend=c("calcium", "Cal_act"), col=c("blue", "red"), lty=c(1,1))
