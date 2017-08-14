# Define simulation parameters
sim_params <- c(timestep = 2,
                endTime = 400)
model_params <- list(vol = 5e-14,
                     init_conc = c(Cal_off = 5,
                                   Cal_on = 0))
# Input calcium time series
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# Simulate Model
output <- calmodulin_detSim(input_df,
                            sim_params,
                            model_params)
# Plot output
plot(output$time, output$Ca, col = "blue", type = "l", xlab = "time [s]", ylab = "concentration [nMol/l]")
lines(output$time, output$Cal_off, col = "grey", type = "l")
lines(output$time, output$Cal_on, col = "red", type = "l")
