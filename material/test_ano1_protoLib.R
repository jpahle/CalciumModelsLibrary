# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- c(timestep = 1,
                endTime = 100)
# Model Parameters (List)
model_params <- list(vol = 1e-11,
                     init_conc = c(Cl_ext = 300,
                                   C = 100,
                                   C_c = 0,
                                   C_1 = 0,
                                   C_1c = 0,
                                   C_2 = 0,
                                   C_2c = 0,
                                   O = 0,
                                   O_c = 0,
                                   O_1 = 0,
                                   O_1c = 0,
                                   O_2 = 0,
                                   O_2c = 0))


# Read Ca timeseries
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert Ca concentration from input table to particle number (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vol"]]
input_df["Ca"] <- input_df["Ca"]/f

# Simulate model
start.time <- Sys.time()

output <- sim_ano(input_df, sim_params, model_params)
output <- as.data.frame(output)

end.time <- Sys.time()

time.taken <- end.time - start.time
time.taken

# Plot output
colnames(output) <- c("time", "calcium", "Cl_ext", "C", "C_c", "C_1", "C_1c", "C_2", "C_2c", "O", "O_c", "O_1", "O_1c", "O_2", "O_2c")
plot(output$time, output$calcium, col="blue", type="l", xlab="time", ylab="concentration")
lines(output$time, output$C, col="red", type="l")
lines(output$time, output$C_c, col="orange", type="l")
lines(output$time, output$O, col="green", type="l")
lines(output$time, output$O_c, col="cyan", type="l")
legend("topright", legend=c("calcium", "C", "C_c", "O", "O_c"),
                   col=c("blue", "red", "orange", "green", "cyan"),
                   lty=c(1,1))
