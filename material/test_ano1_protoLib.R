# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- c(timestep = 1,
                endTime = 10)
# Model Parameters (List)
model_params <- list(vol = 1e-11,
                     init_conc = c(Ca_cyt = 200,
                                   Cl_ext = 30e6,
                                   C = 1,
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
output <- sim_ano(input_df, sim_params, model_params)
output <- as.data.frame(output)

# Plot output
colnames(output) <- c("time", "calcium", "Ca_cyt", "Cl_ext", "C", "C_c", "C_1", "C_1c", "C_2", "C_2c", "O", "O_c", "O_1", "O_1c", "O_2", "O_2c")
plot(output$time, output$calcium, col="blue", xlim=c(90, 160), ylim = c(0, 55), type="l", xlab="time", ylab="concentration")
lines(output$time, output$Ca_cyt, col="black", type = "l")
lines(output$time, output$Cl_ext, col="red", type="l")
lines(output$time, output$C, col="green", type="l")
lines(output$time, output$C_c, col="cyan", type="l")
lines(output$time, output$C_1, col="orange", type="l")
legend("topright", legend=c("calcium", "Ca_cyt", "Cl_ext", "C", "C_c", "C_1"),
                   col=c("blue", "black", "red", "green", "cyan", "orange"),
                   lty=c(1,1))
