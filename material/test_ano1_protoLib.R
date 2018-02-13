# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- c(timestep = 0.05,
                endTime = 100)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 1e-11),
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
# convert part number from input table to concentration (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f

# Simulate model
start.time <- as.numeric(Sys.time())*1000

output <- sim_ano(input_df, sim_params, model_params)

end.time <- as.numeric(Sys.time())*1000

time.taken <- end.time - start.time
cat(time.taken)

# Plot output
par(mar = c(5,5,2,5))
colnames(output) <- c("time", "calcium", "Cl_ext", "C", "C_c", "C_1", "C_1c", "C_2", "C_2c", "O", "O_c", "O_1", "O_1c", "O_2", "O_2c")
plot(output$time, output$calcium, col="blue", type="l", xlim = c(0,35), xlab="time [s]", ylab="activated Channel [nmol/l]")
lines(output$time, output$C_1, col="red", type="l")
axis(side = 4)
mtext(side = 4, line = 3, 'calcium [a.u]')
legend("topright", legend=c("calcium", "activated channel"),
                   col=c("blue", "red"),
                   lty=c(1,1))
