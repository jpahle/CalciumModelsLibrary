# Set rng seed for debugging tests
set.seed(1)

# Simulation Parameters (Vector)
sim_params <- c(timestep = 0.05,
                endTime = 100)
# Model Parameters (List)
model_params <- list(vol = 5e-14,
                     init_conc = c(Prot_inact = 5,
                                   Prot_act = 0))

# Read Ca timeseries
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert Ca concentration from input table to particle number (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*vol
vol <- 5e-14
input_df["Ca"] <- input_df["Ca"]/f

# Simulate model
output <- sim_calmodulin(input_df, sim_params, model_params)
output <- as.data.frame(output)

# Plot output
colnames(output) <- c("time", "calcium", "Prot_inact", "Prot_act")
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), type="l", xlab="time", ylab="concentration")
lines(output$time, output$Prot_act, col="red", type="l")
legend("topright", legend=c("calcium", "Prot_act"), col=c("blue", "red"), lty=c(1,1))
