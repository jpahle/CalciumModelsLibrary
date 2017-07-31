# Define simulation parameters
sim_params <- c(timestep = 0.05,
                endTime = 100)
model_params <- list(vol = 5e-14,
                     init_conc = c(A = 5,
                                   B = 0))

# Input calcium time series
# input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))

# Simulate Model
output <- detSim(exampleModel,
                 sim_params,
                 model_params)

# Plot output
plot(output)
