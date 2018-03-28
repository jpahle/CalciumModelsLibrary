# Set rng seed for debugging tests
#set.seed(1)

# Simulation Parameters (Vector)
sim_params <- list(timestep = 0.05,
                   endTime = 100)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 5e-14),
                     init_conc = c(Prot_inact = 5,
                                   Prot_act = 0))

# Read Ca timeseries
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert part number from input table to concentration (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f

start.time <- as.numeric(Sys.time())*1000

# Simulate model
output <- sim_glycphos(input_df, sim_params, model_params)

end.time <- as.numeric(Sys.time())*1000

time.taken <- end.time - start.time
cat(time.taken)

# Plot output
colnames(output) <- c("time", "Ca", "Prot_inact", "Prot_act")
# define active glycphos fraction
glycphos_active_frac <- with(output, Prot_act/(Prot_inact+Prot_act))
par(mar = c(5,5,2,5))
plot(output$time, output$Ca, col="blue", type="l", xlab="time [s]", ylab="glycogen phos. [nmol/l]")
lines(output$time, output$Prot_act, col="red", type="l")
lines(output$time, glycphos_active_frac, col="green", type="l")
#axis(side = 4)
#mtext(side = 4, line = 3, 'calcium [a.u]')
legend("topright", legend=c("calcium", "glycogen phos.", "active fraction [0,1]"), col=c("blue", "red", "green"), lty=c(1,1,1))
