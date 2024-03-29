# Set rng seed for debugging tests
set.seed(1)

# Simulation Parameters (Vector)
sim_params <- list(timestep = 0.05,
                   endTime = 100)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 5e-14),
                     init_conc = c(Prot_inact = 5,
                                   Prot_act = 0),
                     params    = c(k_on = 0.02,
                                   k_off = 0.1))

# Read Ca timeseries
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert part number from input table to concentration (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f

starttime <- as.numeric(Sys.time())*1000

# Simulate model
output <- sim_calcineurin(input_df, sim_params, model_params)

endtime <- as.numeric(Sys.time())*1000

timetaken <- endtime - starttime
cat(timetaken)

# Plot output
colnames(output) <- c("time", "calcium", "Prot_inact", "Prot_act")
par(mar = c(5,5,2,5))
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), ylim=c(0,15), type="l", xlab="time [s]", ylab="calcineurin [nmol/l]")
lines(output$time, output$Prot_act, col="red", type="l")
axis(side = 4)
mtext(side = 4, line = 3, 'calcium [a.u]')
legend("topright", legend=c("calcium", "calcineurin"), col=c("blue", "red"), lty=c(1,1))
