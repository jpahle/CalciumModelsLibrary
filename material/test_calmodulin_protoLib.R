# Set rng seed for debugging tests
set.seed(1)

# Simulation Parameters (Vector)
sim_params <- list(outputTimes = c(seq(0,20,0.1),30,40,50,60,120))
#sim_params <- list()
#sim_params <- list(endTime = 50, timestep = 2.6)
#sim_params <- list(outputTimes = c(0,1))

# Model Parameters (List)
model_params <- list(vols      = c(vol = 5e-14),
                     init_conc = c(Prot_inact = 5,
                                   Prot_act =   0))

# Read Ca timeseries
input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert part number from input table to concentration (c*f=n since f = Avogadro*Vol)
f <- 6.0221415e14*model_params[["vols"]][["vol"]]
input_df["Ca"] <- input_df["Ca"]/f

# Simulate model
output <- sim_calmodulin(input_df, sim_params, model_params)

# Plot output
par(mar = c(5,5,2,5))
colnames(output) <- c("time", "calcium", "Prot_inact", "Prot_act")
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), ylim=c(0,15), type="l", xlab="time [s]", ylab="calmodulin [nmol/l]")
lines(output$time, output$Prot_act, col="red", type="l")
axis(side = 4)
mtext(side = 4, line = 3, 'calcium [a.u]')
legend("topright", legend=c("calcium", "calmodulin"), col=c("blue", "red"), lty=c(1,1))
