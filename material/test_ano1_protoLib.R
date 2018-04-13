# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- list(timestep = 0.01,
                   endTime = 20)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 1e-11),
                     init_conc = c(Cl_ext = 3e06))

# Read Ca timeseries
#input_df <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))
# convert part number from input table to concentration (c*f=n since f = Avogadro*Vol)
#f <- 6.0221415e14*model_params[["vols"]][["vol"]]
#input_df["Ca"] <- input_df["Ca"]/f

# Sine(baseline, amp, period, phase, duration, resolution)
input_df <- as.data.frame(OscillatorGenerator::Sine(0, 1e06, 5, 5, 20, 0.01))

colnames(input_df) <- c("time", "Ca")

# Simulate model
start.time <- as.numeric(Sys.time())

output <- sim_ano(input_df, sim_params, model_params)

end.time <- as.numeric(Sys.time())

time.taken <- end.time - start.time
cat(time.taken)

# Plot output
par(mar = c(5,5,2,5))
colnames(output) <- c("time", "calcium", "Cl_ext", "C", "C_c", "C_1", "C_1c", "C_2", "C_2c", "O", "O_c", "O_1", "O_1c", "O_2", "O_2c")
# sum of active ano1 species
active_ano_sum <- with(output, O + O_c + O_1 + O_1c + O_2 + O_2c)
# sum of inactive ano1 species
inactive_ano_sum <- with(output, C + C_c + C_1 + C_1c + C_2 + C_2c)
# active ano fraction
active_ano_frac <- 1/(1+(inactive_ano_sum/active_ano_sum))
plot(output$time, active_ano_frac, col="red", type="l", ylim = c(0,0.06), xlab="time [s]", ylab="active ANO1 fraction [nmol/l]")
par(new=T)
plot(output$time, output$calcium, col="blue", type="l", axes=F, xlab=NA, ylab=NA, cex=1.2)
axis(side = 4)
mtext(side = 4, line = 3, "calcium [a.u.]")
legend("topright", legend=c("calcium", "active ANO1 fraction"),
       col=c("blue", "red"),
       lty=c(1,1))
