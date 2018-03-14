# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- list(timestep = 0.1,
                   endTime = 200)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 5e-15),
                     init_conc = c(W_I = 40,
                                   W_B = 0,
                                   W_P = 0,
                                   W_T = 0,
                                   W_A = 0))


# Create calcium input signal (unit: concentration nmol/l):
# increase Ca from 50 to 600 at 100s, hold for 40s, then drop to 50 again
# (from Dupont_camkii.cps)
x <- seq(0, 400, 1)
y <- append(rep(0, 99), rep(600, 41))
y <- append(y, rep(50, 261))
input_df <- data.frame("time" = x, "Ca" = y)

start.time <- as.numeric(Sys.time())*1000

# Simulate model
output <- sim_camkii(input_df, sim_params, model_params)

end.time <- as.numeric(Sys.time())*1000

time.taken <- end.time - start.time
cat(time.taken)

# Plot output
par(mar = c(5,5,2,5))
colnames(output) <- c("time", "calcium", "W_I", "W_B", "W_P", "W_T", "W_A")
plot(output$time, output$calcium, col="blue", xlim=c(90, 160), ylim = c(0, 55), type="l", xlab="time [s]", ylab="CamKII [nmol/l]")
lines(output$time, output$W_I, col="black", type = "l")
lines(output$time, output$W_B, col="red", type="l")
lines(output$time, output$W_P, col="green", type="l")
lines(output$time, output$W_T, col="cyan", type="l")
lines(output$time, output$W_A, col="orange", type="l")
axis(side = 4)
mtext(side = 4, line = 3, 'calcium [a.u]')
legend("topright", legend=c("calcium", "W_I", "W_B", "W_P", "W_T", "W_A"),
                   col=c("blue", "black", "red", "green", "cyan", "orange"),
                   lty=c(1,1))
