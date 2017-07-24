# Set rng seed
set.seed(1)

# Define variables
timestep <- 0.5
vol <- 5e-15
f <- 6.0221415e14*vol
init_conc <- c(40, 0, 0, 0, 0)

# Create calcium input signal:
# increase Ca from 50 to 600 at 100s, hold for 40s, then drop to 50 again
# (from Dupont_camkii.cps)
x <- seq(0, 400, 1)
y <- append(rep(0, 99), rep(600*f, 41))
y <- append(y, rep(50, 261))
input <- data.frame("time" = x, "Ca" = y)

# Simulate model
output <- sim_camkii(input$time, input$Ca/f, timestep, vol, init_conc)
output <- as.data.frame(output)

# Plot output
colnames(output) <- c("time", "calcium", "W_I", "W_B", "W_P", "W_T", "W_A")
plot(output$time, output$calcium, col="blue", xlim=c(100, 160), ylim = c(0, 55), type="l", xlab="time", ylab="concentration")
lines(output$time, output$W_I, col="black", type = "l")
lines(output$time, output$W_B, col="red", type="l")
lines(output$time, output$W_P, col="green", type="l")
lines(output$time, output$W_T, col="cyan", type="l")
lines(output$time, output$W_A, col="orange", type="l")
legend("topright", legend=c("calcium", "W_I", "W_B", "W_P", "W_T", "W_A"),
                   col=c("blue", "black", "red", "green", "cyan", "orange"),
                   lty=c(1,1))
