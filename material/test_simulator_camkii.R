# Create calcium input signal:
# increase Ca from 50 to 600 at 100s, hold for 40s, then drop to 50 again
# (like in Dupont_camkii.cps)
x <- seq(0, 400, 1)
y <- append(rep(0, 99), rep(600*f, 41))
y <- append(y, rep(50, 261))
input <- data.frame("time" = x, "Ca" = y)

# Define initial concentrations of model species (Wi, Wb, Wp, Wt, Wa)
init_conc <- c(40, 0, 0, 0, 0)
# Define simulation parameters
timestep <- 0.05
vol <- 5e-14
f <- 6.0221415e14*vol #for conc to particle number conversion
# Simulate Model
output <- simulator(input$time, input$Ca/f, init_conc, camkii_props, camkii_stM, timestep, vol)

# Plot output
output <- as.data.frame(output)
colnames(output) <- c("time", "calcium", "W_I", "W_B", "W_P", "W_T", "W_A")
plot(output$time, output$calcium, col="blue", xlim=c(0, 400), ylim = c(0, 120), type="l", xlab="time", ylab="concentration")
lines(output$time, output$W_I, col="black", type = "l")
lines(output$time, output$W_B, col="red", type="l")
lines(output$time, output$W_P, col="green", type="l")
lines(output$time, output$W_T, col="cyan", type="l")
lines(output$time, output$W_A, col="orange", type="l")
legend("topright", legend=c("calcium", "W_I", "W_B", "W_P", "W_T", "W_A"),
       col=c("blue", "black", "red", "green", "cyan", "orange"),
       lty=c(1,1))
