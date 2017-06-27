# Define variables (values from copasi file dupont_camkii.cps)
timestep <- 0.5
vol <- 5e-15
f <- 6.0221415e14*vol
a <- -0.22
b <- 1.826
c <- 0.1
k_IB <- 0.01
k_BI <- 0.8
k_PT <- 1
k_TP <- 1e-12
k_TA <- 0.0008
k_AT <- 0.01
k_AA <- 0.29
c_I <- 0
c_B <- 0.75
c_P <- 1
c_T <- 0.8
c_A <- 0.8
camT <- 1000
Kd <- 1000
Vm_phos <- 0.005
Kd_phos <- 0.3
totalC <- 40
Wi_conc <- 40
Wb_conc <- 0
Wp_conc <- 0
Wt_conc <- 0
Wa_conc <- 0
h <- 4;

# Create calcium input signal:
# increase Ca from 50 to 600 at 100s, hold for 40s, then drop to 50 again
# (from Dupont_camkii.cps)
x <- seq(0, 400, 1)
y <- append(rep(0, 99), rep(600*f, 41))
y <- append(y, rep(50, 261))
input <- data.frame("time" = x, "Ca" = y)

# Simulate model
output <- camkii(input$time, input$Ca/f, timestep, vol, a, b, c, k_IB, k_BI, k_PT, k_TP, k_TA, k_AT, k_AA, c_I, c_B, c_P, c_T, c_A, camT, Kd, Vm_phos, Kd_phos, totalC, Wi_conc, Wb_conc, Wp_conc, Wt_conc, Wa_conc, h)

# Plot output
output <- as.data.frame(output)
colnames(output) <- c("time", "W_I", "W_B", "W_P", "W_T", "W_A", "calcium")
plot(output$time, output$calcium, col="blue", xlim=c(100, 160), ylim = c(0, 55), type="l", xlab="time", ylab="concentration")
lines(output$time, output$W_I, col="black", type = "l")
lines(output$time, output$W_B, col="red", type="l")
lines(output$time, output$W_P, col="green", type="l")
lines(output$time, output$W_T, col="cyan", type="l")
lines(output$time, output$W_A, col="orange", type="l")
legend("topright", legend=c("calcium", "W_I", "W_B", "W_P", "W_T", "W_A"),
                   col=c("blue", "black", "red", "green", "cyan", "orange"),
                   lty=c(1,1))
