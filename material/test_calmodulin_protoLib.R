# Set rng seed
set.seed(1)

# Define variables
timestep <- 0.05
vol <- 5e-15
f <- 6.0221415e14*vol
init_conc <- c(5, 0)

# Read Ca timeseries
# Calcium input is expected to be in nM!
input <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))

# Simulate model
output <- sim_pkc(input$time, input$Ca/f, timestep, vol, init_conc)
output <- as.data.frame(output)

# Plot output
colnames(output) <- c("time", "Prot_inact", "Prot_act", "calcium")
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), type="l", xlab="time", ylab="concentration")
lines(output$time, output$Prot_act, col="red", type="l")
legend("topright", legend=c("calcium", "Prot_act"), col=c("blue", "red"), lty=c(1,1))
