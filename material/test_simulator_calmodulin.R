# Read input data
input <- read.table("material/ca5e-14_2.85_1000_0.05s.out", col.names = c("time", "steps", "G_alpha", "PLC", "Ca"))

# Define initial concentrations of model species (E_inact and E_act)
init_conc <- c(5,0)
# Define simulation parameters
timestep <- 0.05
vol <- 5e-14
f <- 6.0221415e14*vol #for conc to particle number conversion
# Simulate Model
output <- simulator(input$time, input$Ca/f, init_conc, calmodulin_props, calmodulin_stM, timestep, vol)

# Plot output
output <- as.data.frame(output)
colnames(output) <- c("time", "calcium", "Prot_inact", "Prot_act")
plot(output$time, output$calcium, col="blue", xlim=c(0, 100), type="l", xlab="time", ylab="concentration")
lines(output$time, output$Prot_act, col="red", type="l")
legend("topright", legend=c("calcium", "Prot_act"), col=c("blue", "red"), lty=c(1,1))
