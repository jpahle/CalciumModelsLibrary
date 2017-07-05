# Read Ca timeseries
# Calcium input is expected to be in nM!
input <- read.table("material/Sine_Input_short.txt", col.names = c("time", "Ca"))

# Define initial concentrations of model species (PKCinact0_conc,.,.,.,.,PKCbasal0_conc,.,.,.,.,.)
init_conc <- c(1000, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0)
# Define simulation parameters
timestep <- 1
vol <- 1e-15
f <- 6.0221415e14*vol #for conc to particle number conversion

# Simulate Model
output <- simulator(input$time, input$Ca/f, init_conc, pkc_props, pkc_stM, timestep, vol)
output <- as.data.frame(output)
colnames(output) <- c("time",
                      "calcium",
                      "PKC_inact",
                      "CaPKC",
                      "DAGCaPKC",
                      "AADAGPKC_inact",
                      "AADAGPKC_act",
                      "PKCbasal",
                      "AAPKC",
                      "CaPKCmemb",
                      "AACaPKC",
                      "DAGPKCmemb",
                      "DAGPKC")

# Define relevant output for plotting
output$PKC_act <- with(output, PKCbasal +
                               CaPKCmemb +
                               AACaPKC +
                               DAGPKCmemb +
                               AADAGPKC_act +
                               AAPKC)
output$active_fraction <- with(output, (PKCbasal +
                                        CaPKCmemb +
                                        AACaPKC +
                                        DAGPKCmemb +
                                        AADAGPKC_act +
                                        AAPKC) /
                                       (PKCbasal +
                                        CaPKCmemb +
                                        AACaPKC +
                                        DAGPKCmemb +
                                        AADAGPKC_act +
                                        AAPKC +
                                        PKC_inact))
# Plot output
plot(output$time, output$calcium, xlim = c(0, 1), col="blue",type="l", xlab="time", ylab="concentration")
lines(output$time, output$PKC_act, col="red", type = "l")
lines(output$time, output$active_fraction, col="green", type="l")
legend("topright", legend=c("calcium", "PKC_act", "active_fraction"),
       col=c("blue", "red", "green"),
       lty=c(1,1))
