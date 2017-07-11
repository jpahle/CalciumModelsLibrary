# Define variables
timestep <- 1
vol <- 1e-15
f <- 6.0221415e14*vol
# initial concentrations of model species
# (PKCinact0_conc,.,.,.,.,PKCbasal0_conc,.,.,.,.,.)
init_conc <- c(1000, 0, 0, 0, 0, 20, 0, 0, 0, 0, 0)

# Read Ca timeseries
# Calcium input is expected to be in nM!
input <- read.table("material/Sine_Input.txt", col.names = c("time", "Ca"))

# Simulate model
output <- pkc(input$time, input$Ca/f, timestep, vol, init_conc, calc_propensities, update_system_state)
output <- as.data.frame(output)
colnames(output) <- c("time",
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
                      "DAGPKC",
                      "calcium"
                      )

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
par(mfrow = c(3,1))
plot(output$time, output$calcium, col="blue",type="l", xlab="time", ylab="concentration", main = "Calcium")
plot(output$time, output$PKC_act, col="red", type = "l", xlab="time", ylab="concentration", main = "sum of active PKC species")
plot(output$time, output$active_fraction, col="green", type="l", xlab="time", ylab="concentration", main = "active PKC fraction of total PKC")
par(mfrow = c(1,1))
