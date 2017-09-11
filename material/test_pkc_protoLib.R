# Set rng seed for debugging tests
set.seed(1)

# Simulation parameters (Vector)
sim_params <- c(timestep = 1,
                endTime = 1000)
# Model Parameters (List)
model_params <- list(vols      = c(vol = 1e-15),
                     init_conc = c(PKC_inact = 1000,
                                   CaPKC = 0,
                                   DAGCaPKC = 0,
                                   AADAGPKC_inact = 0,
                                   AADAGPKC_act = 0,
                                   PKCbasal = 20,
                                   AAPKC = 0,
                                   CaPKCmemb = 0,
                                   AACaPKC = 0,
                                   DAGPKCmemb = 0,
                                   DAGPKC = 0))

# Read Ca timeseries
input_df <- read.table("material/Sine_Input.txt", col.names = c("time", "Ca"))

# Simulate model
output <- sim_pkc(input_df, sim_params, model_params)
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
                      "DAGPKC"
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
plot(output$time, output$active_fraction, col="orange", type="l", xlab="time", ylab="concentration", main = "active PKC fraction of total PKC")
par(mfrow = c(1,1))
