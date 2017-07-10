# Define variables
timestep <- 1
vol <- 1e-15
f <- 6.0221415e14*vol
k1 <- 1
k2 <- 50
k3 <- 1.2e-7
k4 <- 0.1
k5 <- 1.2705
k6 <- 3.5026
k7 <- 1.2e-7
k8 <- 0.1
k9 <- 1
k10 <- 0.1
k11 <- 2
k12 <- 0.2
k13 <- 0.0006
k14 <- 0.5
k15 <- 7.998e-6
k16 <- 8.6348
k17 <- 6e-7
k18 <- 0.1
k19 <- 1.8e-5
k20 <- 2
AA <- 11000
DAG <- 5000
PKCinact0_conc <- 1000
PKCbasal0_conc <- 20

# Read Ca timeseries
# Calcium input is expected to be in nM!
input <- read.table("material/Sine_Input.txt", col.names = c("time", "Ca"))

# Simulate model
output <- pkc(input$time, input$Ca/f, timestep, vol, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, AA, DAG, PKCinact0_conc, PKCbasal0_conc)
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
