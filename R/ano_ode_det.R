#' Specific wrapper function for ano calling the LSODA simulation function (provided by package deSolve)
#'
#' @param input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nMol/l).
#' @param input_sim_params A NumericVector: contains values for the simulation end ("endTime") and its timesteps ("timestep").
#' @param input_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
#' @return the result of calling the lsoda simulation algorithm from deSolve 
#' @examples
#' detSim_ano()
#' @export
detSim_ano <- function(input_df, input_sim_params, input_model_params) {



  ############################################################################
  ################################# - Model - ################################
  # USER INPUT 1 for new models: define default model parameters 

  # define default model parameters
  default_model_params <- list(vols      = c(vol = 1e-11),
                               init_conc = c(Cl_ext = 30e6,
                                             C = 1,
                                             C_c = 0,
                                             C_1 = 0,
                                             C_1c = 0,
                                             C_2 = 0,
                                             C_2c = 0,
                                             O = 0,
                                             O_c = 0,
                                             O_1 = 0,
                                             O_1c = 0,
                                             O_2 = 0,
                                             O_2c = 0),
                               params    = c(Vm = -0.06,
                                             T = 293.15,
                                             a1 = 0.0077,
                                             b1 = 917.1288,
                                             k01 = 0.5979439,
                                             k02 = 2.853,
                                             acl1 = 1.8872,
                                             bcl1 = 5955.783,
                                             kccl1 = 1.143e-12,
                                             kccl2 = 0.0009,
                                             kocl1 = 1.1947e-06,
                                             kocl2 = 3.4987,
                                             za1 = 0,
                                             zb1 = 0.0064,
                                             zk01 = 0,
                                             zk02 = 0.1684,
                                             zacl1 = 0.1111,
                                             zbcl1 = 0.3291,
                                             zkccl1 = 0.1986,
                                             zkccl2 = 0.0427,
                                             zkocl1 = 0.6485,
                                             zkocl2 = 0.03,
                                             l = 41.6411,
                                             L = 0.1284,
                                             m = 0.0102,
                                             M = 0.0632,
                                             h = 0.3367,
                                             H = 14.2956))
  ################################# - Model - ################################
  ############################################################################



  # force input dataframe order (First column = time, second column = Ca)
  input_df <- data.frame(time = input_df$time, Ca = input_df$Ca)
  # get default vectors and user supplied input vectors (if they exist)
  # set user vectors to default first so that if user vectors don't exist: they just contain the default values 
  default_vols <- as.list(default_model_params[["vols"]])
  default_init_conc <- as.list(default_model_params[["init_conc"]]) 
  default_params <- as.list(default_model_params[["params"]])
  user_vols <- default_vols
  if(exists("vols", where = input_model_params)) {
    user_vols <- input_model_params[["vols"]]
  } else {
    print("Default volume(s) have been used.")
  }
  user_init_conc <- default_init_conc
  if(exists("init_conc", where = input_model_params)) {
    user_init_conc <- input_model_params[["init_conc"]] 
  } else {
    print("Default initial condition(s) have been used.")
  }
  user_params <- default_params
  if(exists("params", where = input_model_params)) {
    user_params <- input_model_params[["params"]]
  } else {
    print("Default reaction parameter(s) have been used.")
  }



  # UPDATE VALUES 
  # Replace entries in vectors copied from default_model_params with user-supplied values if necessary
  # 1.) Volumes update: 
  for (param_name in names(user_vols)){
    if(exists(param_name, where = default_vols)){
      default_vols[[param_name]] = user_vols[[param_name]]
    } else {
      stop("No such index! Check input parameter vectors.")
    }
  }
  # 2.) Initial conditions update:
  for (param_name in names(user_init_conc)){
    if(exists(param_name, where = default_init_conc)){
      default_init_conc[[param_name]] = user_init_conc[[param_name]]
    } else {
      stop("No such index! Check input parameter vectors.")
    }
  }
  # 3.) Differential equation parameters update:
  for (param_name in names(user_params)){
    if(exists(param_name, where = default_params)){
      default_params[[param_name]] = user_params[[param_name]]
    } else {
      stop("No such index! Check input parameter vectors.")
    }
  }  



  # DEFINE OUTPUT TIMES
  # 1.) simulation output times can be generated from timestep and endTime (evenly spaced)
  # (use default simulation output parameters if none are supplied by user)
  # default simulation timestep value = 0.01
  timestep = 0.01
  if(exists("timestep", where = input_sim_params)){
    timestep = input_sim_params[["timestep"]]
  } else {
    print("Default simulation timestep (0.01s) has been used.")
  }
  # default simulation end time = 100s
  endTime = 100
  if(exists("endTime", where = input_sim_params)) {
    endTime = input_sim_params[["endTime"]]
  } else {
    print("Default simulation end time (100s) has been used.")
  }
  # 2.) simulation output times can be supplied as vector by user (even or unevenly spaced) 
  # default simulation output times vector = seq([first entry of time column in input_df],100,0.01)
  user_output_times_vector = seq(input_df$time[1],100,0.01)
  output_times_set = FALSE
  if(exists("outputTimes", where = input_sim_params)) {
    user_output_times_vector = input_sim_params[["outputTimes"]]
    output_times_set = TRUE
  }
  if(isTRUE(output_times_set)) {
    # cut off input_df$time to get only values that are lesser than the last value of user_output_times_vector
    input_df_subset <- subset(input_df, time <= user_output_times_vector[length(user_output_times_vector)], select = c(time, Ca))
    # take user supplied vector as simulation output time vector
    output_times <- user_output_times_vector
  } else {
    # cut off input_df$time to get only values that are lesser equal endTime (to have output that only goes to endTime)
    input_df_subset <- subset(input_df, time <= endTime, select = c(time, Ca))
    # create simulation output time vector from timestep and endTime
    output_times <- seq(input_df$time[1], endTime, timestep)
  }



  ############################################################################
  ################################# - Model - ################################
  # USER INPUT 2 for new models: define differential equations  

  # model description function
  ano_ode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      # define model ODEs
      
      # faraday constant in C / mol
      faradayConst <- 96485.3329
      # gas constant in C * V / K / mol
      gasConst <- 8.3144598
      # temperature in Kelvin (293.15K = 20°C)
      temp = 293.15
      Vterm <- faradayConst * Vm / (gasConst * temp)

      dC <- -(a1 * exp(za1 * Vterm) * C - b1 * exp(-zb1 * Vterm) * O) + -(k01 * exp(zk01 * Vterm) * 2 * Ca * C - l/L * k02 * exp(-zk02 * Vterm) * C_1) + -(kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C - kccl2 * exp(-zkccl2 * Vterm) * C_c)

      dC_c <- (kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C - kccl2 * exp(-zkccl2 * Vterm) * C_c) + -(acl1 * exp(zacl1 * Vterm) * C_c - bcl1 * exp(-zbcl1 * Vterm) * O_c) + -(h/H * k01 * exp(zk01 * Vterm) * 2 * Ca * C_c - l/L * k02 * exp(-zk02 * Vterm) * C_1c)

      dC_1 <- (k01 * exp(zk01 * Vterm) * 2 * Ca * C - l/L * k02 * exp(-zk02 * Vterm) * C_1) + -(l * a1 * exp(za1 * Vterm) * C_1 - L * b1 * exp(-zb1 * Vterm) * O_1) + -(k01 * exp(zk01 * Vterm) * Ca * C_1 - l/L * 2 * k02 * exp(-zk02 * Vterm) * C_2) + -(h * kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C_1 - H * kccl2 * exp(-zkccl2 * Vterm) * C_1c)

      dC_1c <- (h/H * k01 * exp(zk01 * Vterm) * 2 * Ca * C_c - l/L * k02 * exp(-zk02 * Vterm) * C_1c) + (h * kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C_1 - H * kccl2 * exp(-zkccl2 * Vterm) * C_1c) + -(H*m*l/M * acl1 * exp(zacl1 * Vterm) * C_1c - h*L * bcl1 * exp(-zbcl1 * Vterm) * O_1c) + -(h/H * k01 * exp(zk01 * Vterm) * Ca * C_1c - l/L * 2 * k02 * exp(-zk02 * Vterm) * C_2c)

      dC_2 <- -(l^2 * a1 * exp(za1 * Vterm) * C_2 - L^2 * b1 * exp(-zb1 * Vterm) * O_2) + -(h^2 * kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C_2 - H^2 * kccl2 * exp(-zkccl2 * Vterm) * C_2c) + (k01 * exp(zk01 * Vterm) * Ca * C_1 - l/L * 2 * k02 * exp(-zk02 * Vterm) * C_2)

      dC_2c <- (h^2 * kccl1 * exp(zkccl1 * Vterm) * Cl_ext * C_2 - H^2 * kccl2 * exp(-zkccl2 * Vterm) * C_2c) + -(H*m*l^2/M^2 * acl1 * exp(zacl1 * Vterm) * C_2c - h^2*L^2 * bcl1 * exp(-zbcl1 * Vterm) * O_2c) + (h/H * k01 * exp(zk01 * Vterm) * Ca * C_1c - l/L * 2 * k02 * exp(-zk02 * Vterm) * C_2c)

      dO <- (a1 * exp(za1 * Vterm) * C - b1 * exp(-zb1 * Vterm) * O) + -(k01 * exp(zk01 * Vterm) * 2 * Ca * O - k02 * exp(-zk02 * Vterm) * O_1) + -(kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O - kocl2 * exp(-zkocl2 * Vterm) * O_c)

      dO_c <- (kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O - kocl2 * exp(-zkocl2 * Vterm) * O_c) + -(m/M * k01 * exp(zk01 * Vterm) * 2 * Ca * O_c - k02 * exp(-zk02 * Vterm) * O_1c) + (acl1 * exp(zacl1 * Vterm) * C_c - bcl1 * exp(-zbcl1 * Vterm) * O_c)

      dO_1 <- (k01 * exp(zk01 * Vterm) * 2 * Ca * O - k02 * exp(-zk02 * Vterm) * O_1) + -(k01 * exp(zk01 * Vterm) * Ca * O_1 - 2 * k02 * exp(-zk02 * Vterm) * O_2) + -(m * kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O_1 - M * kocl2 * exp(-zkocl2 * Vterm) * O_1c) + (l * a1 * exp(za1 * Vterm) * C_1 - L * b1 * exp(-zb1 * Vterm) * O_1)

      dO_1c <- (m/M * k01 * exp(zk01 * Vterm) * 2 * Ca * O_c - k02 * exp(-zk02 * Vterm) * O_1c) + (m * kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O_1 - M * kocl2 * exp(-zkocl2 * Vterm) * O_1c) + -(m/M * k01 * exp(zk01 * Vterm) * Ca * O_1c - 2 * k02 * exp(-zk02 * Vterm) * O_2c) + (H*m*l/M * acl1 * exp(zacl1 * Vterm) * C_1c - h*L * bcl1 * exp(-zbcl1 * Vterm) * O_1c)

      dO_2 <- (l^2 * a1 * exp(za1 * Vterm) * C_2 - L^2 * b1 * exp(-zb1 * Vterm) * O_2) + (k01 * exp(zk01 * Vterm) * Ca * O_1 - 2 * k02 * exp(-zk02 * Vterm) * O_2) + -(m^2 * kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O_2 - M^2 * kocl2 * exp(-zkocl2 * Vterm) * O_2c)

      dO_2c <- (H*m*l^2/M^2 * acl1 * exp(zacl1 * Vterm) * C_2c - h^2*L^2 * bcl1 * exp(-zbcl1 * Vterm) * O_2c) + (m/M * k01 * exp(zk01 * Vterm) * Ca * O_1c - 2 * k02 * exp(-zk02 * Vterm) * O_2c) + (m^2 * kocl1 * exp(zkocl1 * Vterm) * Cl_ext * O_2 - M^2 * kocl2 * exp(-zkocl2 * Vterm) * O_2c)

      #dAF_integral <- (O + O_c + O_1 + O_1c + O_2 + O_2c)

      # return list (=state vector) with differentials (and 0 for Ca and Cl_ext since they are external signals) 
      list(c(0, 0, dC, dC_c, dC_1, dC_1c, dC_2, dC_2c, dO, dO_c, dO_1, dO_1c, dO_2, dO_2c))
    })
  }
  ################################# - Model - ################################
  ############################################################################



  # create simulation times:
  # event times get priority if they are very close to output times
  lsodaTimes = sort(c(
      deSolve::cleanEventTimes(output_times, input_df_subset$time),
      input_df_subset$time
  ))
  # modifies the state vector y according to the event data (current calcium value) for given simulation timepoint t
  # all event data columns: 1 = time, 2 = Ca
  eventFun <- function(t, y, ...) {
    # get Ca value at sim timepoint t
    curr_Ca_val <- input_df_subset[input_df_subset[[1]] == t, 2]
    # set first species in state vector y (which is Ca) to current Ca value
    y[1] <- curr_Ca_val
    # return the new state vector
    y
  }
    
  
  
  ############################################################################
  ############################# - Simulation - ###############################
  # USER INPUT 3 for new models: adapt name of 'func' argument 
  #                              (ano_ode for ano model, etc.)
  
  # simulate model with LSODA
  output <- deSolve::lsoda(y = c(Ca = input_df$Ca[[1]], unlist(default_init_conc)),
                           times = lsodaTimes,
                           func = ano_ode,
                           parms = unlist(default_params),
                           event = list(func = eventFun, time = input_df_subset$time))

  # cut out all event times that were not in output time vector 'times'
  output <- as.data.frame(output)
  output <- output[output$time %in% deSolve::nearestEvent(output_times, lsodaTimes), ]

  # return output matrix
  output
  ############################# - Simulation - ###############################
  ############################################################################
}