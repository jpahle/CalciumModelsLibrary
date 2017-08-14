# Specific wrapper function for calmodulin calling the LSODA simulation function (provided by package deSolve)
calmodulin_detSim <- function(input_df, sim_params, input_model_params) {

  # force input dataframe order (First column = time, second column = Ca)
  input_df <- data.frame(time = input_df$time, Ca = input_df$Ca)
  # cut off input_df$time to get only values that are lesser equal endTime (to have output that only goes to endTime)
  input_df_subset <- subset(input_df, time <= sim_params[["endTime"]], select = c(time, Ca))
  # create simulation output time vector
  output_times <- seq(input_df$time[1], sim_params[["endTime"]], sim_params[["timestep"]])
  
  

  ############################################################################
  ################################# - Model - ################################
  # USER INPUT for new models: define model parameters
  #                            and differential equations

  # define default model parameters
  default_params <- list(k_on = 0.025,
                         k_off = 0.005,
                         Km = 1.0,
                         h = 4.0)

  # model description function
  calmodulin_ode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      # define model ODEs
      dCal_off <- k_off * Cal_on - ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off
      dCal_on  <- ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off - k_off * Cal_on
      # return list (=state vector) with differentials (and 0 for Ca since it is an external signal) 
      list(c(0, dCal_off, dCal_on))
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
  # compare user supplied input model parameters with defaults:
  # if new value for existing parameter is supplied -> overwrite defaults
  for (param_name in names(input_model_params)){
    if(exists(param_name, where = default_params)){
      default_params[param_name] = input_model_params[param_name]
    }
  }
  
  
  
  ############################################################################
  ############################# - Simulation - ###############################
  # USER INPUT for new models: adapt name of 'func' argument 
  #                            (calmodulin_ode for calmodulin model, etc.)
  
# simulate model with LSODA
  output <- deSolve::lsoda(y = c(Ca = input_df$Ca[[1]], input_model_params[["init_conc"]]),
                           times = lsodaTimes,
                           func = calmodulin_ode,
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