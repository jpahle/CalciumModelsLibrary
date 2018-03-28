#' Specific wrapper function for calcineurin calling the LSODA simulation function (provided by package deSolve)
#'
#' @param input_df A Dataframe: the input Calcium time series (with at least two columns: "time" in s and "Ca" in nMol/l).
#' @param input_sim_params A NumericVector: contains values for the simulation end ("endTime") and its timesteps ("timestep").
#' @param input_model_params A List: the model specific parameters. Can contain up to three different vectors named "vols" (model volumes), "init_conc" (initial conditions) and "params" (propensity equation parameters). 
#' @return the result of calling the lsoda simulation algorithm from deSolve 
#' @examples
#' detSim_calcineurin()
#' @export
detSim_calcineurin <- function(input_df, input_sim_params, input_model_params) {



  ############################################################################
  ################################# - Model - ################################
  # USER INPUT 1 for new models: define default model parameters 

  # define default model parameters
  default_model_params <- list(vols      = c(vol = 5e-14),
                               init_conc = c(Prot_inact = 5.0,
                                             Prot_act = 0),
                               params    = c(k_on = 1.0,
                                             k_off = 1.0,
                                             p = 3.0))
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
  calcineurin_ode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      # define model ODEs
      dProt_inact <- k_off * Prot_act - k_on* Ca^p * Prot_inact
      dProt_act  <- k_on* Ca^p * Prot_inact - k_off * Prot_act
      # return list (=state vector) with differentials (and 0 for Ca since it is an external signal) 
      list(c(0, dProt_inact, dProt_act))
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
  #                              (calcineurin_ode for calcineurin model, etc.)
  
  # simulate model with LSODA
  output <- deSolve::lsoda(y = c(Ca = input_df$Ca[[1]], unlist(default_init_conc)),
                           times = lsodaTimes,
                           func = calcineurin_ode,
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