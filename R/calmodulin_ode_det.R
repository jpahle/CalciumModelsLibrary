# Wrapper for general ode simulation function (provided by package deSolve)
calmodulin_detSim <- function(input_df, sim_params, input_model_params) {

  # Force input dataframe order
  input_df <- data.frame(time = input_df$time, Ca = input_df$Ca)
  # TODO: Check that all Ca time points in input_df are unique
  # ...
  # create simulation output time vector
  output_times <- seq(input_df$time[1], sim_params[["endTime"]], sim_params[["timestep"]])



  ############################################################################
  ################################# - Model - ################################
  # USER INPUT for new models: define model parameters
  #                            and differential equations

  # Default model parameters
  default_params <- list(k_on = 0.025,
                         k_off = 0.005,
                         Km = 1.0,
                         h = 4.0)

  # Model description
  calmodulin_ode <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {

      # Model ODEs
      dCal_off <- k_off * Cal_on - ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off
      dCal_on  <- ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off - k_off * Cal_on

      list(c(0, dCal_off, dCal_on))
    })
  }

  ################################# - Model - ################################
  ############################################################################



  # create simulation times
  # event times get priority if they are very close to output times
  lsodaTimes = sort(c(
      deSolve::cleanEventTimes(output_times, input_df$time),
      input_df$time
  ))
  # Modifies the state vector y according to the event data for given time t
  # All event data columns: 1 = time, 2 = Ca
  eventFun <- function(t, y, ...) {
    ca <- input_df[input_df[[1]] == t, 2]

    y[1] <- ca

    y
  }
  # Compare input model parameters with defaults
  # if new value for existing parameter is supplied -> overwrite defaults
  for (param_name in names(input_model_params)){
    if(exists(param_name, where = default_params)){
      default_params[param_name] = input_model_params[param_name]
    }
  }
  # Simulate Model
  # USER INPUT for new models: adapt name of 'func' argument (calmodulin_ode for calmodulin model, etc.)
  output <- deSolve::lsoda(y = c(Ca = input_df$Ca[[1]], input_model_params[["init_conc"]]),
                           times = lsodaTimes,
                           func = calmodulin_ode,
                           parms = unlist(default_params),
                           event = list(func = eventFun, time = input_df$time))

  # Cut out all event times that were not in output time vector 'times'
  # output <- output[output$time %in% deSolve::nearestEvent(times, lsodaTimes), ]

  output
}








