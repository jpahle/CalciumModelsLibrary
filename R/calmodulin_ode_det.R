# Wrapper for general ode simulation function (provided by package deSolve)
calmodulin_detSim <- function(input_df, sim_params, input_model_params) {
  
  # create calcium input function from Ca time series that can be evaluated at specific time points
  Ca_input_fun <- approxfun(input_df$time, input_df$Ca)
  # create simulation output time vector
  times <- seq(input_df$time[1], sim_params[["endTime"]], sim_params[["timestep"]])
 
  
  
  
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
      
      # Ca input
      Ca <- Ca_input_fun(t)
      
      # Model ODEs
      dCal_off <- k_off * Cal_on - ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off
      dCal_on  <- ((k_on * Ca^h)/(Km^h + Ca^h)) * Cal_off - k_off * Cal_on 
      
      res <- c(dCal_off, dCal_on)
      list(res, signal = Ca)
    })
  }

  ################################# - Model - ################################
  ############################################################################
  
  

  
  # Compare input model parameters with defaults
  # if new value for existing parameter is supplied -> overwrite defaults
  for (param_name in names(input_model_params)){
    if(exists(param_name, where = default_params)){
      default_params[param_name] = input_model_params[param_name]
    }
  }
  # Simulate Model
  # USER INPUT for new models: adapt name of 'func' argument (calmodulin_ode for calmodulin model, etc.)
  output <- deSolve::ode(y = input_model_params[["init_conc"]],
                         times = times,
                         func = calmodulin_ode,
                         parms = unlist(default_params))
}








