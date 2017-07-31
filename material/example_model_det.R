########################## - R Wrapper Function - ##########################
# Don't change this generalized block!
# Wrapper for general ode simulation function (provided by package deSolve)
detSim <- function(model, sim_params, input_model_params) {

  # Compare input model parameters with defaults
  # if new value for existing parameter is supplied -> overwrite defaults
  for (param_name in names(input_model_params)){
    if(exists(param_name, where = default_params)){
      default_params[param_name] = input_model_params[param_name]
    }
  }
  # Create time vector
  times = seq(0, sim_params[["endTime"]], by = sim_params[["timestep"]])
  # Simulate Model
  output <- deSolve::ode(y = input_model_params[["init_conc"]],
                times = times,
                func = model,
                parms = unlist(default_params))
}



################################# - Model - ################################
# USER INPUT for new models: define model parameters
#                            and differential equations
# Example Model:
# A -> B with rate constant k1 = 0.1
# B -> A with rate constant k2 = 0.01
# model not dependent on Calcium input

# Default model parameters
default_params <- list(k1 = 0.1, k2 = 0.01)

# Model description
exampleModel <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dA <- -k1 * A + k2 * B
    dB <- k1 * A - k2 * B
    list(c(dA, dB))
  })
}





