#' @title Mute a logfile parameter in workflow evaluation
#'
#' @description Mutes a logfile parameter in the proccess of workflow evaluation; i.e. a logfile$parameter$x is not validated to be part of a workflow script.
#'
#' @param parameter A logfile parameter
#' 
#' @details enviMass workflow function. The parameter must directly be embraced by the function brackets, i.e., mute(logfile$parameters$x) or mute(logfile$parameters$external$x)
#' 

mute<-function(parameter){
	
	eval(parameter)
	
}

