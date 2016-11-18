#' @title Reset selection ins the enviMass UI
#'
#' @export
#'
#' @description Reset selection ins the enviMass UI
#'
#' @param session shiny session ID
#' 
#' @details enviMass workflow UI function. 
#' 

reset_selections<-function(session){

	################################################################
	updateSelectInput(session,"Ion_mode_Cal",selected = "none")
	updateNumericInput(session,"sel_meas_ID", value = 0)
	updateNumericInput(session,"sel_meas", value = 0)
	updateSelectInput(session,"Pos_compound_select",selected = "Choose")	
	################################################################
	
}
