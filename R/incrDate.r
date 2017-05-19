#' @title Increments date
#'
#' @description Increments a date in "YYYY-MM-DD" formats
#'
#' @param Date. Character string. Date in "YYYY-MM-DD" format.
#' @param increment. Numeric. Number of days to increment with.
#' 
#' @details enviMass workflow function
#' 

incrDate<-function(Date,increment){
	
	################################################
	if(!is.character(Date)){stop("Date not a character string.")}
	if(!is.numeric(increment)){stop("increment not numeric")}
	################################################
	new_Date<-as.POSIXct(Date)
	new_Date<-(new_Date+(increment*60*60*26)) # Zeitumstellung in CET-format!!!
	new_Date<-strsplit(as.character(new_Date)," ")[[1]][1]
	################################################	
	return(new_Date);
	
}
