#' @title Find minimum or maximum date
#'
#' @description Extracts the minimum or maximum date from a vector of date character strings
#'
#' @param Date. Character string. vector with dates in "YYYY-MM-DD" format.
#' @param get_min. Logical. Else extract maximum date.
#' 
#' @details enviMass workflow function
#' 

minDate<-function(Date,get_min){
	
	################################################
	if(!is.character(Date)){stop("Date not a character string.")}
	if(!is.logical(get_min)){stop("get_min not logical.")}
	################################################
	if(length(Date)==1){return(Date)}
	################################################	
	the_Date<-as.POSIXct(Date)
	if(get_min){
		the_Date<-min(new_Date)[1]
	}else{
		the_Date<-max(new_Date)[1]	
	}
	the_Date<-strsplit(as.character(the_Date)," ")[[1]][1]
	################################################	
	return(the_Date);
	
}