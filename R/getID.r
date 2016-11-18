#' @title Find smallest available (file, ...) ID
#'
#' @export
#'
#' @description Given a list of numerical IDs, \code{IDs} returns the smallest integer not part of that list. 
#'
#' @param IDs Vector of integers
#' 
#' @details enviMass workflow function
#' 

getID<-function(IDs){

	################################################################
	if(!is.numeric(IDs)){stop("IDs not numeric")}
	################################################################
	IDs<-IDs[order(IDs,decreasing=FALSE)]
	lowestID<-(0)
	for(i in 1:length(IDs)){
	    if(!any(IDs==i)){
		    lowestID<-i;
			break;
		}
	}
	if(lowestID==0){
	    lowestID<-(max(IDs)+1);
	}
	################################################################
	return(lowestID)
	
}
