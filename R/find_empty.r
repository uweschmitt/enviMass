#' @title Find entry gaps in a list obejct
#'
#' @export
#'
#' @description Find entry gaps in a list obejct
#'
#' @param any_list A list object
#'
#' @return Vector of numerics indexing list gaps.
#' 
#' @details  enviMass workflow function.
#' 

find_empty<-function(
	any_list
){

	##############################################################################
	if(!is.list(any_list)){stop("any_list mut be a list")}
	if(length(any_list)==0){return(c())}
	##############################################################################
	positions<-c()
	for(i in 1:length(any_list)){
		if(length(any_list[[i]])==0){
			positions<-c(positions,i)
		}
	}	
	##############################################################################
	return(positions)
  
}


































