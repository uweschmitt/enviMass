#' @title Compare two lists for equality of entries
#'
#' @export
#'
#' @description Compare two lists for equality of entries
#'
#' @param list1. First list
#' @param list2. second list
#'
#' @return Logical whether all pairwise entries are equal
#' 
#' @details  enviMass workflow function.
#' 

comp_list<-function(
	list1,
	list2
){

	##############################################################################
	if(length(list1)!=length(list2)){stop("comp_list: arguments not of same length!")}
	##############################################################################
	for(i in 1:length(list1)){
		if(list1[[i]]!=list2[[i]]){
			return(FALSE)
		}
	}
	##############################################################################
	return(TRUE)
  
}


































