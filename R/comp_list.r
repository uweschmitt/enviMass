#' @title Compare two lists for matching elements in the first list with at least one in the second list
#'
#' @export
#'
#' @description Compare two lists for equality of entries
#'
#' @param list1. First list
#' @param list2. Second list
#' @param as_pairs. Strictly pairwise matching?
#'
#' @return Logical whether all pairwise entries are equal
#' 
#' @details  enviMass workflow function.
#' 

comp_list<-function(
	list1,
	list2,
	as_pairs=TRUE
){

	##############################################################################
	if(length(list1)!=length(list2) & as_pairs){stop("comp_list: arguments not of same length!")}
	##############################################################################
	if(as_pairs){
		for(i in 1:length(list1)){
			if(list1[[i]]!=list2[[i]]){
				return(FALSE)
			}
		}
	}else{
		for(i in 1:length(list1)){
			found<-FALSE
			for(j in 1:length(list2)){
				if(list1[[i]]==list2[[j]]){
					found<-TRUE;
					break; # innermost loop
				}
			}
			if(!found){
				return(FALSE)
			}
		}
	}
	##############################################################################
	return(TRUE)
  
}


































