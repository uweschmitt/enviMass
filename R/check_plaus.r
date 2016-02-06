#' @title Check measured pattern plausibility
#'
#' @description Check measured pattern plausibility
#'
#' @param cent_peak_mat
#' @param profileList
#' @param RT_tol_inside
#' 
#' @details enviMass workflow function
#' 
	check_plaus<-function(cent_peak_mat,profileList,RT_tol_inside){
		# if only one row with (centroid,peak) left = all is plausible
		# - plausibility of this one (centroid,peak) must be given at this stage, can be inherited further!
		# unique centroids combined?
		if(any(duplicated(cent_peak_mat[,1]))){return(FALSE)}
		# unique peaks combined?
		if(any(duplicated(cent_peak_mat[,2]))){return(FALSE)}
		# all within small RT window ?
		rangeRT<-range(profileList[[2]][cent_peak_mat[,2],3])
		if((rangeRT[2]-rangeRT[1])>RT_tol_inside){return(FALSE)}
		# does intensity pattern match?
		if(length(cent_peak_mat[,1])>1){ # ... which needs more than one peak
	# FINISH			
	# use intercept = 0?
		}
		# else all is plausible
		return(TRUE)
	}
