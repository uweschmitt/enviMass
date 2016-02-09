#' @title Check measured pattern plausibility
#'
#' @description Check measured pattern plausibility
#'
#' @param cent_peak_combi
#' @param pattern_centro
#' @param profileList
#' @param RT_tol_inside
#' 
#' @details enviMass workflow function
#' 

	check_plaus<-function(cent_peak_combi,pattern_centro,profileList,RT_tol_inside,int_tol){
		# if only one row with (centroid,peak) left = all is plausible
		# - plausibility of this one (centroid,peak) must be given at this stage, can be inherited further!
		# unique centroids combined?
		if(any(duplicated(cent_peak_combi[,1]))){return(FALSE)}
		# unique peaks combined?
		if(any(duplicated(cent_peak_combi[,2]))){return(FALSE)}
		# all within small RT window ?
		rangeRT<-range(profileList[[2]][cent_peak_combi[,2],3])
		if((rangeRT[2]-rangeRT[1])>RT_tol_inside){return(FALSE)}
		# does intensity pattern match?
		if(length(cent_peak_combi[,1])>1){ # ... which needs more than one peak
			for(n in 1:(length(cent_peak_combi[,1])-1)){
				for(m in (n+1):length(cent_peak_combi[,1])){
					ratio_int<-(profileList[[2]][cent_peak_combi[n,2],2]/profileList[[2]][cent_peak_combi[m,2],2])
					ratio_int_theo_high<-(
						(pattern_centro[cent_peak_combi[n,1],1]+(pattern_centro[cent_peak_combi[n,1],1]*int_tol/100))/
						(pattern_centro[cent_peak_combi[m,1],1]-(pattern_centro[cent_peak_combi[m,1],1]*int_tol/100))
					)
					if(ratio_int_theo_high<ratio_int){return(FALSE)}
					ratio_int_theo_low<-(
						(pattern_centro[cent_peak_combi[n,1],1]-(pattern_centro[cent_peak_combi[n,1],1]*int_tol/100))/
						(pattern_centro[cent_peak_combi[m,1],1]+(pattern_centro[cent_peak_combi[m,1],1]*int_tol/100))
					)
					if(ratio_int_theo_low>ratio_int){return(FALSE)}	
				}
			}		
		}
		# else all is plausible
		return(TRUE)
	}
