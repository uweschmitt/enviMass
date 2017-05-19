#' @title Determine if a profile contains blind peaks
#'
#' @export
#'
#' @description \code{intensup} runs a trend detection and blind subtraction for the list of time profiles.
#'
#' @param profileList A profile list.
#'
#' @return Updated \code{profileList[[7]]}.
#' 
#' @details enviMass workflow function for trend detection and blind subtraction.
#' 

in_blind<-function(
	profileList	
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
	############################################################################
	for(i in 1:length(profileList[[7]][,1])){
		n<-profileList[[7]][i,1]
		m<-profileList[[7]][i,2]
		for(j in n:m){
			if( (profileList[[9]][profileList[[4]]==as.character(profileList[[2]][j,6])])=="blank"){
				profileList[[7]][i,8]<-1
				profileList[[7]][i,11]<-(profileList[[7]][i,11]+1)
			}
		}
		profileList[[7]][i,10]<-((m-n)-profileList[[7]][i,11])
	}
    ############################################################################
    return(profileList)

}
