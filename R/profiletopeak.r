#' @title Summarise profile information
#'
#' @description \code{profiletopeak} converts information from the profile list into a matrix, partly resembling
#' a matrix of peaks (often called peaklist).
#'
#' @param profileList A profile list.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#'
#' @return Matrix. See \code{colnames} of return for more information.
#' 
#' @details enviMass workflow function
#' 
#' @export

profiletopeak<-function(profileList,progbar){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
	############################################################################
	# define "fake" peaklist ###################################################
	peaklist<-matrix(ncol=17,nrow=length(profileList[[7]][,8]),0)
	colnames(peaklist)<-c(
		"mean_m/z","mean_intensity","mean_RT","max_intensity", #1
		"in_blind?","above_blind?","var_mz","min_RT","max_RT", #5
		"profileID","number_peaks_total","past_incident","current_incident", #10
		"current_intensity","Component","Homologue","Mass defect" #14
	)
    if(progbar==TRUE){  prog<-winProgressBar("Convert profiles to peaklist",min=1,max=length(profileList[[7]][,8]));
						setWinProgressBar(prog, 0, title = "Convert profiles to peaklist...", label = NULL);}
	for(k in 1:length(profileList[[7]][,8])){
		if(progbar==TRUE){setWinProgressBar(prog, k, title = "Convert profiles to peaklist...", label = NULL)}
		peaklist[k,1]<-profileList[[7]][k,14]
		peaklist[k,2]<-mean(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),2])		
		peaklist[k,3]<-profileList[[7]][k,15]
		peaklist[k,4]<-max(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),2])	
		peaklist[k,5]<-profileList[[7]][k,8] # in blind
		peaklist[k,6]<-profileList[[7]][k,9] # mean sample above mean blind?
		peaklist[k,7]<-var(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),1])	
		peaklist[k,8]<-min(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),3])	
		peaklist[k,9]<-max(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),3])	
		peaklist[k,10]<-profileList[[7]][k,4]	
		peaklist[k,11]<-profileList[[7]][k,3]	
		peaklist[k,12]<-profileList[[7]][k,6] # overall trend
		peaklist[k,13]<-profileList[[7]][k,5] # current trend
		peaklist[k,14]<-profileList[[7]][k,17] # current intensity
	}
	if(progbar==TRUE){ close(prog); }
    ############################################################################
	return(peaklist);
	
}


