#' @title Summarise profile information
#'
#' @description \code{profiletopeak2}
#'
#' @param profileList A profile list (profiled, not just loaded or partitioned).
#' @param logfile enviMass project logfile.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#' @param get_what To be completed ...
#' @param at_ID To be completed ...
#'
#' @return Matrix or integer
#' 
#' @details enviMass workflow function
#' 

profiletopeak2<-function(
	profileList,
	logfile,
	progbar,
	get_what="peaks",
	at_ID=FALSE
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
	############################################################################
	# get latest sample & blind ID #############################################
	atPOSIX<-profileList[[3]];
    sampletype<-profileList[[9]];
    sampleID<-profileList[[4]];
    atdate<-c();
    attime<-c();
    for(i in 1:length(atPOSIX)){
          atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
          attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
    }
    attime<-as.difftime(attime);
    atdate<-as.Date(atdate);
    ord<-order(as.numeric(atdate),as.numeric(attime),sampleID);
    atPOSIXsort<-atPOSIX[ord];
    atdate<-atdate[ord];
    attime<-attime[ord];
	sampleID<-sampleID[ord];
	sampletype<-sampletype[ord];
    timeset<-matrix(nrow=length(atPOSIX),ncol=3,0);
    for(i in 1:length(sampleID)){
      if(sampletype[i]=="sample"){
        timeset[i,2]<-as.numeric(sampleID[i]);
      }
      if(sampletype[i]=="blank"){
        timeset[i,3]<-as.numeric(sampleID[i]);
      }
    }
	############################################################################	
	get_ID<-timeset[timeset[,2]!=0,2]
	if(at_ID=="FALSE"){
		ID_sample<-get_ID[length(get_ID)]
		if(get_what=="ID_sample"){
			return(ID_sample)
		}
	}else{
		ID_sample<-at_ID;
	}
	############################################################################
	get_ID<-timeset[timeset[,3]!=0,3]
	ID_blank<-get_ID[length(get_ID)]
	if(get_what=="ID_blank"){
		return(ID_blank)
	}
	############################################################################
	# retrieve data ############################################################
	peak_prof<-matrix(ncol=14,nrow=length(profileList[[7]][,1]),0)
	colnames(peak_prof)<-c(
		"mean_mz","mean_rt","mean_int",
		"sample_mz","sample_rt","int_sample",
		"blind_mz","blind_rt","int_blind",
		"peak_counts","profile_ID","component_ID",
		"old_trend","new_trend")
    if(progbar==TRUE){  prog<-winProgressBar("Convert profiles to peaklist",min=1,max=length(profileList[[7]][,8]));
						setWinProgressBar(prog, 0, title = "Convert profiles to peaklist...", label = NULL);}
	for(k in 1:length(profileList[[7]][,1])){
		if(progbar==TRUE){setWinProgressBar(prog, k, title = "Convert profiles to peaklist...", label = NULL)}
		peak_prof[k,1]<-profileList[[7]][k,14][[1]]
		peak_prof[k,2]<-profileList[[7]][k,15][[1]]
		peak_prof[k,3]<-profileList[[7]][k,16][[1]]
		if(any(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),6]==ID_sample)){
			getit<-(profileList[[7]][k,1]:profileList[[7]][k,2])
			getit<-getit[profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),6]==ID_sample]
			peak_prof[k,4]<-profileList[[2]][getit,1]
			peak_prof[k,5]<-profileList[[2]][getit,2]
			peak_prof[k,6]<-profileList[[2]][getit,3]
		}
		if(any(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),6]==ID_blank)){
			getit<-(profileList[[7]][k,1]:profileList[[7]][k,2])
			getit<-getit[profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),6]==ID_blank]
			peak_prof[k,7]<-profileList[[2]][getit,1]
			peak_prof[k,8]<-profileList[[2]][getit,2]
			peak_prof[k,9]<-profileList[[2]][getit,3]
		}
		peak_prof[k,10]<-profileList[[7]][k,3][[1]]
		peak_prof[k,11]<-profileList[[7]][k,4]
	}
	if(progbar==TRUE){ close(prog); }
	############################################################################
	return(peak_prof);
}
