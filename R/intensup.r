#' @title Trend detection and blind subtraction for time profiles.
#'
#' @export
#'
#' @description \code{intensup} runs a trend detection and blind subtraction for the list of time profiles.
#'
#' @param profileList A profile list.
#' @param from Logical or integer of index.
#' @param to Logical or integer of index.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#' @param blindsub Logical. Run blind subtraction?
#' @param blindfold Numerical. Blind definition; above blind, if \code{blindfold} larger in intensity.
#' @param lags Vector of numericals.
#' @param threshold Numerical. A trend is reported if its intensity is \code{threshold} above the mean intensity plus the intensity deviation of other trends.
#' @param notrend Logical. Report global trend intensity as maximum intensity after blind subtraction.
#'
#' @return Updated \code{profileList[[7]]}.
#' 
#' @details enviMass workflow function for trend detection and blind subtraction.
#' 

intensup<-function(
	profileList,
	from=FALSE,
	to=FALSE,
	progbar=FALSE,
	blindsub=TRUE, # do a blind subtraction?
	blindfold=100, # how much higher in intensity than blind?
	lags=c(5,14),  # time lags
	threshold=3,   # trend threshold: 
	notrend=FALSE  # no global threshold, but global maximum above blind	
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
    if(!from){m=1}else{m=from}
    if(!to){n=length(profileList[[7]][,1])}else{n=to}
	if(blindsub!=FALSE){if(!is.numeric(blindfold) || (blindfold<0)){stop("Invalid blindfold argument; aborted.")}}
    if(blindsub!=FALSE){subit=1;subrat=blindfold;}else{subit=2;subrat=0;}
	if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
	if(!is.logical(notrend)){stop("notrend must be logical.")}
	############################################################################
    # set matrix to sort & store data from a profile ###########################
    atPOSIX<-profileList[[3]];
    sampletype<-profileList[[9]];
    sampleID<-profileList[[4]];
	# filter out other file types such as spiked ones
	keep<-((sampletype=="sample")|(sampletype=="blank"))
	atPOSIX<-atPOSIX[keep]
	sampletype<-sampletype[keep]
	sampleID<-sampleID[keep]
	#
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
    timeset<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),0);
    for(i in 1:length(sampleID)){
      if(sampletype[i]=="sample"){
        timeset[i,2]<-as.numeric(sampleID[i]);
      }
      if(sampletype[i]=="blank"){
        timeset[i,3]<-as.numeric(sampleID[i]);
      }
    }
    numtime<-(as.numeric(atdate)+as.numeric(attime/(24*60*60)))
	colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))			
	leng<-max(seq(1,length(timeset[,1]),1)[timeset[,2]!=0])	
	latestID<-timeset[leng,2][[1]]
	############################################################################
	# check & adjust lags ######################################################
	if(any(lags>(max(numtime)-min(numtime)+1))){
		lags<-lags[lags<=(max(numtime)-min(numtime)+1)]
		cat("WARNING: at least one lag longer than covered time period - omitted!\n")
		if(length(lags)==0){
			stop("...no lags left; aborted.")
		}		
	}
	# INSERT ... also check any of the inter-sample distances ##################
	if(max(diff(numtime[timeset[,2]!=0],lag=1))>max(lags)){
		warning("\n At least one gap in sample time series larger than largest lag!")
	}
	############################################################################
    profileList[[7]][,5:7]<-0;
    if(progbar==TRUE){  prog<-winProgressBar("Extract intensity differences...",min=m,max=n);
						setWinProgressBar(prog, 0, title = "Extract intensity differences...", label = NULL);}
    for(k in m:n){
      if(progbar==TRUE){setWinProgressBar(prog, k, title = "Extract intensity differences...", label = NULL)}
      if(profileList[[7]][k,3]>1){
        ########################################################################
        # fill timeset #########################################################
        timeset[,4:length(timeset[1,])]<-0;
        timeset[,c(4,5)] <-.Call("fill_timeset",
                                as.numeric(timeset),
                                as.numeric(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),6]), # sampleIDs
                                as.numeric(profileList[[2]][(profileList[[7]][k,1]:profileList[[7]][k,2]),2]), # intensities
                                as.integer(length(timeset[,1])),
                                PACKAGE="enviMass"
                            )	
		########################################################################
        # interpolate & subtract blind,  #######################################
        ########################################################################
        # subtract blank #######################################################
        ########################################################################
        if(any(timeset[,4]>0)){ # any non-blind peak present?
			what<-1 # !=1 -> get raw output, i.e., peak series
			that<-.Call("meandel",
				as.numeric(timeset),
				as.integer(subit),
				as.numeric(subrat),
				as.numeric(numtime),
				as.integer(what),
				as.numeric(lags),
				as.numeric(threshold),
				as.integer(notrend),
				PACKAGE="enviMass"
			)
			if(what!=1){ 
				that<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),that);
				colnames(that)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))			
				# plot smoothed series ... and abort ###########################
				plot.new();
				plot.window(xlim=c(min(numtime),max(numtime)),ylim=c(min(that[that[,2]!=0,4]),max(that[that[,2]!=0,4])));
				box();axis(1);axis(2);
				title(xlab="Time",ylab="Intensity")
				points(numtime[that[,2]!=0],that[that[,2]!=0,4],col="red",type="l");
				for(i in 1:length(lags)){
					points(numtime[that[,2]!=0],that[that[,2]!=0,(5+i)],col="darkgrey",type="l");
					points(that[that[,(5+i+(length(lags)*2))]!=0,(5+i+(length(lags)*2))],that[that[,(5+i+(length(lags)*2))]!=0,(5+i+length(lags))],col="blue",pch=19);					
				}
				stop(" YOU wanted the smoothed series...\n")
			}else{		
				profileList[[7]][k,5]<-max(that[4,]); # current>abs.dev
				profileList[[7]][k,6]<-max(that[5,]); # global>abs.dev
				profileList[[7]][k,7]<-max(that[3,]); # abs.dev
				if(any(timeset[,5]>0)){ 			  # in blind?
					profileList[[7]][k,8]<-1 # in blind
					profileList[[7]][k,11]<-length(timeset[timeset[,5]!=0,5]) 	# number_peaks_blind
					profileList[[7]][k,13]<-mean(timeset[timeset[,5]!=0,5]) 	# mean_int_blind
				}else{
					profileList[[7]][k,8]<-0 # in blind
					profileList[[7]][k,11]<-0 # number_peaks_blind
					profileList[[7]][k,13]<-0 # mean_int_blind					
				}
				#############################################
				# Replicates: mean sample above mean blind ?
				if(any(timeset[,5]>0)){
					if( mean(timeset[,4][timeset[,2]!=0])>=(mean(timeset[,5][timeset[,3]!=0])*blindfold) ){
						profileList[[7]][k,9]<-1 
					}else{
						profileList[[7]][k,9]<-0 
					}
				}else{
					profileList[[7]][k,9]<-1
				}				
				#############################################
				profileList[[7]][k,10]<-length(timeset[timeset[,4]!=0,4]) 	# number_peaks_sample
				profileList[[7]][k,12]<-mean(timeset[timeset[,4]!=0,4]) 	# mean_int_sample	
			}
		}else{
			# profileList[[7]][k,7]<-() # abs.dev = not of interest for blind
			profileList[[7]][k,8]<-1 # only in blind = not above (single peak)
			profileList[[7]][k,9]<-0	
			profileList[[7]][k,10]<-0 # number_peaks_sample
			profileList[[7]][k,11]<-length(timeset[timeset[,5]!=0,5])# number_peaks_blind
			profileList[[7]][k,12]<-0 # mean_int_sample
			profileList[[7]][k,13]<-mean(timeset[timeset[,5]!=0,5]) # mean_int_blind								
		}
		profileList[[7]][k,17]<-timeset[leng,4]
		########################################################################
      }else{
		profileList[[7]][k,7]<-0;
		if( any( timeset[,3]== (profileList[[2]][profileList[[7]][k,1],6]) ) ){ # only in blind
			profileList[[7]][k,8]<-1 # in blind = not above (single peak)
			profileList[[7]][k,9]<-0	
			profileList[[7]][k,10]<-0 # number_peaks_sample
			profileList[[7]][k,11]<-1 # number_peaks_blind
			profileList[[7]][k,12]<-0 # mean_int_sample
			profileList[[7]][k,13]<-(profileList[[2]][profileList[[7]][k,1],2]) # mean_int_blind						
		}else{ # only in sample
			profileList[[7]][k,8]<-0 
			profileList[[7]][k,9]<-1 # not in blind = above (single peak)		
			profileList[[7]][k,6]<-(profileList[[2]][profileList[[7]][k,1],2])
			if(any(profileList[[2]][profileList[[7]][k,1],6]==latestID)){
			  profileList[[7]][k,5]<-(profileList[[2]][profileList[[7]][k,1],2])
			  profileList[[7]][k,17]<-(profileList[[2]][profileList[[7]][k,1],2])
			}
			profileList[[7]][k,10]<-1 # number_peaks_sample
			profileList[[7]][k,11]<-0 # number_peaks_blind
			profileList[[7]][k,12]<-(profileList[[2]][profileList[[7]][k,1],2]) # mean_int_sample
			profileList[[7]][k,13]<-0 # mean_int_blind							
		}
      }
    }
    if(progbar==TRUE){close(prog);}
	profileList[[1]][[4]]<-TRUE;
    ############################################################################
    return(profileList)

}
