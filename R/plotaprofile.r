#' @title Plot profiles
#'
#' @export
#'
#' @description Plot profiles
#'
#' @param profileList A profile list.
#' @param profileID ID of profile to be plotted
#' @param logint Logical. Use log10 intensities?
#' @param blindsub Logical. Plot blind subtraction?
#' @param blindfold Numerical.
#' @param lags Vector of numericals.
#' @param threshold Numerical
#' @param add Logical. Plot parameter. Add to current plot?
#' @param textit Logical. Plot parameter.
#' @param simple Logical. Plot parameter.
#' @param colorit Logical. Plot parameter.
#' @param use_lwd Logical. Plot parameter.
#'
#' @return A dataset.
#' 
#' @details enviMass workflow function for profile plotting
#' 
#' @seealso  \code{startprofiles}, \code{agglomer}, \code{partcluster}, \code{intensup} 


plotaprofile<-function(
	profileList,
	profileID,
	logint=FALSE,
	blindsub=TRUE,
	blindfold=100,
	lags=c(5,14),
	threshold=3,
	add=FALSE,
	textit=TRUE,
	simple=FALSE,
	colorit=FALSE,
	use_lwd=FALSE
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
	if(blindsub!=FALSE){if(!is.numeric(blindfold) || (blindfold<0)){stop("Invalid blindfold argument; aborted.")}}
    if(blindsub!=FALSE){subit=1;subrat=blindfold;}else{subit=2;subrat=0;}
	if(!is.numeric(lags)){stop("lags argument must be numeric; aborted.")}
	if(colorit=="FALSE"){colorit<-"darkgreen"}
	if(use_lwd=="FALSE"){use_lwd<-2}
	############################################################################
    # set matrix to sort & store data from a profile ###########################
    atPOSIX<-profileList[[3]];
    sampletype<-profileList[[9]];
    sampleID<-profileList[[4]];
    atdate<-c();
    attime<-c();
    for(i in 1:length(atPOSIX)){
          atdate<-c(atdate, strsplit(atPOSIX[i]," ")[[1]][1]);
          attime<-c(attime, strsplit(atPOSIX[i]," ")[[1]][2]);
    }
    attime2<-as.difftime(attime);
    atdate<-as.Date(atdate);
    ord<-order(as.numeric(atdate),as.numeric(attime2),sampleID);
    atPOSIXsort<-atPOSIX[ord];
    atdate<-atdate[ord];
    attime2<-attime2[ord];
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
    latestID<-timeset[length(timeset[,2]),2]
	numtime<-(as.numeric(atdate)+as.numeric(attime2/24))
	colnames(timeset)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))	
    ############################################################################
    timeset[,4:5]<-0;
    timeset[,c(4,5)] <-.Call("fill_timeset",
                                as.numeric(timeset),
                                as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),6]), # sampleIDs
                                as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),2]), # intensities
                                as.integer(length(timeset[,1])),
                                PACKAGE="enviMass"
							)
 	what<-2 # !=1 -> get raw output, i.e., smoothed series
	notrend=1 # !=1 -> no global trends, only maximum intensity AFTER blind subtraction
	that1<-.Call("meandel",
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
   	that1<-matrix(nrow=length(atPOSIX),ncol=(5+(4*length(lags))),that1);
	colnames(that1)<-c("above blank?","sampleID","blankID","sample_int","blank_int",rep("lag_int",length(lags)),rep("del_int",length(lags)),rep("max_time",length(lags)),rep("blind_int",length(lags)))	
	what<-1 # !=1 -> get maximum intensities
	that2<-.Call("meandel",
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
	############################################################################
	# retrieve peak IDs ########################################################
	sampleID<-as.integer(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),6])
	peakID<-as.integer(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),4])
	atpeakID<-c()
	for(i in 1:length(timeset[,2])){
		if(timeset[i,2]!=0){
			if(any(sampleID==timeset[i,2])){
				atpeakID<-c(atpeakID,peakID[sampleID==timeset[i,2]])
			}else{
				atpeakID<-c(atpeakID,0)
			}
		}else{
			if(timeset[i,3]!=0){
				if(any(sampleID==timeset[i,3])){
					atpeakID<-c(atpeakID,peakID[sampleID==timeset[i,3]])
				}else{
					atpeakID<-c(atpeakID,0)
				}
			}
		}
	}	
	dataset<-data.frame(
		rev(as.character(atdate)),
		rev(as.character(attime)),
		rev(timeset[,2]),
		rev(timeset[,4]),
		rev(timeset[,3]),
		rev(timeset[,5]),
		rev(atpeakID),
		stringsAsFactors = FALSE
	)	
	dataset[,4]<-format(dataset[,4],scientific=TRUE,digits=2)
	dataset[,6]<-format(dataset[,6],scientific=TRUE,digits=2)
	dataset[,3]<-as.integer(dataset[,3])
	dataset[,5]<-as.integer(dataset[,5])
	dataset[,7]<-as.integer(dataset[,7])	
	names(dataset)<-c("Date","Time","ID sample","Intensity sample","ID blind","Intensity blind","ID peak")
	############################################################################
    dated<-as.POSIXct(atPOSIXsort)
    timelimit<-c(min(dated),max(dated))
    if(logint){
		if(!add){	
			plot.new()
			plot.window(xlim=c(0,10),ylim=c(0,10))
			if(textit){
				text(8.5,9.5,labels="Sample intensity",col="darkgreen",pos=4)
				text(8.5,9,labels="Blank intensity",col="red",pos=4)      
				text(8.5,8.5,labels=paste("mean m/z = ",round(mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),1]),digits=4),sep=""),col="black",pos=4)
				text(8.5,8,labels=paste("mean RT = ",round(mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),3]),digits=1),sep=""),col="black",pos=4)      
				#text(8.5,7.5,labels=paste("Partit. ID = ",round(unique(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),7]),digits=0),sep=""),col="black",pos=4)      
				text(8.5,7,labels=paste("Profile ID = ",round(unique(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),8]),digits=0),sep=""),col="black",pos=4)      	
			}
			plot.window(xlim=c(timelimit),ylim=c(0,max(log10(timeset[,4:5]))))
			axis(1,at=dated,labels=dated,col="grey",cex.axis=1)
			axis(2);
			box();
			title(xlab="Time",ylab="log10(intensity)")
			if(!simple){
				abline(h=log10(that2[6,]+(that2[3,]*threshold)),col="red",lty=2)
				abline(h=log10(that2[6,]),col="darkblue",lty=2)	
				points(dated[that1[,2]!=0],log10(that1[that1[,2]!=0,4]),type="l",col=colorit,lwd=use_lwd)
				points(dated[that1[,3]!=0],log10(that1[that1[,3]!=0,5]),type="l",col="red",lwd=2)
				for(i in 1:length(lags)){
					points(dated[that1[,2]!=0],log10(that1[that1[,2]!=0,(5+i)]),col="darkgrey",type="l");
				}
			}	
			points(dated[that1[,2]!=0],log10(that1[that1[,2]!=0,4]),type="l",col=colorit,lwd=use_lwd)			
		}else{	
			points(dated[that1[,2]!=0],log10(that1[that1[,2]!=0,4]),type="l",col=colorit,lwd=use_lwd)
		}
	}else{
		if(!add){	
			plot.new()
			plot.window(xlim=c(0,10),ylim=c(0,10))
			if(textit){
				text(8.5,9.5,labels="Sample intensity",col="darkgreen",pos=4)
				text(8.5,9,labels="Blank intensity",col="red",pos=4)      
				text(8.5,8.5,labels=paste("mean m/z = ",round(mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),1]),digits=4),sep=""),col="black",pos=4)
				text(8.5,8,labels=paste("mean RT = ",round(mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),3]),digits=1),sep=""),col="black",pos=4)      
				#text(8.5,7.5,labels=paste("Partit. ID = ",round(unique(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),7]),digits=0),sep=""),col="black",pos=4)      
				text(8.5,7,labels=paste("Profile ID = ",round(unique(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),8]),digits=0),sep=""),col="black",pos=4)      	
			}
			plot.window(xlim=c(timelimit),ylim=c(0,max((timeset[,4:5]))))
			axis(1,at=dated,labels=dated,col="grey",cex.axis=1)
			axis(2);
			box();
			title(xlab="Time",ylab="Intensity")
			if(!simple){
				abline(h=(that2[6,]+(that2[3,]*threshold)),col="red",lty=2)
				abline(h=(that2[6,]),col="darkblue",lty=2)	
				points(dated[that1[,3]!=0],(that1[that1[,3]!=0,5]),type="l",col="red",lwd=2)
				for(i in 1:length(lags)){
					points(dated[that1[,2]!=0],(that1[that1[,2]!=0,(5+i)]),col="darkgrey",type="l");
				}			
			}		
			points(dated[that1[,2]!=0],(that1[that1[,2]!=0,4]),type="l",col=colorit,lwd=use_lwd)
		}else{	
			points(dated[that1[,2]!=0],(that1[that1[,2]!=0,4]),type="l",col=colorit,lwd=use_lwd)
		}
    }
    ############################################################################
    return(dataset)

}


