#' @title Get ID of file with latest date & time
#'
#' @export
#'
#' @description Given any measurements dataframe, extracts the ID of the newest file by date&time
#'
#' @param measurements Dataframe listing available enviMass files
#' 
#' @details enviMass workflow function. 
#' 

plot_components<-function(
		profileList,
		prof_IDs,
		links_profiles,
		what="relations",
		xlim=FALSE,
		ylim=FALSE,
		await_input=FALSE,
		skipit=TRUE,
		min_peaks=3
	){

	################################################################
	# make relation plot ###########################################
	if(any(what=="relations")){
		plot.new()
		if(xlim[1]==FALSE){
			xlim2<-c(min(profileList[["index_prof"]][prof_IDs,"mean_mz"]),max(profileList[["index_prof"]][prof_IDs,"mean_mz"]))
		}else{
			xlim2<-xlim
		}
		if(ylim[1]==FALSE){
			ylim2<-c(min(profileList[["index_prof"]][prof_IDs,"mean_RT"]),max(profileList[["index_prof"]][prof_IDs,"mean_RT"]))
		}else{
			ylim2<-ylim
		}		
		plot.window(xlim=xlim2,ylim=ylim2)			
		box();axis(1);axis(2);
		title(xlab="m/z",ylab="RT [s]")
		legend("topleft",
			legend=c("EIC","isot.","adduc.","homol."),
			col=c("red","blue","green","gray"),
			lwd=2.5
		)			
		# mark first profile
		points(
			profileList[["index_prof"]][prof_IDs[1],"mean_mz"],profileList[["index_prof"]][prof_IDs[1],"mean_RT"],
			pch=19,cex=4,col="yellow")	
		# add relations
		for(k in 1:length(prof_IDs)){
			fork<-profileList[["index_prof"]][prof_IDs[k],"links"]
			if(fork==0){next}
			##################################################################
			if(length(links_profiles[[fork]]$EIC[,1])>0){
				for(m in 1:length(length(links_profiles[[fork]]$EIC[,1]))){
					form<-links_profiles[[fork]]$EIC[m,1]
					if(any(prof_IDs==form)){colorit<-"red"}else{colorit<-"lightgrey";if(skipit){next}}
					if(links_profiles[[fork]]$EIC[m,"link counts"]>min_peaks){lty2<-1}else{lty2<-2}
					lines(
						x=c(profileList[["index_prof"]][form,"mean_mz"],profileList[["index_prof"]][prof_IDs[k],"mean_mz"]),
						y=c(profileList[["index_prof"]][form,"mean_RT"],profileList[["index_prof"]][prof_IDs[k],"mean_RT"]),
						col=colorit,lwd=1.5,lty=lty2
					)
				}
			}
			##################################################################
			if(length(links_profiles[[fork]]$isot[,1])>0){
				for(m in 1:length(links_profiles[[fork]]$isot[,1])){
					form<-links_profiles[[fork]]$isot[m,1]
					if(any(prof_IDs==form)){colorit<-"blue"}else{colorit<-"lightgrey";if(skipit){next}}
					if(links_profiles[[fork]]$isot[m,"link counts"]>min_peaks){lty2<-1}else{lty2<-2}
					lines(
						x=c(profileList[["index_prof"]][form,"mean_mz"],profileList[["index_prof"]][prof_IDs[k],"mean_mz"]),
						y=c(profileList[["index_prof"]][form,"mean_RT"],profileList[["index_prof"]][prof_IDs[k],"mean_RT"]),
						col=colorit,lwd=1.5,lty=lty2
					)
				}
			}				
			##################################################################
			if(length(links_profiles[[fork]]$adduc[,1])>0){
				for(m in 1:length(links_profiles[[fork]]$adduc[,1])){
					form<-links_profiles[[fork]]$adduc[m,1]
					if(any(prof_IDs==form)){colorit<-"green"}else{colorit<-"lightgrey";if(skipit){next}}
					if(links_profiles[[fork]]$adduc[m,"link counts"]>min_peaks){lty2<-1}else{lty2<-2}
					lines(
						x=c(profileList[["index_prof"]][form,"mean_mz"],profileList[["index_prof"]][prof_IDs[k],"mean_mz"]),
						y=c(profileList[["index_prof"]][form,"mean_RT"],profileList[["index_prof"]][prof_IDs[k],"mean_RT"]),
						col=colorit,lwd=1.5,lty=lty2
					)
				}
			}
			##################################################################		
			if(length(links_profiles[[fork]]$homol[,1])>0){
				for(m in 1:length(links_profiles[[fork]]$homol[,1])){
					form<-links_profiles[[fork]]$homol[m,1]
					if(any(prof_IDs==form)){colorit<-"black"}else{colorit<-"lightgrey";if(skipit){next}}
					if(links_profiles[[fork]]$homol[m,"link counts"]>min_peaks){lty2<-1}else{lty2<-2}
					lines(
						x=c(profileList[["index_prof"]][form,"mean_mz"],profileList[["index_prof"]][prof_IDs[k],"mean_mz"]),
						y=c(profileList[["index_prof"]][form,"mean_RT"],profileList[["index_prof"]][prof_IDs[k],"mean_RT"]),
						col=colorit,lwd=1.5,lty=lty2
					)
				}
			}
			##################################################################						
		}
		# plot profile points
		points(
			profileList[["index_prof"]][prof_IDs,"mean_mz"],profileList[["index_prof"]][prof_IDs,"mean_RT"],
			pch=19,cex=.8,col="black")	
	}
	################################################################
	if(any(what=="profiles")){
		those<-prof_IDs
		samp<-as.numeric(profileList[[4]])
		samp<-samp[order(samp)]
		maxint<-0
		for(j in 1:length(those)){
			int_max<-max(profileList[[2]][profileList[[7]][those[j],1]:profileList[[7]][those[j],2],2])
			if(int_max>maxint){
				maxint<-int_max
				this<-j
			}
		}
		plot.new()
		plot.window(
			xlim=c(min(samp),max(samp)),ylim=c(0,maxint)
		)
		box();axis(1);axis(2);
		title(xlab="Samples",ylab="Intensity")
		for(j in 1:length(those)){
			int<-profileList[[2]][profileList[[7]][those[j],1]:profileList[[7]][those[j],2],2]				
			sam<-profileList[[2]][profileList[[7]][those[j],1]:profileList[[7]][those[j],2],6]
			all_int<-rep(0,length(samp))
			at<-match(sam,samp)
			all_int[at]<-int
			lines(samp[all_int!=0],all_int[all_int!=0],col="grey")
			if(j==this){
				lines(samp[all_int!=0],(all_int)[all_int!=0],col="red",lwd=2.5)
			}
		}	
	}
	################################################################
	if(await_input){invisible(readline(prompt="Press [enter] to continue"))}
	
}
