#' @title Plot mass against intensity.
#'
#' @export
#'
#' @description \code{mass_int} Plot mass against intensity..
#'
#' @param profileList A profile list.
#' @param profileID ID of profile to be plotted
#'
#' @details  enviMass workflow function. 
#' 

mass_int<-function(
	profileList,
	profileID
){

	###############################################################################
    mz<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),1]))
	int<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),2]))
	plot(mz,int,pch=19,cex=.7,xlab="m/z",ylab="Intensity")
	abline(v=mean(mz),lwd=2,lty=2,col="lightblue")
	###############################################################################

}