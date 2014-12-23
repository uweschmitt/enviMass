#' @title Mass density of a profile
#'
#' @export
#'
#' @description Plot mass density of a profile
#'
#' @param profileList A profile list.
#' @param profileID ID of profile to be plotted
#' @param bootstrap Add bootstrap confidene intervals?
#' @param boot_size Size of bootstrap sample
#' @param use_weights Weight by intensity?
#' 
#' @details enviMass workflow function for profile plotting
#' 

mass_dens<-function(
	profileList,
	profileID,
	bootstrap=FALSE,
	boot_size=200,
	use_weights=FALSE
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
	if(!is.logical(bootstrap)){stop("bootstrap must be logical")}
	if(!is.logical(use_weights)){stop("use_weights must be logical")}
    ############################################################################
    mz<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),1]))
	if(use_weights){
		wei<-(as.numeric(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==profileID,1]:profileList[[7]][profileList[[7]][,4]==profileID,2]),2]))
		wei<-(wei/sum(wei))
	}else{
		wei<-rep(1/length(mz),length(mz))
	}
	d<-density(	mz,
				bw="SJ",
				weights=wei
	)
	############################################################################
	if(bootstrap){
		ind<-seq(1,length(mz),1)
		fit2 <- replicate( boot_size,{ 
				use <- sample(ind, replace=TRUE); 
				x<-mz[use]
				wei_2<-wei[use]
				wei_2<-(wei_2/sum(wei_2))
				density(x, from=min(d$x), to=max(d$x),weights=wei_2)$y  
			}) 
		fit3 <- apply(fit2, 1, quantile, c(0.025,0.975) ) 
		plot(d,main="",lwd=2,ylim=c(0,max(max(fit3[1,]),max(fit3[2,]))));
		polygon( c(d$x, rev(d$x)), c(fit3[1,], rev(fit3[2,])), col='lightgrey', border=F) 
		lines(d);box()
	}else{
		plot(d,main="",lwd=2);
	}
    ############################################################################
	rug(mz,col="darkgreen")
	at_mean<-mean(mz)
	abline(v=at_mean,lwd=2,col="lightblue",lty=2)
	dn <- cumsum(d$y)/sum(d$y)
	li <- which(dn>=0.025)[1]
	ui <- which(dn>=0.975)[1]
	abline(v=d$x[c(li,ui)],col="red",lty=2)
	dif<-((d$x[ui]-d$x[li])/2)
	abline(v=d$x[li]+dif,lwd=2,col="darkblue")
	stringed<-paste(
			as.character(round((d$x[li]+dif),digits=5))," (+/-", as.character(round((dif/at_mean*1E6),digits=2))," ppm)"
	,sep="")
	############################################################################
	return(stringed)
	
}

