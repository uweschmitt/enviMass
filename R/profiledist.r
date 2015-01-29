#' @title Plot trend distributions
#'
#' @description \code{profiledist} boxplots current and global trend distributions
#'
#' @param profileList profileList
#' @param ret Logical. Should outlier ranking be returned?
#'
#' @return 	The boxplot in grey shows the intensity distributions of all trends of concern, listing the IDs, mean masses (m/z) and mean retention time (RT)
#' of the profiles with the most intense trends on the right. Colored points are used to elucidate the current 
#' trend intensities from the latest input file. Red dots signify profiles with intensities in the outlier range 
#' of the global (past and latest) trend intensities; green dots symbolize those below. 
#' 
#' @details enviMass workflow function
#' 


profiledist<-function(profileList,ret=FALSE){

    ############################################################################
    if(!profileList[[1]][[4]]){stop("profileList not checked for temporal trends; aborted.")}
	if(is.logical(ret)){stop("ret must be logical")}
	############################################################################
	plot.new()
	par(mar=c(2,4,1,1))
	if(any(profileList[[7]][,6]!=0)){
		minit<-log10(min(c(profileList[[7]][profileList[[7]][,5]>0,5],profileList[[7]][profileList[[7]][,6]>0,6])));
		maxit<-log10(max(c(profileList[[7]][profileList[[7]][,5]>0,5],profileList[[7]][profileList[[7]][,6]>0,6])));
	}else{
		plot.window(xlim=c(0,10),ylim=c(0,10));
		text(5,5,labels="No past and/or current incidents detected")
		return("Nothing to plot!")
	}
	plot.window(xlim=c(-3,10),ylim=c(minit,maxit));
	title(ylab="log10 Intensity");
	mtext("Global trend intensities",side=1,at=-1,col="darkgrey",line=0)
	mtext("Current trend intensities",side=1,at=8,col="black",line=0)	
	# old incidents ############################################################
	this<-boxplot(log10(profileList[[7]][profileList[[7]][,6]>0,6]),plot=FALSE)
	atit<-rep(3.5,length(this$out));
	atit<-(atit+sample(seq(-0.25,0.25,0.002),length(atit),replace=TRUE))	
	points(atit,this$out, pch=19,cex=0.5,col="darkgrey")
	boxplot(log10(profileList[[7]][profileList[[7]][,6]>0,6]),add=TRUE,at=3.5,outline=FALSE,border="darkgrey",boxwex=1)
	box(col="white")
	axis(2);
	allit<-log10(profileList[[7]][profileList[[7]][,6]>0,6])
	allitID<-profileList[[7]][profileList[[7]][,6]>0,4]
	if(length(allit)>0){
		allitID<-allitID[order(allit,decreasing=FALSE)]
		allit<-allit[order(allit,decreasing=FALSE)]
		doat<-maxit
		i=length(allit)
		while((doat>minit)&(i>0)){
			mz<-mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==allitID[i],1]:profileList[[7]][profileList[[7]][,4]==allitID[i],2]),1])
			mz<-round(mz,digits=4)
			RT<-mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==allitID[i],1]:profileList[[7]][profileList[[7]][,4]==allitID[i],2]),3])
			RT<-round(RT,digits=1)
			text(-3.1,doat,labels=paste("ID: ",as.character(allitID[i])," - m/z: ",as.character(mz)," - RT: ",as.character(RT),sep=""),pos=4,col="darkgrey",cex=1)
			lines(c(1.6,3.1),c(doat,allit[i]),col="darkgrey")
			doat<-(doat-((maxit-minit)/15))
			i=i-1
		}
		below<-(allit[allit<this$stats[5,1]]);
		above<-(allit[allit>this$stats[5,1]]);
	}else{
		text(0,this$stats[3,1],labels="No past incidents",pos=4,col="darkgreen")
	}
	# new incidents ############################################################
	allit<-log10(profileList[[7]][profileList[[7]][,5]>0,5])

    scores <- allit - this$stats[5, 1];
	allitID<-profileList[[7]][profileList[[7]][,5]>0,4]

	if(length(allit)>0){
		allitID<-allitID[order(allit,decreasing=FALSE)]
		allit<-allit[order(allit,decreasing=FALSE)]
		scores<-scores[order(allit,decreasing=FALSE)]
		doat<-maxit
		i=length(allit)
		while((doat>minit)&(i>0)){
			mz<-mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==allitID[i],1]:profileList[[7]][profileList[[7]][,4]==allitID[i],2]),1])
			mz<-round(mz,digits=4)
			RT<-mean(profileList[[2]][(profileList[[7]][profileList[[7]][,4]==allitID[i],1]:profileList[[7]][profileList[[7]][,4]==allitID[i],2]),3])
			RT<-round(RT,digits=1)
			if(scores[i]>0){
				lines(c(4.35,5.7),c(allit[i],doat),col="red")
				text(5.85,doat,labels=paste("ID: ",as.character(allitID[i])," - m/z: ",as.character(mz)," - RT: ",as.character(RT),sep=""),pos=4,col="red",cex=1)
			}else{
				lines(c(4.35,5.7),c(allit[i],doat),col="darkgreen")
				text(5.85,doat,labels=paste("ID: ",as.character(allitID[i])," - m/z: ",as.character(mz)," - RT: ",as.character(RT),sep=""),pos=4,col="darkgreen",cex=1)
			}
			doat<-(doat-((maxit-minit)/15))
			i=i-1
		}
		below<-(allit[scores < 0]);
		above<-(allit[scores > 0]);
		points(rep(4,length(above)),above,pch=19,col="red",cex=1);
		points(rep(4,length(below)),below,pch=19,col="darkgreen",cex=1);
	}else{
		text(6,this$stats[3,1],labels="No current incidents",pos=4,col="darkgreen")
	}
    ############################################################################

	if(ret){
		ranking <- list(allitID, scores)
		names(ranking) <- c("profile_id", "score");
		return(ranking);
	}
}


