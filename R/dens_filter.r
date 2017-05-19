#' @title Convert dates
#'
#' @description \code{dens_filter} gives estimates for ppm deviations and log10 cutoff intensities
#'
#' @param MSlist
#' @param plotit logical
#' @param n Number of descents
#' @param m Minimum peak count in a descents
#' @param endpoints Which endpoints of the greedy descent to plot? Those before ("b") or after ("a") the density filtering?
#' 
#' @details enviMass workflow function
#' 


dens_filter <-
function(MSlist,plotit=FALSE,n=2000,m=5,endpoints="a"){

  	##############################################################################
  	if(!length(MSlist)==8){stop("This is not an MSlist object")}
  	if(!MSlist[[1]][[1]]){stop("MSlist empty or invalid. Use readMSdata to upload raw .mzML data first.")}
  	if(!is.logical(plotit)){stop("Argument plotit musst be logical")}
  	##############################################################################
  	# order, structure ###########################################################
  	intens_order<-order(MSlist[["Scans"]][[2]][,"intensity"],decreasing=TRUE)
  	scans<-.Call("indexed",
    	as.numeric(MSlist[["Scans"]][[2]][,"RT"]),
    	as.integer(length(MSlist[["Scans"]][[1]])),
    	as.integer(3),
    	PACKAGE="enviMass" 
  	)
  	scan_num<-rep(0,dim(MSlist[["Scans"]][[2]])[1])
  	for(k in 1:(dim(scans)[1])){
    	scan_num[scans[k,1]:scans[k,2]]<-k
  	}
  	len<-length(MSlist[["Scans"]][[1]])
  	##############################################################################
  	# sample #####################################################################
  	if(plotit){
		plot.new()
		plot.window(ylim=c(log10(min(MSlist[["Scans"]][[2]][,"intensity"])),
			log10(max(MSlist[["Scans"]][[2]][,"intensity"]))),xlim=c(0,200))
		box();axis(1);axis(2)
		title(xlab="ppm-deviation from apex",ylab="log10 intensity")
	}
  	dead<-rep(NA,length(intens_order))
  	dead_ppm<-rep(NA,length(intens_order))
  	in_dead<-1
  	listit<-list()
  	at_list<-1
  	for(j in 1:length(intens_order)){
    	if(at_list>n){break} 
    	apex<-intens_order[j]
    	if(!is.na(dead[apex])){next}
    	at_scan<-scan_num[apex]
    	# RT-tail descend ##########################################################
    	doit<-TRUE
    	at_peak<-apex
    	found_peaks<-apex
    	while(doit & (at_scan<len)){
      		doit<-FALSE
      		at_scan<-(at_scan+1)
      		at_1<-(scans[at_scan,1]:scans[at_scan,2])
      		at_1<-at_1[abs(MSlist[["Scans"]][[2]][at_peak,"m/z"]-MSlist[["Scans"]][[2]][at_1,"m/z"])<.5]
      		if(length(at_1)){
        		at_1<-at_1[MSlist[["Scans"]][[2]][at_1,"intensity"]<=MSlist[["Scans"]][[2]][at_peak,"intensity"]]
        		if(length(at_1)){
          			at_1<-at_1[which.min(abs(MSlist[["Scans"]][[2]][at_peak,"m/z"]-MSlist[["Scans"]][[2]][at_1,"m/z"]))]
          			found_peaks<-c(found_peaks,at_1)
          			at_peak<-at_1
          			doit<-TRUE
        		}
      		} 
    	}
    	dead[found_peaks]<-1
    	lengp<-length(found_peaks)
    	if(lengp<m){next}
    	dead_ppm[in_dead:(in_dead+lengp-2)]<-(abs(MSlist[["Scans"]][[2]][found_peaks[-1],"m/z"]-MSlist[["Scans"]][[2]][found_peaks[1],"m/z"])/MSlist[["Scans"]][[2]][found_peaks[1],"m/z"]*1E6)
    	in_dead<-(in_dead+lengp-1)
    	listit[[at_list]]<-found_peaks
    	at_list<-(at_list+1)
		############################################################################	
		if(plotit){
			points(
				(abs(MSlist[["Scans"]][[2]][found_peaks[-1],"m/z"]-MSlist[["Scans"]][[2]][found_peaks[1],"m/z"])/MSlist[["Scans"]][[2]][found_peaks[1],"m/z"]*1E6),
				log10(MSlist[["Scans"]][[2]][found_peaks[-lengp],"intensity"]),
				pch=19,cex=.2,col="lightgrey"
			)
      if(endpoints=="b"){
  			points(
  				(abs(MSlist[["Scans"]][[2]][found_peaks[length(found_peaks)],"m/z"]-MSlist[["Scans"]][[2]][found_peaks[1],"m/z"])/MSlist[["Scans"]][[2]][found_peaks[1],"m/z"]*1E6),
  				log10(MSlist[["Scans"]][[2]][found_peaks[lengp],"intensity"]),
  				pch=19,cex=.5,col="black"
  			)
      }
		}	
  	}  
  	if(length(listit)<n){stop("\n ABORTED: sampling size n not reached")}
  	dead_ppm<-dead_ppm[1:(in_dead-1)]
  	##############################################################################
  	# get density cutoff #########################################################
  	dens<-density(dead_ppm[dead_ppm<200],adjust=.8)
  	x<-dens$x[dens$x>0]
  	y<-dens$y[dens$x>0]
  	for(k in 2:length(x)){
    	if((y[k-1]>=y[k])&(y[k]<=y[k+1])){
      		break
    	}
  	}
  	cut_ppm<-x[k]
  	if(plotit & (endpoints=="a")){
		  plot.window(xlim=c(0,200),ylim=c(min(y),max(y)))
		  points(x,y,type="l",col="darkgreen")
		  points(x[k],y[k],pch=19,cex=1.2,col="darkgreen")
  	}	
  	##############################################################################
  	# filter lists ###############################################################
  	store_y<-c()
  	store_x<-c()
  	store_y_min<-c()
  	store_x_min<-c()
  	#keep<-rep(FALSE,length(listit))
  	for(j in 1:length(listit)){
  		those<-((abs(MSlist[["Scans"]][[2]][listit[[j]],"m/z"]-MSlist[["Scans"]][[2]][listit[[j]][1],"m/z"])/
				MSlist[["Scans"]][[2]][listit[[j]][1],"m/z"]*1E6))
		if(!all(those[1:m]<cut_ppm)){
			next
		}else{
			for(k in 2:length(those)){
				if(those[k]>cut_ppm){
					break	
				}
			}
			k<-(k-1)
			store_x<-c(store_x,those[2:k])
  			store_y<-c(store_y,log10(MSlist[["Scans"]][[2]][listit[[j]][2:k],"intensity"]))
  			store_x_min<-c(store_x_min,those[k])
  			store_y_min<-c(store_y_min,log10(MSlist[["Scans"]][[2]][listit[[j]][k],"intensity"]))
			#keep[j]<-TRUE
		}
  	}
  	##############################################################################
  	# do some plotting ###########################################################
  	if(plotit){
      if(endpoints=="a"){
    		plot.window(ylim=c(log10(min(MSlist[[4]][[2]][,"intensity"])),
    			log10(max(MSlist[[4]][[2]][,"intensity"]))),xlim=c(0,200))
    		points(store_x_min,store_y_min,pch=19,cex=.5,col="black")
        abline(v=quantile(store_x_min, .95),col="red")
    		abline(h=quantile(store_y_min, .05),col="red")
      }   
      box();
  	}
  	##############################################################################  
  	return(c(quantile(store_x_min, .95),quantile(store_y_min, .05)))

}
