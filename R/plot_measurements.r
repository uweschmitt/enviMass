#' @title Overview plot of measurements
#'
#' @export
#'
#' @description Overview plot of measurements in a project
#'
#' @param logfile
#' @param ranges
#' 
#' @details enviMass workflow function for profile plotting
#' 


plot_measurements<-function(
	logfile,
	ranges_overview
){

    ############################################################################
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");	
	if(any(measurements[,"ID"]!="-")){
		dated<-measurements[,"Date"]
		timed<-measurements[,"Time"]
		datetime<-c()
		for(i in 1:length(timed)){
		  datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
		}
		atPOSIX<-as.POSIXct(datetime);
		attime2<-as.difftime(timed);
		atdate<-as.Date(dated);
		ord<-order(as.numeric(atdate),as.numeric(attime2));
		atPOSIXsort<-atPOSIX[ord];
		dated<-as.POSIXct(atPOSIXsort)
		############################################################################
		par(mar=c(2,4.5,.8,.2))
		plot.new()
		if(length(ranges_overview$x)==0){
			use_x_lim<-c(min(atPOSIX),max(atPOSIX))
			dated2<-pretty(dated)
		}else{
			use_x_lim<-ranges_overview$x	
			dated2<-dated[
				(dated>=ranges_overview$x[1])&
				(dated<=ranges_overview$x[2])
			]
			if(length(dated2)>0){
				dated2<-pretty(dated2)
			}
		}
		if(length(ranges_overview$y)==0){
			use_y_lim<-c(0,14)
		}else{
			use_y_lim<-ranges_overview$y
		}
		plot.window(
			xlim=use_x_lim,ylim=use_y_lim
		)
		text(min(atPOSIX),13,labels="+ Positive ionization",pos=4,cex=1,col="darkgray")
		text(min(atPOSIX),6,labels="- Negative ionization",pos=4,cex=1,col="darkgray")
		abline(h=7,col="black")
		mtext("Spiked",side=2,at=2,las=1,col="blue",cex=1)
		mtext("Calibration",side=2,at=3,las=1,col="red",cex=1)
		mtext("Blank/blind",side=2,at=4,las=1,col="darkgreen",cex=1)
		mtext("Sample",side=2,at=5,las=1,col="black",cex=1)
		mtext("Spiked",side=2,at=9,las=1,col="blue",cex=1)
		mtext("Calibration",side=2,at=10,las=1,col="red",cex=1)
		mtext("Blank/blind",side=2,at=11,las=1,col="darkgreen",cex=1)
		mtext("Sample",side=2,at=12,las=1,col="black",cex=1)	
		if((!is.null(ranges_overview$x))||(!is.null(ranges_overview$y))){
			mtext("Now zoomed in",side=3,col="gray",line=-1)
		}
		# time regions for calibration files
		if(any((measurements$Type=="calibration")&(measurements$Mode=="positive"))){ # add bars for calibration groups, positive!
			measurements2<-measurements[
				(measurements$Type=="calibration")&(measurements$Mode=="positive")
			,]
			for(i in 1:length(measurements2[,"ID"])){
				start_dated<-measurements2[i,"Date"]
				start_timed<-measurements2[i,"Time"]
				start_datetime<-paste(start_dated,start_timed,"CET",sep=" ")
				atPOSIX_start<-as.POSIXct(start_datetime)
				end_dated<-measurements2[i,"date_end"]
				end_timed<-measurements2[i,"time_end"]
				end_datetime<-paste(end_dated,end_timed,"CET",sep=" ")
				cat(end_datetime)
				atPOSIX_end<-as.POSIXct(end_datetime)
				lines(y=c(10,10),x=c(atPOSIX_start,atPOSIX_end),col="red")
				text(mean(c(atPOSIX_start,atPOSIX_end)),10,as.character(measurements2[i,]$tag2),pos=3,col="red")
				lines(y=c(9.8,10.2),x=c(atPOSIX_start,atPOSIX_start),col="red")
				lines(y=c(9.8,10.2),x=c(atPOSIX_end,atPOSIX_end),col="red")
			}
		}
		if(any((measurements$Type=="calibration")&(measurements$Mode=="negative"))){ # add bars for calibration groups, negative!
			measurements2<-measurements[
				(measurements$Type=="calibration")&(measurements$Mode=="negative")
			,]
			for(i in 1:length(measurements2[,"ID"])){
				start_dated<-measurements2[i,"Date"]
				start_timed<-measurements2[i,"Time"]
				start_datetime<-paste(start_dated,start_timed,"CET",sep=" ")
				atPOSIX_start<-as.POSIXct(start_datetime)
				end_dated<-measurements2[i,"date_end"]
				end_timed<-measurements2[i,"time_end"]
				end_datetime<-paste(end_dated,end_timed,"CET",sep=" ")
				atPOSIX_end<-as.POSIXct(end_datetime)
				lines(y=c(3,3),x=c(atPOSIX_start,atPOSIX_end),col="red")
				text(mean(c(atPOSIX_start,atPOSIX_end)),3,as.character(measurements2[i,]$tag2),pos=3,col="red")
				lines(y=c(2.8,3.2),x=c(atPOSIX_start,atPOSIX_start),col="red")
				lines(y=c(2.8,3.2),x=c(atPOSIX_end,atPOSIX_end),col="red")
			}
		}
		# add points for each file
		for(i in 1:length(measurements[,1])){
			at_x<-atPOSIX[i]
			if(
				(at_x<use_x_lim[1]) |
				(at_x>use_x_lim[2])
			){next;}
			if(measurements[i,]$Mode=="positive"){
				if(measurements[i,]$Type=="sample"){at_y<-12;use_col<-"black";with_pch<-15}
				if(measurements[i,]$Type=="blank"){at_y<-11;use_col<-"darkgreen";with_pch<-15}		
				if(measurements[i,]$Type=="calibration"){at_y<-10;use_col<-"black";with_pch<-0}		
				if(measurements[i,]$Type=="spiked"){at_y<-9;use_col<-"blue";with_pch<-15}	
				if(measurements[i,]$Type=="other"){at_y<-8;use_col<-"black";with_pch<-15}				
			}else{
				if(measurements[i,]$Type=="sample"){at_y<-5;use_col<-"black";with_pch<-15}		
				if(measurements[i,]$Type=="blank"){at_y<-4;use_col<-"darkgreen";with_pch<-15}	
				if(measurements[i,]$Type=="calibration"){at_y<-3;use_col<-"black";with_pch<-0}			
				if(measurements[i,]$Type=="spiked"){at_y<-2;use_col<-"blue";with_pch<-15}		
				if(measurements[i,]$Type=="other"){at_y<-1;use_col<-"black";with_pch<-15}			
			}
			if(measurements[i,]$include=="FALSE"){
				use_col<-"gray"
			}
			points(at_x,at_y,col=use_col,pch=with_pch,cex=.8)
		}
		#title(main="Select and double-click to zoom in, double-click again to zoom out.",cex.main=1)
		if(length(dated2)>0){
			axis(1,at=dated2,labels=dated2,col="grey",cex.axis=.9) # former at=dated
		}else{

		}
		############################################################################
	}else{
		plot.new()
		plot.window(xlim=c(0,1),ylim=c(0,1))
		text(.5,.5,labels="No files available")
	}
}


