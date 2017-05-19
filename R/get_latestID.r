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

get_latestID<-function(measurements){

	################################################################
    measurements<-measurements[measurements[,"include"]=="TRUE",]
	measurements<-measurements[((measurements[,"Type"]=="sample") |(measurements[,"Type"]=="blank")),]
	################################################################
    dated<-measurements[,"Date"]
    timed<-measurements[,"Time"]
    datetime<-c()
    for(i in 1:length(timed)){
      datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
    }
	atPOSIX<-as.POSIXct(datetime);
	sampleID<-measurements[,"ID"];
	ord<-order(atPOSIX,decreasing=TRUE);
	sampleID<-sampleID[ord];
	latestID<-sampleID[1]
	################################################################
	return(latestID)
	
}
