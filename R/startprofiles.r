#' @title Initialize profile list object
#'
#' @export
#'
#' @description \code{startprofiles} initializes a list object containing all peaks from the files in an enviMass project, associated metadata and placeholder.
#'
#' @param logfile logfile object of an enviMass project.
#' @param frac Numerical, \code{0<frac=<1}. Fraction of files to use; oldest are omitted.
#' @param sets Integer. Number of latest files to include.
#' @param progbar Logical. Should a progress bar be shown? Only for Windows.
#' @param ion_mode Character string, either "positive" or "negative".
#' @param until Integer, ID of file. All peaks of files up to the date of this file will be included.
#' @param selective Logical. Should only peaklist with measurements$profiled==TRUE be inluded?
#' @param types. File types to include in profiling, e.g., "sample", "blind", "calibration" or "spiked". For "spiked", all related files to subtract from are also included.
#'
#' @return profile list
#' 
#' @details  enviMass workflow function
#' 
#' @seealso \code{agglomer}, \code{partcluster} 

startprofiles<-function(
	logfile,
	frac=FALSE,
	sets=FALSE,
	progbar=FALSE,
	ion_mode="positive",
	until=FALSE,
	selective=FALSE,
	types=FALSE
){

    ############################################################################
    # check inputs #############################################################
    if(frac!=FALSE){if((frac<0)||(frac>=1)){stop("\n invalid frac argument. FALSE or 0<frac<=1")}}
	if((sets!=FALSE) & (!is.numeric(sets))){stop("\n sets must be FALSE or numeric; aborted.")}
	if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
	if(any(objects()=="peaklist")){rm(peaklist)}
	if(!is.logical(selective) & (selective!="FALSE" & selective!="TRUE")){stop("\n Argument selective must be logical")}
    ############################################################################
    # set up raw format ########################################################
    profiles<-list(0)
    profiles[[1]]<-data.frame(TRUE,FALSE,FALSE,FALSE)    # state
    colnames(profiles[[1]])<-c("peaks?","agglom?","profiling","trends?")
    profiles[[2]]<-0  # peaks
    profiles[[3]]<-0  # datetime
    profiles[[4]]<-0  # time
    profiles[[5]]<-0  # place
    profiles[[6]]<-0  # index_agglom
    profiles[[7]]<-0  # index_prof
    profiles[[8]]<-0  # parameters
    profiles[[9]]<-0  # sample type
    names(profiles)<-c("state","peaks","datetime","sampleID","place",
    "index_agglom","index_prof","parameters","type")
    ############################################################################
    # read in data #############################################################
    measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
    measurements<-measurements[measurements[,names(measurements)=="include"]=="TRUE",]
	if(selective=="TRUE"){
		measurements<-measurements[measurements[,names(measurements)=="profiled"]=="TRUE",]
	}
	measurements<-measurements[measurements[,names(measurements)=="Mode"]==ion_mode,]	
	# adjust time period, sort #################################################
    dated<-measurements[,"Date"]
    timed<-measurements[,"Time"]
    datetime<-c()
    for(i in 1:length(timed)){
      datetime<-c(datetime,paste(dated[i],timed[i],"CET",sep=" "))
    }
	atPOSIX<-as.POSIXct(datetime);
	sampleID<-measurements[,"ID"];
	locus<-measurements[,"Place"];
	typus<-measurements[,"Type"];
	ord<-order(atPOSIX,decreasing=TRUE);
	atPOSIX<-atPOSIX[ord];
	datetime<-datetime[ord];
	sampleID<-sampleID[ord];
	locus<-locus[ord];
	typus<-typus[ord];
	# filter by types ##########################################################
	if(types[1]!="FALSE"){
		remain<-rep(TRUE,length(sampleID))
		remain[is.na(match(typus,types))]<-FALSE
		if(any(types=="spiked") & any(measurements[,"Type"]=="spiked")){ # include subtraction files, too
			subtr_files<-measurements[measurements[,"ID"]==sampleID[typus=="spiked"],]$tag2
			remain[match(subtr_files,sampleID)]<-TRUE	
		}
		datetime<-datetime[remain]
		sampleID<-sampleID[remain]
		locus<-locus[remain]		
		typus<-typus[remain]	
		atPOSIX<-atPOSIX[remain]
	}
	# filter by latest & counts ################################################
	if(until!="FALSE" & any(sampleID==until)){
		untilPOSIX<-atPOSIX[sampleID==until]
		from<-(1:length(sampleID))[atPOSIX<=untilPOSIX][1]
	}else{
		if(until!="FALSE"){warning("\n problem in startprofiles: until used in conjunction with a removed sampleID. revise or debug?")}
		from<-1
	}
	if(sets!=FALSE){
		if((from+sets-1)>length(sampleID)){
			to<-length(sampleID)
		}else{
			to<-(from+sets-1)
		}
	}else{
		to<-length(sampleID)
	}
	datetime<-datetime[from:to];
	sampleID<-sampleID[from:to];
	locus<-locus[from:to];		
	typus<-typus[from:to];	
	# save to list #############################################################
    profiles[[3]]<-datetime;
    profiles[[4]]<-sampleID;
    profiles[[5]]<-locus;
    profiles[[8]]<-0;
    profiles[[9]]<-typus;	
    leng<-length(measurements[,"include"]);
    at<-c(0);
	############################################################################
	# get length of required matrix to store peaks #############################
	if(progbar==TRUE){  prog<-winProgressBar("Retrieve matrix length",min=0,max=length(sampleID));
                      setWinProgressBar(prog, 0, title = "Retrieve matrix length", label = NULL);
					  progi=0;}
	for(i in 1:leng){ 
		if(any(sampleID==as.numeric(measurements[i,1]))){
			if(selective=="TRUE"){
				if(measurements[i,names(measurements)=="profiled"]=="FALSE"){
					next;
				}
			}
			if(types[1]!="FALSE"){
				if(is.na(match(measurements[i,names(measurements)=="Type"],types))){
					next;
				}
			}
			if(progbar==TRUE){
				progi=progi+1;
				setWinProgressBar(prog, progi, title = "Retrieve matrix length", label = NULL)
			}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			if(logfile$parameters$blind_omit=="yes"){
				peaklist<<-peaklist[(peaklist[,colnames(peaklist)=="keep_2"]==1),,drop=FALSE]
			}
			peaklist<<-peaklist[(peaklist[,colnames(peaklist)=="keep"]==1),,drop=FALSE] # replicates
			if(length(peaklist[,1])==0){next}
			if(frac!=FALSE){
				at<-c(at+(floor(length(peaklist[,1])*frac)))
			}else{
				at<-c(at+length(peaklist[,1]))
			}
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
	}
	peaks<-matrix(nrow=(at),ncol=9,0)
    colnames(peaks)<-c("m/z","intensity","RT","peakIDs","links","sampleIDs","partitionIDs","profileIDs","in_blind")
    if(progbar==TRUE){close(prog);}
	da1<-c(1)
	############################################################################
	if(progbar==TRUE){  prog<-winProgressBar("Read peaklists",min=0,max=length(sampleID));
                      setWinProgressBar(prog, 0, title = "Read peaklists", label = NULL);
					  progi=0;}
    for(i in 1:leng){
		if(any(sampleID==as.numeric(measurements[i,1]))){
			if(selective==TRUE){
				if(measurements[i,names(measurements)=="profiled"]=="FALSE"){
					next;
				}
			}
			if(types[1]!="FALSE"){
				if(is.na(match(measurements[i,names(measurements)=="Type"],types))){
					next;
				}
			}
			if(progbar==TRUE){
				progi=progi+1;
				setWinProgressBar(prog, i, title = "Read peaklists", label = NULL);
			}
			load(file=file.path(logfile[[1]],"peaklist",as.character(measurements[i,1])),
				verbose=FALSE,envir=as.environment(".GlobalEnv"));
			if(logfile$parameters$blind_omit=="yes"){
				peaklist<<-peaklist[(peaklist[,colnames(peaklist)=="keep_2"]==1),,drop=FALSE]
			}
			peaklist<<-peaklist[(peaklist[,colnames(peaklist)=="keep"]==1),,drop=FALSE] # replicates
			if(length(peaklist[,1])==0){next}		
			if(frac!=FALSE){
				peaklist<<-peaklist[order(peaklist[,2],decreasing=TRUE),];
				that<-c((floor(length(peaklist[,1])*frac)));
			}else{
				that<-c(length(peaklist[,1]))
			}
			da2<-c(da1+that-1)				
			use_columns<-c("m/z","max_int","RT","peak_ID")
			if(logfile$workflow[names(logfile$workflow)=="recal"] == "yes"){
				use_columns[1]<-"m/z_corr"
			}
			if(logfile$workflow[names(logfile$workflow)=="align"] == "yes"){
				use_columns[3]<-"RT_corr"
			}
			peaks[da1:da2,]<-as.matrix(cbind( peaklist[1:that,use_columns], 	# must use the peakID as listed!
					rep(0,that),rep(as.numeric(measurements[i,1]),that),
					rep(0,that),rep(0,that),peaklist[1:that,colnames(peaklist)=="keep_2"])
			);
			da1<-c(da2+1);
			rm(peaklist,envir=as.environment(".GlobalEnv"))
		}
    }
	if(progbar==TRUE){close(prog);}	
    peaks<-peaks[order(peaks[,1],decreasing=FALSE),]
    profiles[[2]]<-peaks;
    rm(peaks)
    ############################################################################
    return(profiles)

}

