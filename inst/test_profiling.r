########################################################################
# On blind subtraction #################################################
cat("\n Checking profiling ...")
if(
	(logfile$workflow[names(logfile$workflow)=="profiling"])=="yes"
){

	####################################################################
	profileList_pos<-startprofiles(
		logfile,
		frac=FALSE,
		sets=as.numeric(logfile$parameters$prof_maxfiles),
		progbar=logfile$parameters$progressBar,
		ion_mode="positive",
		until=logfile$parameters$upto_file,
		selective=logfile$parameters$prof_select,
		types=c("sample","blank","spiked")
	)						
	profileList_pos<-agglomer(
		profileList_pos,
		dmass=(as.numeric(logfile$parameters$prof_dmz)+1),
		ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
		dret=(as.numeric(logfile$parameters$prof_drt)+10)
	)		
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,"include"]=="TRUE",]
	measurements<-measurements[measurements[,"Mode"]=="positive",]
	replicates<-measurements$tag3
	IDs<-measurements$ID
	profileList_pos<-partcluster(
		profileList=profileList_pos,
		dmass=as.numeric(logfile$parameters$prof_dmz),
		ppm=as.logical(as.character(logfile$parameters$prof_ppm)),
		dret=as.numeric(logfile$parameters$prof_drt),
		from=FALSE,
		to=FALSE,
		progbar=logfile$parameters$progressBar,
		plot_it=FALSE,
		replicates=replicates,
		IDs,
		with_test=TRUE
	)
	####################################################################	
	if(sum(profileList_pos[[2]][,"partitionIDs"]==0)){stop("issue_1")}	
	those<-match(
		profileList_pos[[2]][,"partitionIDs"],
		seq(length(profileList_pos["index_agglom"][[1]][,1]))
	)	
	if(any(is.na(those))){stop("issue_2")}
	for(i in 1:length(profileList_pos["index_agglom"][[1]][,1])){
		from<-profileList_pos["index_agglom"][[1]][i,1]
		to<-profileList_pos["index_agglom"][[1]][i,2]		
		if(!all(profileList_pos[[2]][from:to,"partitionIDs"]==i)){stop("issue_3")}
	}
	if(sum(profileList_pos[[2]][,"profileIDs"]==0)){stop("issue_4")}	
	those<-match(
		profileList_pos[[2]][,"profileIDs"],
		profileList_pos[[7]][,"profile_ID"]
	)	
	if(any(is.na(those))){stop("issue_5")}
	those<-match(
		profileList_pos[[7]][,"profile_ID"],
		profileList_pos[[2]][,"profileIDs"]
	)	
	if(any(is.na(those))){stop("issue_6")}	
	those<-match(
		profileList_pos[[2]][,"profileIDs"],
		seq(length(profileList_pos[[7]][,"profile_ID"]))
	)	
	if(any(is.na(those))){stop("issue_7")}	
	if(!all(profileList_pos[[7]][,"profile_ID"]==seq(length(profileList_pos[[7]][,"profile_ID"])))){stop("issue_5")}
	for(i in 1:length(profileList_pos["index_prof"][[1]][,1])){
		from<-profileList_pos["index_prof"][[1]][i,1]
		to<-profileList_pos["index_prof"][[1]][i,2]		
		if(!all(profileList_pos[[2]][from:to,"profileIDs"]==i)){stop("issue_8")}
	}
	####################################################################
	
}
cat("\n Checking profiling ... done.")
########################################################################


