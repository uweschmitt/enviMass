#' @title Reset workflow functionalities
#'
#' @export
#'
#' @description Reset downstream workflow functionalities
#'
#' @param down Name of node (logfile$Tasks_to_redo) to be altered
#' @param added Name of node (logfile$Tasks_to_redo) to be added
#' @param except Name of node (logfile$Tasks_to_redo) to be excluded - dangerous
#' @param single_file File-wise handler
#' @param check_node Only evaluate if concerned node is enabled in workflow?
#' @param
#' @param
#' @param
#' 
#' @details enviMass workflow function
#' 

workflow_set<-function(down,added=FALSE,except=FALSE,single_file=FALSE,check_node=FALSE,...){
	
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in workflow_set.r!")}
	########################################################################################
	if(!is.logical(added) & !is.logical(except)){
		if((any(!is.na(match(added,except)))) || (any(!is.na(match(except,added))))){
			stop("workflow_set: added or except but not both.")
		}
	}
	if(is.na(match(down,names(logfile$Tasks_to_redo)))){
		stop(paste("workflow_set: unknown down argument",down))
	}
	if(length(down)>1){
		stop("workflow_set: which down?")
	}
	########################################################################################
	# leave funtion if check_node=TRUE (=parameters changed) but node not run ##############
	if(check_node){
		if(logfile$workflow[names(logfile$workflow)==down]=="no"){
			return(NULL);
		}
		if(logfile$workflow[names(logfile$workflow)==down]=="FALSE"){
			return(NULL);
		}		
	}
	########################################################################################
	depend<-logfile[[11]]
	########################################################################################
	# retrieve tasks to redo ###############################################################
	# redefined dependencies:
	if(!is.logical(added)){
		for(i in 1:length(added)){
			depend[rownames(depend)==added[i],colnames(depend)==down]<-1
		}
	}
	if(!is.logical(except)){
		for(i in 1:length(except)){
			depend[rownames(depend)==except[i],colnames(depend)==down]<-0
		}
	}
	# direct dependencies:
	work_stream<-rownames(depend)[depend[,colnames(depend)==down]==1]
	# collect indirect downstream dependencies
	doit<-TRUE
	while(doit){
		doit<-FALSE
		new_stream<-work_stream
		for(i in 1:length(work_stream)){
			new_nodes<-rownames(depend)[depend[,colnames(depend)==work_stream[i]]==1]
			new_nodes<-new_nodes[is.na(match(new_nodes,new_stream))]
			if(length(new_nodes)>0){
				new_stream<-c(new_stream,new_nodes)
				doit<-TRUE
			}
		}
		work_stream<-unique(new_stream)
	}
	########################################################################################
	# update Tasks_to_redo #################################################################
	########################################################################################
	if(any(work_stream=="peakpicking")){
		logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="peakpicking"]<<-TRUE;
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(!single_file){
			measurements[,10]<-FALSE;
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		rm(measurements1);	
	}
	########################################################################################
	if(any(work_stream=="qc")){
		if(logfile$workflow[names(logfile$workflow)=="qc"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="qc"]<<-TRUE;			
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="qc"]<<-TRUE;
		}
	}	
	########################################################################################
	if(any(work_stream=="pattern")){ 
		logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="pattern"]<<-TRUE;
	}	
	########################################################################################
	if(any(work_stream=="recal")){
		if(logfile$workflow[names(logfile$workflow)=="recal"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="recal"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,12]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="recal"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,12]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);		
		}
	}	
	########################################################################################
	if(any(work_stream=="allign")){
		if(logfile$workflow[names(logfile$workflow)=="align"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="allign"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,13]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="allign"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,13]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);		
		}
	}	
	########################################################################################
	if(any(work_stream=="norm")){
		if(logfile$workflow[names(logfile$workflow)=="norm"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="norm"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,14]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="norm"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,14]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);		
		}
	}	
	########################################################################################
	if(any(work_stream=="replicates")){
		if(logfile$workflow[names(logfile$workflow)=="replicates"]=="yes"){	
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="replicates"]<<-TRUE;
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="replicates"]<<-TRUE;
		}
	}	
	########################################################################################
	if(any(work_stream=="blinds")){
		if(logfile$workflow[names(logfile$workflow)=="blinds"]=="yes"){	
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="blinds"]<<-TRUE;
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="blinds"]<<-TRUE;
		}
	}	
	########################################################################################
	if(any(work_stream=="profiling")){
		if(logfile$workflow[names(logfile$workflow)=="profiling"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="profiling"]<<-TRUE;				
		}
		if(!check_node){		
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="profiling"]<<-TRUE;	
		}
	}	
	########################################################################################
	if(any(work_stream=="IS_screen")){
		if(logfile$workflow[names(logfile$workflow)=="IS_screen"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_screen"]<<-TRUE;	
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_screen"]<<-TRUE;			
		}
	}	
	########################################################################################
	if(any(work_stream=="target_screen")){
		if(logfile$workflow[names(logfile$workflow)=="target_screen"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="target_screen"]<<-TRUE;	
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="target_screen"]<<-TRUE;			
		}
	}	
	########################################################################################
	if(any(work_stream=="IS_normaliz")){
		if(logfile$workflow[names(logfile$workflow)=="IS_normaliz"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_normaliz"]<<-TRUE;
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_normaliz"]<<-TRUE;		
		}
	}	
	########################################################################################
	if(any(work_stream=="trendblind")){
		if(logfile$workflow[names(logfile$workflow)=="trenddetect"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="trendblind"]<<-TRUE;
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="trendblind"]<<-TRUE;		
		}		
	}	
	########################################################################################
	if(any(work_stream=="LOD")){
		if(logfile$workflow[names(logfile$workflow)=="LOD"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="LOD"]<<-TRUE;
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="LOD"]<<-TRUE;		
		}		
	}	
	########################################################################################
	if(any(work_stream=="quantification")){
		if(logfile$workflow[names(logfile$workflow)=="quantification"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="quantification"]<<-TRUE;
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="quantification"]<<-TRUE;		
		}		
	}		
	########################################################################################	
	
	########################################################################################
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    logfile<<-logfile;
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in workflow_set.r!")}
	########################################################################################	
				
}
