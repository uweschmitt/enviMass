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

workflow_set<-function(down,added=FALSE,except=FALSE,single_file=FALSE,check_node=FALSE){
	
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in workflow_set.r!")}
	########################################################################################
	if(!is.logical(added) & !is.logical(except)){
		if((any(!is.na(match(added,except))))||(any(!is.na(match(except,added))))){
			stop("workflow_set: added or except but not both.")
		}
	}
	if(is.na(match(down,names(logfile$Tasks_to_redo)))){
		stop("workflow_set: unknown down argument.")
	}
	if(length(down)>1){
		stop("workflow_set: which down?")
	}
	########################################################################################
	# define workflow order of logfile$Tasks_to_redo by server.calculation.r ###############
	# dependencies must simply go after their parent node ################################## 
	# order here actually irrelevant, because calculation order set in server_calculation  #
	work_order<-c(1,2,8,3,5,4,9,6,12,10,11,15,13,7)
	work_names<-(names(logfile$Tasks_to_redo)[work_order])
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
	# define matrix of downstream workflow dependencies ####################################
	# requires only a definition of direct ones - inderect ones will be retrieved below ####
	# below specified in a row-wise fashion (but stored columnwise): #######################
	depend<-matrix(ncol=length(work_names),nrow=length(work_names),0)
	colnames(depend)<-work_names
	rownames(depend)<-work_names					# peakpick	QC	pattern	recal	allign	normalize	replicates	profiling	IS_screen	target_screen	norm_prof	trendblind	LOD		quantification
	depend[,colnames(depend)=="peakpick"]<-			c(0,		1,	0,		1,		1,		1,			1,			1,			1,			1,				1,			1,			1,		1)
	depend[,colnames(depend)=="QC"]<-				c(0,		0,	0,		1,		1,		0,			1,			1,			1,			1,				1,			1,			1,		1)
	depend[,colnames(depend)=="pattern"]<-			c(0,		0,	0,		1,		1,		0,			0,			0,			1,			1,				1,			0,			0,		1)
	depend[,colnames(depend)=="recal"]<-			c(0,		0,	0,		0,		1,		0,			1,			1,			1,			1,				1,			0,			0,		1)
	depend[,colnames(depend)=="allign"]<-			c(0,		0,	0,		0,		0,		0,			1,			1,			1,			1,				1,			0,			0,		0)
	depend[,colnames(depend)=="normalize"]<-		c(0,		0,	0,		0,		0,		0,			0,			1,			1,			1,				0,			1,			1,		1)
	depend[,colnames(depend)=="replicates"]<-		c(0,		0,	0,		0,		0,		0,			0,			1,			1,			1,				1,			1,			1,		1)
	depend[,colnames(depend)=="profiling"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			1,			1,				1,			1,			0,		1)
	depend[,colnames(depend)=="IS_screen"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			1,				1,			0,			0,		1)
	depend[,colnames(depend)=="target_screen"]<-	c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		1)
	depend[,colnames(depend)=="norm_prof"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			1,			0,		0)
	depend[,colnames(depend)=="trendblind"]<-		c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			0,			0,		1)
	depend[,colnames(depend)=="LOD"]<-				c(0,		0,	0,		0,		0,		0,			0,			0,			1,			1,				1,			1,			0,		1)
	depend[,colnames(depend)=="quantification"]<-	c(0,		0,	0,		0,		0,		0,			0,			0,			0,			0,				0,			1,			0,		0)	
	diag(depend)<-1
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
	if(any(work_stream=="peakpick")){
		logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="peakpick"]<<-TRUE;
		measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
		if(!single_file){
			measurements[,10]<-FALSE;
		}
		write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		rm(measurements1);	
	}
	########################################################################################
	if(any(work_stream=="QC")){
		if(logfile$workflow[names(logfile$workflow)=="qc"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="QC"]<<-TRUE;			
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="QC"]<<-TRUE;
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
	if(any(work_stream=="normalize")){
		if(logfile$workflow[names(logfile$workflow)=="norm"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="normalize"]<<-TRUE;			
			measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
			measurements[,14]<-FALSE;
			write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="normalize"]<<-TRUE;			
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
	if(any(work_stream=="profiling")){
		if(logfile$workflow[names(logfile$workflow)=="profiled"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="profiling"]<<-TRUE;				
		}
		if(!check_node){		
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="profiling"]<<-TRUE;	
		}
	}	
	########################################################################################
	if(any(work_stream=="IS_screen")){
		if(logfile$workflow[names(logfile$workflow)=="screen_IS"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_screen"]<<-TRUE;	
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="IS_screen"]<<-TRUE;			
		}
	}	
	########################################################################################
	if(any(work_stream=="target_screen")){
		if(logfile$workflow[names(logfile$workflow)=="screen_target"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="target_screen"]<<-TRUE;	
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="target_screen"]<<-TRUE;			
		}
	}	
	########################################################################################
	if(any(work_stream=="norm_prof")){
		if(logfile$workflow[names(logfile$workflow)=="profnorm"]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="norm_prof"]<<-TRUE;
		}
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)=="norm_prof"]<<-TRUE;		
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
