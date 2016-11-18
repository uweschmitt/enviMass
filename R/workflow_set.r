#' @title Reset workflow functionalities
#'
#' @export
#'
#' @description Reset downstream workflow functionalities
#'
#' @param down Name of node (logfile$Tasks_to_redo) to be altered
#' @param added Name of node (logfile$Tasks_to_redo) to be added - additionally
#' @param except Name of node (logfile$Tasks_to_redo) to be excluded - dangerous
#' @param single_file File-wise handler
#' @param check_node Only evaluate whether the concerned node is enabled in workflow?
#' @param single_node Apply single_file == TRUE resets of data.frame-measurements-entries only to down node?
#' 
#' @details enviMass workflow function
#' 

workflow_set<-function(down,added=FALSE,except=FALSE,single_file=FALSE,check_node=FALSE,single_node=FALSE,...){

	########################################################################################	
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in workflow_set.r!")}
	########################################################################################
	if(!is.logical(added) & !is.logical(except)){
		if((any(!is.na(match(added,except)))) || (any(!is.na(match(except,added))))){
			stop("workflow_set: added or except nodes but not both for a node.")
		}
	}
	if(is.na(match(down,names(logfile$Tasks_to_redo)))){
		stop(paste("workflow_set: unknown down argument",down))
	}
	if(length(down)>1){
		stop("workflow_set: from which node down? Please specify only one node at a time.")
	}
	########################################################################################
	
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
	depend<-logfile$workflow_depend
	diag(depend)<-1
	########################################################################################
	# retrieve tasks to redo ###############################################################
	# redefined dependencies:
	if(!is.logical(added)){
		for(i in 1:length(added)){
			depend[rownames(depend)==added[i],colnames(depend)==down]<-1
		}
	}
	# direct dependencies:
	work_stream<-rownames(depend)[depend[,colnames(depend)==down]>0]
	# collect indirect downstream dependencies
	doit<-TRUE
	while(doit){
		doit<-FALSE
		new_stream<-work_stream
		for(i in 1:length(work_stream)){
			new_nodes<-rownames(depend)[depend[,colnames(depend)==work_stream[i]]>0]
			new_nodes<-new_nodes[is.na(match(new_nodes,new_stream))]
			if(length(new_nodes)>0){
				new_stream<-c(new_stream,new_nodes)
				doit<-TRUE
			}
		}
		work_stream<-unique(new_stream)
	}
	# explicitly excluded dependencies:
	if(!is.logical(except)){
		if(any(work_stream==except)){
			those<-match(work_stream,except)	
			work_stream<-work_stream[is.na(those)]
		}
	}
	########################################################################################
	
	########################################################################################
	# update Tasks_to_redo #################################################################
	########################################################################################
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	while(length(work_stream)>0){
		if(logfile$workflow[names(logfile$workflow)==work_stream[1]]=="yes"){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==work_stream[1]]<<-TRUE;			
			if(	
				((!single_file) & (any(names(measurements)==work_stream[1]))) || 
				((single_node) & (work_stream[1]==down) & (any(names(measurements)==work_stream[1])))
			){
				measurements[,names(measurements)==work_stream[1]]<-FALSE;
			}
		}	
		if(!check_node){
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==work_stream[1]]<<-TRUE;
			if( 
				((!single_file) & (any(names(measurements)==work_stream[1]))) || 
				((single_node) & (work_stream[1]==down) & (any(names(measurements)==work_stream[1])))			
			){
				measurements[,names(measurements)==work_stream[1]]<-FALSE;
			}
		}
		work_stream<-work_stream[-1]
	}
	write.csv(measurements,file=file.path(logfile[[1]],"dataframes","measurements"),row.names=FALSE);
	rm(measurements);
	########################################################################################
		
	########################################################################################
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    logfile<<-logfile;
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in workflow_set.r!")}
	########################################################################################	
				
}
