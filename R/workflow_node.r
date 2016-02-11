#' @title Execute workflow node
#'
#' @export
#'
#' @description Executes workflow nodes, depending on their state
#'
#' @param name_workflow
#' @param name_summary
#' @param name_redo
#' @param name_output
#' @param path_do
#' @param path_undo
#' @param session
#' @param output
#' @param input
#' 
#' @details enviMass workflow function
#' 

workflow_node<-function(name_workflow,name_summary,name_redo,name_output,path_do=FALSE,path_undo=FALSE,session,output,input){

	if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected #1 in workflow_node.r at",name_output))}
	######################################################################################
	if(  
		logfile$workflow[names(logfile$workflow)==name_workflow]=="yes" && 
		!(
			(logfile$summary[(logfile$summary[,1]==name_summary),2]=="TRUE") &&
			(logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==name_redo]=="FALSE")
		)
	){
		if(!is.logical(path_do)){
			source(path_do,local=TRUE);
		}
		logfile$summary[(logfile$summary[,1]==name_summary),2]<<-"TRUE";
		logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==name_redo]<<-"FALSE";
		save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		summa[(logfile$summary[,1]==name_summary),2]<-"done"
		summa[(logfile$summary[,1]==name_summary),2]<<-"done"
		output$summa_html<-renderText(summary_html(summa));
		cat(paste(name_output,"done \n"));
		output$dowhat<-renderText(paste(name_output,"done ... wait"))
	######################################################################################
	}else{
		if(logfile$workflow[names(logfile$workflow)==name_workflow]=="no"){
			if(
				!(
					(logfile$summary[(logfile$summary[,1]==name_summary),2]=="FALSE") &
					(logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==name_redo]=="FALSE")
				)
			){
				if(!is.logical(path_undo)){
					source(path_undo,local=TRUE)
					summa[(logfile$summary[,1]==name_summary),2]<-"removed"
					summa[(logfile$summary[,1]==name_summary),2]<<-"removed"
					output$summa_html<-renderText(summary_html(summa));		
					cat(paste(name_output,"removed \n"));
					output$dowhat<-renderText(paste(name_output,"removed .... wait"))
				}
			}else{
				summa[(logfile$summary[,1]==name_summary),2]<-"skipped"
				summa[(logfile$summary[,1]==name_summary),2]<<-"skipped"
				output$summa_html<-renderText(summary_html(summa));		
				cat(paste(name_output,"removed \n"));
				output$dowhat<-renderText(paste(name_output,"removed .... wait"))			
			}
			logfile$Tasks_to_redo[names(logfile$Tasks_to_redo)==name_redo]<<-"FALSE";
			logfile$summary[(logfile$summary[,1]==name_summary),2]<<-"FALSE";
			save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
		}else{
			summa[(logfile$summary[,1]==name_summary),2]<-"ok"
			summa[(logfile$summary[,1]==name_summary),2]<<-"ok"
			output$summa_html<-renderText(summary_html(summa));
			cat(paste(name_output,"done before \n"));
			output$dowhat<-renderText(paste(name_output,"done before .... wait"))
		}
	}
	######################################################################################
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	if(any(ls()=="logfile")){stop(paste("\n illegal logfile detected #2 in workflow_node.r at",name_output))}	
}



