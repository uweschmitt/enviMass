#' @title Generate HTML table output of project state 
#'
#' @export
#'
#' @description Given a list of numerical IDs, \code{IDs} returns the smallest integer not part of that list. 
#'
#' @param IDs Vector of integers
#' 
#' @details enviMass workflow function
#' 

summary_html<-function(in_logfile_summary){
	
	use_summary<-in_logfile_summary[c(1,2,3,4,5,16,11,7,8,14,12,13,9,17,18,19,10,15),]
	use_summary[,1]<-as.character(use_summary[,1])
	use_summary[,2]<-as.character(use_summary[,2])			
	mat<-matrix(ncol=2,nrow=length(use_summary[,1]),"")
	mat[use_summary[,2]=="FALSE",2]<-"background-color: darkred;"
	mat[use_summary[,2]=="TRUE",2]<-"background-color: darkgreen;"	
	if(use_summary[1,2]=="TRUE"){mat[use_summary[,1]=="Data available?",2]<-"background-color: green;"}
	mat[use_summary[,2]=="ok",2]<-"background-color: lightgreen"
	mat[use_summary[,2]=="...",2]<-"background-color: orange"
	mat[use_summary[,2]=="done",2]<-"background-color: green"	
	mat[use_summary[,2]=="skipped",2]<-"background-color: grey"
	mat[use_summary[,2]=="removed",2]<-"background-color: blue"
	use_summary[,2]<-""
	use_summary[use_summary[,1]=="Data available?",1]<-"Files included"
	use_summary[use_summary[,1]=="peakpicking",1]<-"Peak picking"	
	use_summary[use_summary[,1]=="qc",1]<-"Quality control"	
	use_summary[use_summary[,1]=="pattern",1]<-"Compound patterns"	
	use_summary[use_summary[,1]=="recal",1]<-"Mass recalibration"				
	use_summary[use_summary[,1]=="norm",1]<-"Median intensity normalization "				
	use_summary[use_summary[,1]=="profiling",1]<-"Profile extraction"	
	use_summary[use_summary[,1]=="IS_normaliz",1]<-"IS-based intensity normalization"	
	use_summary[use_summary[,1]=="trendblind",1]<-"Trend detection"	
	use_summary[use_summary[,1]=="LOD",1]<-"LOD interpolation"
	use_summary[use_summary[,1]=="quantification",1]<-"Quantification"
	use_summary[use_summary[,1]=="blinds",1]<-"Blind subtraction #1"	
	use_summary[use_summary[,1]=="replicates",1]<-"Replicate filter"	
	use_summary[use_summary[,1]=="IS_screen",1]<-"IS screening"
	use_summary[use_summary[,1]=="target_screen",1]<-"Target screening"
	use_summary[use_summary[,1]=="IS_subtr",1]<-"IS subtraction"
	use_summary[use_summary[,1]=="target_subtr",1]<-"Target subtraction"
	use_summary[use_summary[,1]=="blind_subtr",1]<-"Blind subtraction #2"
	for(i in 1:length(mat[,1])){
		if(i%%2==0){
			mat[i,1]<-c("background-color: lightgrey")
		}else{
			mat[i,1]<-c("background-color: white")
		}
	}
	summary_report<-htmlTable::htmlTable(use_summary,
		rnames=FALSE,
		header=c("Tasks","Status"),
		align="left",ctable=TRUE,
		css.cell = mat
	)
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in summary_html.r!")}
	return(summary_report)
				
}
