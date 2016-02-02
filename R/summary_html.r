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
	
	use_summary<-in_logfile_summary[c(1,2,3,4,5,11,7,8,9,10,12,13),]
	use_summary[,1]<-as.character(use_summary[,1])
	use_summary[,2]<-as.character(use_summary[,2])			
	mat<-matrix(ncol=2,nrow=length(use_summary[,1]),"")
	mat[use_summary[,2]=="FALSE",2]<-"background-color: darkred;"
	mat[use_summary[,2]=="TRUE",2]<-"background-color: darkgreen;"	
	if(use_summary[1,2]=="TRUE"){mat[use_summary[,1]=="Data available?",2]<-"background-color: green;"}
	mat[use_summary[,2]=="ok",2]<-"background-color: lightgreen;"
	mat[use_summary[,2]=="...",2]<-"background-color: orange"
	mat[use_summary[,2]=="done",2]<-"background-color: green"	
	mat[use_summary[,2]=="skipped",2]<-"background-color: grey"
	use_summary[,2]<-""
	use_summary[use_summary[,1]=="Data available?",1]<-"Files included"
	use_summary[use_summary[,1]=="Peak pick?",1]<-"Peak picking"	
	use_summary[use_summary[,1]=="QC?",1]<-"Quality control"	
	use_summary[use_summary[,1]=="Isotope pattern?",1]<-"Compound patterns"	
	use_summary[use_summary[,1]=="m/z recal.?",1]<-"Mass recalibration"				
	use_summary[use_summary[,1]=="Intensity norm.?",1]<-"Median intensity normalization "				
	use_summary[use_summary[,1]=="Profiled?",1]<-"Profile extraction"	
	use_summary[use_summary[,1]=="IS norm.?",1]<-"IS-based intensity normalization"	
	use_summary[use_summary[,1]=="Trend+Blind?",1]<-"Trend detection, blind subtraction "				
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
	return(summary_report)
				
}
