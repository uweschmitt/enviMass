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

summary_html<-function(use_summary){
	
	use_summary[,1]<-as.character(use_summary[,1])
	use_summary[,2]<-as.character(use_summary[,2])		
	len_init<-length(use_summary[,2])
	add<-length(unique(use_summary[,2]))
	states<-as.character(unique(use_summary[,2]))
	states[states=="FALSE"]<-"excluded / not run"
	states[states=="TRUE"]<-"included"
	states[states=="ok"]<-"done before"	
	states[states=="..."]<-"scheduled"		
	add_tab<-as.data.frame(matrix(nrow=(add+1),ncol=2,""),col.names=c("Tasks","Done?"),stringsAsFactors=FALSE)
	add_tab[,1]<-c("---------------------------------",states)
	names(add_tab)<-c("Tasks","Done?")
	use_summary<-rbind(
		use_summary,add_tab
	)
	mat<-matrix(ncol=2,nrow=length(use_summary[,1]),"")
	mat[use_summary[,2]=="FALSE",2]<-"background-color: grey"
	mat[use_summary[,2]=="TRUE",2]<-"background-color: green"	
	mat[use_summary[,2]=="ok",2]<-"background-color: lightgreen"
	mat[use_summary[,2]=="...",2]<-"background-color: orange"
	mat[use_summary[,2]=="done",2]<-"background-color: green"	
	mat[use_summary[,2]=="skipped",2]<-"background-color: grey"
	mat[use_summary[,2]=="removed",2]<-"background-color: lightblue"
	mat[use_summary[,1]=="not run",1]<-"background-color: grey"
	mat[use_summary[,1]=="included",1]<-"background-color: green"	
	mat[use_summary[,1]=="done before",1]<-"background-color: lightgreen"
	mat[use_summary[,1]=="scheduled",1]<-"background-color: orange"
	mat[use_summary[,1]=="done",1]<-"background-color: green"
	mat[use_summary[,1]=="skipped",1]<-"background-color: grey"	
	mat[use_summary[,1]=="removed",1]<-"background-color: lightblue"	
	
	use_summary[,2]<-""
	use_summary[use_summary[,1]=="peakpicking",1]<-"Peak picking"	
	use_summary[use_summary[,1]=="qc",1]<-"Quality control"	
	use_summary[use_summary[,1]=="pattern",1]<-"Compound patterns"	
	use_summary[use_summary[,1]=="recal",1]<-"Mass recalibration"				
	use_summary[use_summary[,1]=="norm",1]<-"Median intensity normalization "				
	use_summary[use_summary[,1]=="profiling",1]<-"Profile extraction"	
	use_summary[use_summary[,1]=="IS_normaliz",1]<-"IS-based intensity normalization"	
	use_summary[use_summary[,1]=="trendblind",1]<-"Trend detection"	
	use_summary[use_summary[,1]=="LOD",1]<-"LOD interpolation"
	use_summary[use_summary[,1]=="calibration",1]<-"Calibration"
	use_summary[use_summary[,1]=="quantification",1]<-"Quantification"	
	use_summary[use_summary[,1]=="blind",1]<-"Blind subtraction #1"	
	use_summary[use_summary[,1]=="replicates",1]<-"Replicate filter"	
	use_summary[use_summary[,1]=="IS_screen",1]<-"IS screening"
	use_summary[use_summary[,1]=="target_screen",1]<-"Target screening"
	use_summary[use_summary[,1]=="subtr",1]<-"Profile filtering"
	use_summary[use_summary[,1]=="isotopologues",1]<-"Isotopol. grouping"	
	use_summary[use_summary[,1]=="adducts",1]<-"Adduct grouping"	
	use_summary[use_summary[,1]=="homologues",1]<-"Homologue detection"		
	use_summary[use_summary[,1]=="recovery",1]<-"Recovery"			
	use_summary[use_summary[,1]=="components_files",1]<-"File componentization"		
	use_summary[use_summary[,1]=="components_profiles",1]<-"-"		

	for(i in 1:len_init){
		if(i%%2==0){
			mat[i,1]<-c("background-color: lightgrey")
		}else{
			mat[i,1]<-c("background-color: white")
		}
	}
	#mat[(len_init+1),]<-c("background-color: white")
	summary_report<-htmlTable::htmlTable(
		use_summary,
		rnames=FALSE,
		header=c("Tasks","Status"),
		align="left",ctable=TRUE,
		css.cell = mat
	)
	if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in summary_html.r!")}
	return(summary_report)
				
}
