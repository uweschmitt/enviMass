#' @title Update input/output values in script with logfile settings
#'
#' @export
#'
#' @description Updates shiny elements
#'
#' @param path String of file.r to evaluate
#' @param logfile enviMass logfile
#' 
#' @details enviMass workflow function to update shiny ui.r script with logfile parameters of workflow settings, Do be used inside a shiny session.
#' 

#path<-"D:/MS/R_packages/enviMass_devel/inst/webMass/ui_mainPanel.r"

workflow_get<-function(path,logfile,session){

	##########################################################################################
	if(!file.exists(path)){stop("\n path in workflow_get invalid; debug!")}
	if(!length(logfile$parameters)){stop("\n logfile problem in workflow_get invalid; no parameters available. Debug!")}
	if(!length(logfile$workflow)){stop("\n logfile problem in workflow_get invalid; no workflow nodes available. Debug!")}	
	##########################################################################################
	textit<-readLines(
		file.path(path)
	)
	##########################################################################################
	done_parameters<-rep(FALSE,length(logfile$parameters))
	done_workflow<-rep(FALSE,length(logfile$workflow))	
	##################################################################################
	types<-c("numericInput","sliderInput","selectInput","checkboxInput","checkboxGroupInput","textInput","radioButtons")
	types<-paste(types,"(",sep="")
	for(i in 1:length(textit)){
		for(j in 1:length(types)){
			if(grepl(types[j],textit[i], fixed = TRUE)){
				at<-gregexpr(pattern=types[j],textit[i], fixed = TRUE)
				for(k in 1:length(at[[1]][1])){
					# extract input name between quotes
					from<-(at[[1]][1][k]+nchar(types[j])+1)
					for(to in from:(from+10000)){
						if( substr(textit[i],to,to)=="\""){
							to<-(to-1)
							break;
						}
					}
					if(grepl("inputId",substr(textit[i],from-1,to+1))){
						stop("\n shiny input declared with explicit inputId - please revise; omit inputId!")
					}
					that<-substr(textit[i],from,to)
					if(any(names(logfile$parameters)==that)){
						if(types[j]=="numericInput("){
							do<-"updateNumericInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",value=as.numeric(logfile$parameters$",that,"))",sep="")
							))
						}
						if(types[j]=="sliderInput("){
							do<-"updateSliderInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",value=as.numeric(logfile$parameters$",that,"))",sep="")
							))					
						}						
						if(types[j]=="selectInput("){
							do<-"updateSelectInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",selected=as.character(logfile$parameters$",that,"))",sep="")
							))					
						}						
						if(types[j]=="checkboxInput("){
							do<-"updateCheckboxInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",value=as.logical(logfile$parameters$",that,"))",sep="")
							))					
						}					
						if(types[j]=="checkboxGroupInput("){
							do<-"updateCheckboxGroupInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",value=as.logical(logfile$parameters$",that,"))",sep="")
							))					
						}
						if(types[j]=="textInput("){
							do<-"updateTextInput"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",value=as.character(logfile$parameters$",that,"))",sep="")
							))					
						}
						if(types[j]=="radioButtons("){
							do<-"updateRadioButtons"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",selected=as.character(logfile$parameters$",that,"))",sep="")
							))												
						}
						cat("\n Updated parameter: ");cat(that);
						
						done_parameters[names(logfile$parameters)==that]<-TRUE
						
						
						next;
					}
					if(any(names(logfile$workflow)==that)){
						if(types[j]=="radioButtons("){
							do<-"updateRadioButtons"
							eval(parse(text=
								paste(do,"(session,\"",that,"\",selected = as.character(logfile$workflow[names(logfile$workflow)==\"",that,"\"]))",sep="")
							))												
						}else{
							stop("\n Only used radioButtons for workflow settings?")
						}
						done_workflow[names(logfile$workflow)==that]<-TRUE
						cat("\n Updated workflow: ");cat(that);
					}
				}		
			}
		}			
	}
	if(any(!done_workflow)){
		those<-paste(names(logfile$workflow)[!done_workflow],collapse=", ")
		cat("\n Fixed workflow steps: ");cat(those)
	}
	if(any(!done_parameters)){
		those<-paste(names(logfile$parameters)[!done_parameters],collapse=", ")
		cat("\n Fixed parameters: ");cat(those)
	}
	cat("\n")
	##########################################################################################
	return("Finished update")	
	##########################################################################################
	
}




























