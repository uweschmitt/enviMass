
if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_variabels_in.r!")}

#########################################################################
# get parameter and workflow settings from logfile ######################
path<-"ui_mainPanel.r"
enviMass:::workflow_get(path,logfile,session)
#########################################################################
	
#########################################################################
# adducts ###############################################################
updateCheckboxGroupInput(session, "adducts_pos", selected=as.character(logfile$adducts_pos))
updateCheckboxGroupInput(session, "adducts_neg", selected=as.character(logfile$adducts_neg))
updateCheckboxGroupInput(session, "adducts_pos_group", selected=as.character(logfile$adducts_pos_group))
updateCheckboxGroupInput(session, "adducts_neg_group", selected=as.character(logfile$adducts_neg_group))                              
#########################################################################

#########################################################################
# PW Path ###############################################################
updateTextInput(session,inputId="PWpath",value=logfile$PW)
#########################################################################

#########################################################################
# subtraction files, positive: ##########################################
if(any( (measurements[,"ID"]!="-") & (measurements[,"Mode"]=="positive") & (measurements[,"Type"]!="sample"))){
	IDs_pos<-measurements[
		(measurements[,"Mode"]=="positive") & (measurements[,"Type"]!="sample")
	,1]
	names_pos<-measurements[
		(measurements[,"Mode"]=="positive") & (measurements[,"Type"]!="sample")
	,2]
	IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
	if(any(logfile[["Positive_subtraction_files"]]!="FALSE")){
		select_pos<-logfile[["Positive_subtraction_files"]]
		select_pos<-select_pos[select_pos!="FALSE"]
		# include changes from file additions / removals
		select_pos<-select_pos[!is.na(match(select_pos,IDs_pos))]
	}else{
		select_pos<-NULL
	}
	updateCheckboxGroupInput(session,inputId="files_pos_select_subtract", label="", choices=IDs_pos, selected = select_pos)
}
#########################################################################
# subtraction files, negative: ##########################################
if(any( (measurements[,"ID"]!="-") & (measurements[,"Mode"]=="negative") & (measurements[,"Type"]!="sample"))){
	IDs_neg<-measurements[
		(measurements[,"Mode"]=="negative") & (measurements[,"Type"]!="sample")
	,1]
	names_neg<-measurements[
		(measurements[,"Mode"]=="negative") & (measurements[,"Type"]!="sample")
	,2]
	IDs_neg<-paste(IDs_neg,names_neg,sep=" - ")
	if(any(logfile[["Negative_subtraction_files"]]!="FALSE")){
		select_neg<-logfile[["Negative_subtraction_files"]]
		select_neg<-select_neg[select_neg!="FALSE"]
		select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
	}else{
		select_neg<-NULL
	}
	updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
}
#########################################################################

if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_variabels_in.r!")}

