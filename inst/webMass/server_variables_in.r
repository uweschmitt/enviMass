if(any(ls()=="logfile")){stop("\n illegal logfile detected #1 in server_variabels_in.r!")}
################################################################################
# get variable settings from logfile ###########################################
updateTextInput(session, "PWpath",value = as.character(logfile[[4]]))
# PEAK PICKING #################################################################
updateSliderInput(session, "peak_drtgap",value = as.numeric(logfile$parameters[[2]]))
updateSliderInput(session, "peak_dmzdens",value = as.numeric(logfile$parameters[[3]]))
updateNumericInput(session, "peak_minpeak",value = as.numeric(logfile$parameters[[4]]))
updateNumericInput(session, "peak_drtsmall2",value = as.numeric(logfile$parameters[[5]]))
updateSliderInput(session, "peak_drtfill",value = as.numeric(logfile$parameters[[6]]))
updateNumericInput(session, "peak_drtdens2",value = as.numeric(logfile$parameters[[7]]))
updateSliderInput(session, "peak_minint",value = as.numeric(logfile$parameters[[8]]))
updateSliderInput(session, "peak_SN",value = as.numeric(logfile$parameters[[9]]))
updateSliderInput(session, "peak_SB",value = as.numeric(logfile$parameters[[10]]))
updateSliderInput(session, "peak_recurs",value = as.numeric(logfile$parameters[[11]]))
updateSliderInput(session, "peak_ended",value = as.numeric(logfile$parameters[[12]]))
updateSliderInput(session, "peak_weight",value = as.numeric(logfile$parameters[[13]]))
updateNumericInput(session, "peak_maxint",value = as.numeric(logfile$parameters[[14]]))
# PROGRESS BAR #################################################################
updateSelectInput(session, "progressbar", "Show progress bar (Windows only)", choices = c("TRUE","FALSE"),  selected = as.character(logfile$parameters[[21]]))
# RESOLUTION ###################################################################
updateSelectInput(session, "resolution", "Instrument resolution:", choices =  names(resolution_list), selected= as.character(logfile$parameters[[22]])) 
# RECALIBRATION ################################################################              
updateSelectInput(session, "recal_what", "Reference compounds:", c("Internal standards","Target compounds","both"),selected= as.character(logfile$parameters[[30]]))
updateNumericInput(session, "recal_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[[31]]))                
updateSelectInput(session, "recal_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[32]]))
updateNumericInput(session, "recal_drt", "RT tolerance [s]", value = as.numeric(logfile$parameters[[33]]))   
# REPLICATE INTERSECTION ########################################################
updateNumericInput(session, "replicate_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[[15]]))   
updateSelectInput(session, "replicate_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"),  selected = as.character(logfile$parameters[[16]]))
#updateSelectInput(session, "replicate_recalib", "... and corrected by recalibration results (if available)", choices = c("TRUE"="TRUE","FALSE"="FALSE"),  selected = as.character(logfile$parameters[[17]]))
updateNumericInput(session, "replicate_delRT", "RT tolerance of a compound peaks across replicate samples [s]", value = as.numeric(logfile$parameters[[18]]))   
updateSliderInput(session, "replicate_IS_dInt",value = as.numeric(logfile$parameters[[19]]))
# PROFILING ####################################################################
updateSliderInput(session, "prof_sets",value = as.numeric(logfile$parameters[[38]]))
updateTextInput(session, "upto_file","Up to file with ID:",value = as.character(logfile$parameters$upto_file)) # too risky!
updateNumericInput(session, "prof_dmz", "Peak-to-peak deviation within profiles: m/z tolerance ...", value = as.numeric(logfile$parameters[[39]]))                
updateSelectInput(session, "prof_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[40]]))
updateNumericInput(session, "prof_drt", "RT tolerance [s]", value = as.numeric(logfile$parameters[[41]]))                
updateSelectInput(session, "prof_select", "Omit files with table entry prof?=FALSE from profiling?",choices = c("TRUE"="TRUE","FALSE"="FALSE"),selected= as.character(logfile$parameters[[90]]))
# TREND ########################################################################
updateRadioButtons(session, "trend_blind", label = "Additional blind interpolation and subtraction per profile?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$parameters[[36]]) )
updateTextInput(session, "trend_lags","Time lags of trends [days], comma-separated:",value = as.character(logfile$parameters[[34]]))
updateNumericInput(session, "trend_thres", "Trend vs. mean+variance intensity threshold:", value = as.numeric(logfile$parameters[[35]]))   
updateSelectInput(session, "notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), selected= as.character(logfile$parameters[[29]]))
# BLIND ########################################################################
updateNumericInput(session, "blind_fold", "...if intensity ratio sample/blind <", value = as.numeric(logfile$parameters[[37]]))   
updateNumericInput(session, "blind_dmz", "Peak mass deviation among measurements (+/-) ...", value = as.numeric(logfile$parameters[[82]]))  
updateSelectInput(session, "blind_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[83]]))
updateNumericInput(session, "blind_drt", "Peak deviation within profiles: RT tolerance [s]", value = as.numeric(logfile$parameters[[84]]))  
updateCheckboxInput(session,"subtract_pos_bydate", label =  "Subtract with the next blank/blind file preceding a sample by its date & time?", value=as.logical(logfile$parameters[[85]]))
updateCheckboxInput(session,"subtract_pos_byfile", label = "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", value=as.logical(logfile$parameters[[86]]))							
updateCheckboxInput(session,"subtract_neg_bydate", label = "Subtract with the next blank/blind file preceding a sample by its date & time?", value=as.logical(logfile$parameters[[87]]))
updateCheckboxInput(session,"subtract_neg_byfile",label =  "Additional non-sample files to subtract each sample file with (i.e. not preceding by date only), choose file ID:", value=as.logical(logfile$parameters[[88]]))				
updateRadioButtons(session, "blind_omit", label = "Remove affected peaks?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$parameters[[89]]) )
# IS SCREENING #########################################################
updateNumericInput(session, "screen_IS_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", value = as.numeric(logfile$parameters[42]))   
updateNumericInput(session, "screen_IS_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", value = as.numeric(logfile$parameters[43])) 
updateNumericInput(session, "screen_IS_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[45])) 
updateSelectInput(session, "screen_IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[46]]))
updateSliderInput(session, "screen_IS_dInt", value = as.numeric(logfile$parameters[[47]])) # Intensity tolerance %
updateNumericInput(session, "screen_IS_Intcut", "Lower intensity threhold", value = as.numeric(logfile$parameters[48])) 
updateNumericInput(session, "screen_IS_w1", "Cutoff score [0,1]", value = as.numeric(logfile$parameters[49])) 
# TARGET SCREENING ######################################################
updateNumericInput(session, "screen_target_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", value = as.numeric(logfile$parameters[55]))   
updateNumericInput(session, "screen_target_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", value = as.numeric(logfile$parameters[56])) 
updateNumericInput(session, "screen_target_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[58])) 
updateSelectInput(session, "screen_target_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[59]]))
updateSliderInput(session, "screen_target_dInt", value = as.numeric(logfile$parameters[[60]])) # Intensity tolerance %
updateNumericInput(session, "screen_target_Intcut", "Lower intensity threhold", value = as.numeric(logfile$parameters[61])) 
updateNumericInput(session, "screen_target_w1", "Cutoff score [0,1]", value = as.numeric(logfile$parameters[62]))  
# IS-based NORMALIZATION ###############################################		
updateSliderInput(session, "profnorm_cover_files", value = as.numeric(logfile$parameters[[70]])) # Minimum percentage of files covered by each IS profile %
updateNumericInput(session, "profnorm_cover_isccount", "Minimum number of IS profiles", value = as.numeric(logfile$parameters[71])) 
updateCheckboxInput(session, "profnorm_use_blank", label = "Show median deviation of blank/blind profiles?",value = as.logical(logfile$parameters[72]))
updateCheckboxInput(session, "profnorm_use_blank_sample", label = "Use subsampling",value = as.logical(logfile$parameters[73]))
updateNumericInput(session, "profnorm_use_blank_samplecount", "Number of blank/blind profiles in subsample", value = as.numeric(logfile$parameters[74])) 
updateCheckboxInput(session, "profnorm_use_nonblank", label = "Show median deviation of sample (i.e., non-blank) profiles?",value = as.logical(logfile$parameters[75]))
updateCheckboxInput(session, "profnorm_use_nonblank_sample", label = "Use subsampling",value = as.logical(logfile$parameters[76]))
updateNumericInput(session, "profnorm_use_nonblank_samplecount", "Number of sample profiles in subsample", value = as.numeric(logfile$parameters[77])) 
updateSliderInput(session, "profnorm_threshold", value = as.numeric(logfile$parameters[[78]])) 





# COMPONENTIZATION ######################################################


#updateCheckboxGroupInput(session, "isos", "Select relevant isotopes:", choices =  as.character(unique(isotopes[,1])),selected=c("C","S","Br","Cl"))

# HOMOLOGUE SERIES ######################################################

	

# adducts ###############################################################
updateCheckboxGroupInput(session, "adducts_pos", "Positive ions:", choices =  as.character(adducts[adducts[,6]=="positive",1]),selected=as.character(logfile[[7]]))
updateCheckboxGroupInput(session, "adducts_neg", "Negative ions:", choices =  as.character(adducts[adducts[,6]=="negative",1]),selected=as.character(logfile[[8]]))                              

#########################################################################
# subtraction files, positive: ##########################################
if(any( (measurements[,1]!="-") & (measurements[,4]=="positive") & (measurements[,3]!="sample"))){
	IDs_pos<-measurements[
		(measurements[,4]=="positive") & (measurements[,3]!="sample")
	,1]
	names_pos<-measurements[
		(measurements[,4]=="positive") & (measurements[,3]!="sample")
	,2]
	IDs_pos<-paste(IDs_pos,names_pos,sep=" - ")
	if(any(logfile[[13]]!="FALSE")){
		select_pos<-logfile[[13]]
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
if(any( (measurements[,1]!="-") & (measurements[,4]=="negative") & (measurements[,3]!="sample"))){
	IDs_neg<-measurements[
		(measurements[,4]=="negative") & (measurements[,3]!="sample")
	,1]
	names_neg<-measurements[
		(measurements[,4]=="negative") & (measurements[,3]!="sample")
	,2]
	IDs_neg<-paste(IDs_neg,names_neg,sep=" - ")
	if(any(logfile[[14]]!="FALSE")){
		select_neg<-logfile[[14]]
		select_neg<-select_neg[select_neg!="FALSE"]
		select_neg<-select_neg[!is.na(match(select_neg,IDs_neg))]
	}else{
		select_neg<-NULL
	}
	updateCheckboxGroupInput(session,inputId="files_neg_select_subtract", label="", choices=IDs_neg, selected = select_neg)
}
#########################################################################

#################################################################################
# WORKFLOW SETTINGS #############################################################
updateRadioButtons(session, "qc", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="qc"]) )
updateRadioButtons(session, "recal", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="recal"]) )
#updateRadioButtons(session, "RTalign", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[3]]))
updateRadioButtons(session, "intnorm", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="norm"]) )
updateRadioButtons(session, "profnorm", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="IS_normaliz"]) )
updateRadioButtons(session, "profiling", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="profiling"]) )
updateRadioButtons(session, "trenddetect", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="trendblind"]) )
updateRadioButtons(session, "replicates", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="replicates"]) )
updateRadioButtons(session, "blind_filter", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="blinds"]) )
updateRadioButtons(session, "LOD_interpol", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="LOD"]) )
updateRadioButtons(session, "quantif", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="quantification"]) )
updateRadioButtons(session, "screen_IS", label = "Screen internal standards?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="IS_screen"]) )
updateRadioButtons(session, "screen_target", label = "Screen targets/suspects?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="target_screen"]) )
updateRadioButtons(session, "subtr_IS", label = "Subtract internal standards?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="IS_subtr"]) )
updateRadioButtons(session, "subtr_target", label = "Subtract internal targets?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="target_subtr"]) )
updateRadioButtons(session, "subtr_blind", label = "Subtract?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile$workflow[names(logfile$workflow)=="blind_subtr"]) )
################################################################################
################################################################################
if(any(ls()=="logfile")){stop("\n illegal logfile detected #2 in server_variabels_in.r!")}









