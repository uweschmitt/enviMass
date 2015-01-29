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
# PROFILING ####################################################################
updateSliderInput(session, "prof_sets",value = as.numeric(logfile$parameters[[38]]))
updateTextInput(session, "upto_file","Up to file with ID:",value = as.character(logfile$parameters$upto_file)) # too risky!
updateNumericInput(session, "prof_dmz", "Peak-to-peak deviation within profiles: m/z tolerance ...", value = as.numeric(logfile$parameters[[39]]))                
updateSelectInput(session, "prof_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[40]]))
updateNumericInput(session, "prof_drt", "RT tolerance [s]", value = as.numeric(logfile$parameters[[41]]))                
# TREND ########################################################################
updateTextInput(session, "trend_lags","Time lags of trends [days], comma-separated:",value = as.character(logfile$parameters[[34]]))
updateNumericInput(session, "trend_thres", "Trend vs. mean+variance intensity threshold:", value = as.numeric(logfile$parameters[[35]]))   
updateSelectInput(session, "notrend", "Do not show global trend - instead, report it as maximum intensity above blind", choices = c("TRUE"="TRUE","FALSE"="FALSE"), selected= as.character(logfile$parameters[[29]]))
# BLIND ########################################################################
updateSelectInput(session, "blind_do", "Run a blind subtraction...", c("yes"="yes","no"="no"),selected= as.character(logfile$parameters[[36]]))
updateNumericInput(session, "blind_fold", "...if intensity ratio sample/blind <", value = as.numeric(logfile$parameters[[37]]))   
# IS SCREENING #########################################################
updateNumericInput(session, "screen_IS_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", value = as.numeric(logfile$parameters[42]))   
updateNumericInput(session, "screen_IS_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", value = as.numeric(logfile$parameters[43])) 
updateNumericInput(session, "screen_IS_dRTblank", "RT tolerance of peaks in blank/blind relative to their expected RT [s]", value = as.numeric(logfile$parameters[44])) 
updateNumericInput(session, "screen_IS_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[45])) 
updateSelectInput(session, "screen_IS_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[46]]))
updateSliderInput(session, "screen_IS_dInt", value = as.numeric(logfile$parameters[[47]])) # Intensity tolerance %
updateNumericInput(session, "screen_IS_Intcut", "Lower intensity threhold", value = as.numeric(logfile$parameters[48])) 
updateNumericInput(session, "screen_IS_w1", "Score weight for mass matching", value = as.numeric(logfile$parameters[49])) 
updateNumericInput(session, "screen_IS_w2", "Score weight for relative intensity matching", value = as.numeric(logfile$parameters[50])) 
updateNumericInput(session, "screen_IS_w3", "Score weight for occurrence in blank/blind", value = as.numeric(logfile$parameters[51])) 
# TARGET SCREENING ######################################################
updateNumericInput(session, "screen_target_delRT", "RT tolerance of peaks in sample relative to their expected RT [s]", value = as.numeric(logfile$parameters[55]))   
updateNumericInput(session, "screen_target_dRTwithin", "RT tolerance of peaks within an isotope pattern [s]", value = as.numeric(logfile$parameters[56])) 
updateNumericInput(session, "screen_target_dRTblank", "RT tolerance of peaks in blank/blind relative to their expected RT [s]", value = as.numeric(logfile$parameters[57])) 
updateNumericInput(session, "screen_target_dmz", "m/z tolerance ...", value = as.numeric(logfile$parameters[58])) 
updateSelectInput(session, "screen_target_ppm", "... given in:", choices = c("ppm"="TRUE","absolute"="FALSE"), selected= as.character(logfile$parameters[[59]]))
updateSliderInput(session, "screen_target_dInt", value = as.numeric(logfile$parameters[[60]])) # Intensity tolerance %
updateNumericInput(session, "screen_target_Intcut", "Lower intensity threhold", value = as.numeric(logfile$parameters[61])) 
updateNumericInput(session, "screen_target_w1", "Score weight for mass matching", value = as.numeric(logfile$parameters[62])) 
updateNumericInput(session, "screen_target_w2", "Score weight for relative intensity matching", value = as.numeric(logfile$parameters[63])) 
updateNumericInput(session, "screen_target_w3", "Score weight for occurrence in blank/blind", value = as.numeric(logfile$parameters[64])) 
# IS-based NORMALIZATION ###############################################		
updateSliderInput(session, "profnorm_cover_files", value = as.numeric(logfile$parameters[[70]])) # Minimum percentage of files covered by each IS profile %
updateNumericInput(session, "profnorm_cover_isccount", "Minimum number of IS profiles", value = as.numeric(logfile$parameters[71])) 
updateCheckboxInput(session, "profnorm_use_blank", label = "Show median deviation of blank/blind profiles?",value = as.character(logfile$parameters[72]))
updateCheckboxInput(session, "profnorm_use_blank_sample", label = "Use subsampling",value = as.character(logfile$parameters[73]))
updateNumericInput(session, "profnorm_use_blank_samplecount", "Number of blank/blind profiles in subsample", value = as.numeric(logfile$parameters[74])) 
updateCheckboxInput(session, "profnorm_use_nonblank", label = "Show median deviation of sample (i.e., non-blank) profiles?",value = as.character(logfile$parameters[75]))
updateCheckboxInput(session, "profnorm_use_nonblank_sample", label = "Use subsampling",value = as.character(logfile$parameters[76]))
updateNumericInput(session, "profnorm_use_nonblank_samplecount", "Number of sample profiles in subsample", value = as.numeric(logfile$parameters[77])) 
updateSliderInput(session, "profnorm_threshold", value = as.numeric(logfile$parameters[[78]])) 





# COMPONENTIZATION ######################################################


#updateCheckboxGroupInput(session, "isos", "Select relevant isotopes:", choices =  as.character(unique(isotopes[,1])),selected=c("C","S","Br","Cl"))

# HOMOLOGUE SERIES ######################################################

	

# adducts ###############################################################
updateCheckboxGroupInput(session, "adducts_pos", "Positive ions:", choices =  as.character(adducts[adducts[,6]=="positive",1]),selected=as.character(logfile[[7]]))
updateCheckboxGroupInput(session, "adducts_neg", "Negative ions:", choices =  as.character(adducts[adducts[,6]=="negative",1]),selected=as.character(logfile[[8]]))                              



#########################################################################






#################################################################################
# WORKFLOW SETTINGS #############################################################
updateRadioButtons(session, "qc", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[1]]))
updateRadioButtons(session, "recal", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[2]]))
#updateRadioButtons(session, "RTalign", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[3]]))
updateRadioButtons(session, "intnorm", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[4]]))
updateCheckboxInput(session, "screen_IS_sam", value =  as.logical(as.character(logfile[[6]][11])))
updateCheckboxInput(session, "screen_target_sam", value =  as.logical(as.character(logfile[[6]][12])))
updateRadioButtons(session, "profnorm", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[15]]))
updateRadioButtons(session, "profiled", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[9]]))
updateRadioButtons(session, "trenddetect", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[10]]))
updateCheckboxInput(session, "Comp_isotop", value =  as.logical(as.character(logfile[[6]][7])))
updateCheckboxInput(session, "Comp_add", value = as.logical(as.character(logfile[[6]][8])))
updateCheckboxInput(session, "screen_IS_comp", value = as.logical(as.character(logfile[[6]][13])))
updateCheckboxInput(session, "screen_target_comp", value = as.logical(as.character(logfile[[6]][14])))
updateRadioButtons(session, "homol", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[16]]))
updateRadioButtons(session, "massdef", label = "Include?", choices =  c("yes"="yes","no"="no"), selected = as.character(logfile[[6]][[17]]))

################################################################################
################################################################################










