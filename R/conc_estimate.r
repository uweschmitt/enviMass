


conc_estimate<-function(
	intensity,
	at_sample,
	cal_models,
	res_IS_X_screen,
	IS_table,
	get_all=TRUE
){

	#######################################################################################
	at_sample<-as.character(at_sample)
	res_names<-strsplit(names(res_IS_X_screen),"_")
	res_IS_names<-unlist(lapply(res_names, `[[`, 1))
	res_IS_adduct<-unlist(lapply(res_names, `[[`, 2))
	results<-list();at<-1;
	cat("\n")
	#######################################################################################
	for(i in 1:length(cal_models)){
		at_IS<-strsplit(names(cal_models)[i],"_")[[1]][2] # debug?
		at_adduct_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_adduct"] 	# get relevant IS adduct
		at_IS_entry<-which((res_IS_names==at_IS)&(res_IS_adduct==at_adduct_IS))
		if(length(at_IS_entry)!=1){next}
		at_IS_entry_sample<-which(names(res_IS_pos_screen[[at_IS_entry]])==at_sample)
		if(length(at_IS_entry_sample)!=1){next}
		at_peak_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_peak"] 	# get relevant IS peak
		rule_IS<-IS_table[IS_table[,"ID"]==at_IS,"Quant_rule"]	
		# get relevant intensities for IS peak ############################################	
		int_IS<-c()
		int_IS_rank<-c()
		for(k in 1:length(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]]) ){ # all matches for this file
			use_IS_peak<-which(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Peaks[,1]==at_peak_IS)
			if(length(use_IS_peak)!=1){next}									
			int_IS<-c(int_IS,
				(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
			)
			if(rule_IS=="most intense peak"){
				int_IS_rank<-c(int_IS_rank,
					(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$Intensity[use_IS_peak])
				)
			}
			if(rule_IS=="closest RT"){
				int_IS_rank<-c(int_IS_rank,
					1/abs(RT_IS-(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$RT[use_IS_peak]))
				)											
			}											
			if(rule_IS=="closest m/z"){
				int_IS_rank<-c(int_IS_rank,
					1/abs(res_IS_pos_screen[[at_IS_entry]][[at_IS_entry_sample]][[k]]$ppm[use_IS_peak])
				)											
			}													
		}
		if(length(int_IS)==0){next}
		# calculate concentrations ####################################
		get_conc<-c()
		for(b in 1:length(int_IS)){
			int_rat<-(intensity/int_IS[b])
			if(int_rat<cal_models[[i]]$low_bound){cat(",");next}		
			if(int_rat>cal_models[[i]]$up_bound){cat(",");next}											
			if(cal_models[[i]]$call=="resp ~ 0 + lin"){ # linear, 0-intercept
				new_conc<-(
					cal_models[[i]]$coefficients[[1]]*int_rat
				)
			}
			if(cal_models[[i]]$call=="resp ~ lin"){ # linear, with intercept
				new_conc<-(
					cal_models[[i]]$coefficients[[1]]+(cal_models[[i]]$coefficients[[2]]*int_rat)
				)
			}
			if(cal_models[[i]]$call=="resp ~ 0 + lin + quad"){ # quadratic, 0-intercept
				new_conc<-(
					(cal_models[[i]]$coefficients[[1]]*int_rat)+(cal_models[[i]]$coefficients[[2]]*(int_rat^2))
				)
			}
			if(cal_models[[i]]$call=="resp ~ lin + quad"){ # quadratic, 0-intercept
				new_conc<-(
					(cal_models[[i]]$coefficients[[1]])+
					(cal_models[[i]]$coefficients[[2]]*int_rat)+
					(cal_models[[i]]$coefficients[[3]]*(int_rat^2))														
				)
			}
			cat("*")							
			get_conc<-c(get_conc,as.vector(new_conc))				
		}
		# make entries into results list
		if(length(get_conc)>0){
			if(get_all){
				results[[at]]<-get_conc;
			}else{
				results[[at]]<-get_conc[1];			
			}
			names(results)[at]<-names(cal_models)[i]
			at<-(at+1);
		}
	}
	#######################################################################################
	return(results)
	#######################################################################################
	
}






