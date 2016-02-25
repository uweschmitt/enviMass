#' @title Compile a data.frame from screening results
#'
#' @description Check measured pattern plausibility
#'
#' @param screened_listed
#' @param pattern
#' @param measurements_table
#' @param compound_table
#' @param cut_score
#' 
#' @details enviMass workflow function
#' 

	get_screening_results<-function(
		screened_listed,
		pattern,
		profileList,
		measurements_table,
		compound_table,
		cut_score
	){

		IDs<-as.numeric(measurements_table[,1]) 
		num_samples_all<-rep(0,length(screened_listed))
		num_blanks_all<-rep(0,length(screened_listed))
		max_score_sample_all<-rep(0,length(screened_listed))
		max_score_blank_all<-rep(0,length(screened_listed))
		num_peaks_sample_all<-rep(0,length(screened_listed))
		num_peaks_blank_all<-rep(0,length(screened_listed))			
		mean_int_ratio<-rep(0,length(screened_listed))
		IDed<-rep("")
		named<-rep("")
		adducted<-rep("")	
		at_len<-1		
		max_len<-10000
		at_matrix<-matrix(nrow=10000,ncol=7,0)
		colnames(at_matrix)<-c("m/z","log Intensity","RT","m/z deviation [ppm]","RT deviation","above_cutscore","Time sequence")
		for(i in 1:length(screened_listed)){
			IDed[i]<-strsplit(names(pattern)[i],"_")[[1]][1]
			named[i]<-compound_table[compound_table[,1]==strsplit(names(pattern)[i],"_")[[1]][1],2]
			adducted[i]<-strsplit(names(pattern)[i],"_")[[1]][2]
			num_samples<-(0)
			num_blanks<-(0)
			max_score_sample<-(0)
			max_score_blank<-(0)
			num_peaks_sample<-(0)
			num_peaks_blank<-(0)				
			centro_sample<-list()
			centro_blank<-list()
			for(j in 1:length(pattern[[i]][,1])){centro_sample[[j]]<-numeric(0);centro_blank[[j]]<-numeric(0);}
			if(length(screened_listed[[i]])>0){
				for(m in 1:length(screened_listed[[i]])){
					if(length(screened_listed[[i]][[m]])>0){
						is_sample<-(measurements_table[IDs==m,3]=="sample")					
						if(!is_sample){ # could still be doted or blind or ...
							is_blank<-(measurements_table[IDs==m,3]=="blank")
						}else{
							is_blank<-FALSE
						}
						if(!is_sample & !is_blank){next}
						max_score<-0
						max_num_peaks<-0
						for(k in 1:length(screened_listed[[i]][[m]])){
							if(length(screened_listed[[i]][[m]][[k]])>0){
								local_score<-0
								if(!is.na(screened_listed[[i]][[m]][[k]]$score_1)){
									local_score<-(local_score+screened_listed[[i]][[m]][[k]]$score_1)
								}
								if( (local_score>=1) || (is.na(screened_listed[[i]][[m]][[k]]$score_1)) ){
									if(!is.na(screened_listed[[i]][[m]][[k]]$score_2)){
										local_score<-(local_score+screened_listed[[i]][[m]][[k]]$score_2)
									}
								}									
								if(local_score>max_score){
									max_score<-local_score
								}
								if(length(screened_listed[[i]][[m]][[k]]$Peaks[,1])>max_num_peaks){
									max_num_peaks<-length(screened_listed[[i]][[m]][[k]]$Peaks[,1])
								}
								if(is_sample & (local_score>=cut_score)){
									for(d in 1:length(screened_listed[[i]][[m]][[k]][[1]][,1])){
										centro_sample[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]]<-c(
											centro_sample[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]],
											profileList[[2]][screened_listed[[i]][[m]][[k]][[1]][d,2],2]	
										)
									}
								}
								if(is_blank & (local_score>=cut_score)){
									for(d in 1:length(screened_listed[[i]][[m]][[k]][[1]][,1])){
										centro_blank[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]]<-c(
											centro_blank[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]],
											profileList[[2]][screened_listed[[i]][[m]][[k]][[1]][d,2],2]	
										)
									}
								}
								local_len<-length(screened_listed[[i]][[m]][[k]][[7]])
								if((at_len+local_len)>max_len){
									at_matrix<-rbind(
										at_matrix, 
										matrix(nrow=10000,ncol=7,0)
									)
									max_len<-(max_len+1000)
								}
								at_matrix[at_len:(at_len+local_len-1),1]<-screened_listed[[i]][[m]][[k]][[7]]
								at_matrix[at_len:(at_len+local_len-1),2]<-screened_listed[[i]][[m]][[k]][[8]]
								at_matrix[at_len:(at_len+local_len-1),3]<-screened_listed[[i]][[m]][[k]][[9]]
								at_matrix[at_len:(at_len+local_len-1),4]<-screened_listed[[i]][[m]][[k]][[4]]
								at_matrix[at_len:(at_len+local_len-1),5]<-screened_listed[[i]][[m]][[k]][[5]]
								if(local_score>=cut_score){
									at_matrix[at_len:(at_len+local_len-1),6]<-1
								}		
								at_matrix[at_len:(at_len+local_len-1),7]<-(
									as.numeric(as.Date(measurements_table[IDs==m,6]))+
									as.numeric(as.difftime(measurements_table[IDs==m,7])/24)
								)
								at_len<-(at_len+local_len)
							}
						}
						if(is_sample){
							if(max_score>=cut_score){
								num_samples<-(num_samples+1)
							}
							if(max_score>max_score_sample){
								max_score_sample<-max_score
							}						
							if(max_num_peaks>num_peaks_sample){
								num_peaks_sample<-max_num_peaks
							}
						}
						if(is_blank){
							if(max_score>=cut_score){
								num_blanks<-(num_blanks+1)
							}
							if(max_score>max_score_blank){
								max_score_blank<-max_score
							}
							if(max_num_peaks>num_peaks_blank){
								num_peaks_blank<-max_num_peaks
							}	
						}						
					}
				}
				ratios<-c()
				wei<-c()
				for(j in 1:length(centro_sample)){
					if( (length(centro_sample[[j]])>0) & (length(centro_blank[[j]])>0) ){
						ratios<-c(ratios,( mean(centro_sample[[j]])/mean(centro_blank[[j]]) ) )
						wei<-c(wei,((length(centro_sample[[j]])>0)+(length(centro_blank[[j]])>0)))
					}
				}
				if(length(ratios)>0){			
					mean_int_ratio[[i]]<-mean(x=ratios,w=wei)
				}
				num_samples_all[i]<-num_samples
				num_blanks_all[i]<-num_blanks
				max_score_sample_all[i]<-max_score_sample
				max_score_blank_all[i]<-max_score_blank
				num_peaks_sample_all[i]<-num_peaks_sample
				num_peaks_blank_all[i]<-num_peaks_blank	
			}
		}
		##########################################################################################
		# Table with adducts per compound itemized ###############################################
		results_table_1<-data.frame(
			IDed,named,adducted,
			num_samples_all,
			round(max_score_sample_all,digits=2),
			num_peaks_sample_all,
			num_blanks_all,
			round(max_score_blank_all,digits=2),
			num_peaks_blank_all,
			round(mean_int_ratio,digits=1),
			stringsAsFactors=FALSE
		)
		names(results_table_1)<-c(
			"ID","compound","adduct",
			"Sample matches",
			"Max. sample score",
			"Max. sample peaks",
			"Blank matches",
			"Max. blank score",
			"Max. blank peaks",
			"Int. ratio sample/blank"
		)
		##########################################################################################
		# Table with adducts per compound summarized #############################################		
		ID_comp<-unique(IDed)
		adduct_sum<-rep("",length(ID_comp))
		named_sum<-rep("",length(ID_comp))
		max_score_sample_all_sum<-rep(0,length(ID_comp))
		max_score_blank_all_sum<-rep(0,length(ID_comp))		
		num_peaks_sample_all_sum<-rep(0,length(ID_comp))
		num_peaks_blank_all_sum<-rep(0,length(ID_comp))		
		for(i in 1:length(ID_comp)){
			those<-(IDed==ID_comp[i])
			those[!((max_score_sample_all[those]>0) | (max_score_blank_all[those]>0))]<-FALSE
			if(any(those)){	
				named_sum[i]<-unique(named[those])
				adduct_sum[i]<-paste(adducted[those],collapse=", ")
				max_score_sample_all_sum[i]<-max(round(max_score_sample_all[those],digits=2))
				max_score_blank_all_sum[i]<-max(round(max_score_blank_all[those],digits=2))			
				num_peaks_sample_all_sum[i]<-max(num_peaks_sample_all[those])
				num_peaks_blank_all_sum[i]<-max(num_peaks_blank_all[those])
			}else{
				those<-(IDed==ID_comp[i])
				named_sum[i]<-unique(named[those])
			}	
		}
		results_table_2<-data.frame(
			ID_comp,named_sum,adduct_sum,
			max_score_sample_all_sum,
			num_peaks_sample_all_sum,
			max_score_blank_all_sum,
			num_peaks_blank_all_sum
		)
		names(results_table_2)<-c(		
			"ID","compound","adducts",
			"Max. sample score",
			"Max. sample peaks",
			"Max. blank score",
			"Max. blank peaks"
		)
		##########################################################################################
		results<-list()
		results[[1]]<-results_table_1
		results[[2]]<-results_table_2
		if(at_len>0){
			at_matrix<-at_matrix[1:(at_len-1),]
			results[[3]]<-at_matrix
		}else{
			results_table_3<-numeric(0)
		}	
		return(results)
	}

	
	
	
	
	
	
	
	
	
	
	
	