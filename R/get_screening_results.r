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

	get_screening_results<-function(screened_listed,pattern,profileList,measurements_table,compound_table,cut_score){

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
								if(!is.na(screened_listed[[i]][[m]][[k]]$score_2)){
									local_score<-(local_score+screened_listed[[i]][[m]][[k]]$score_2)
								}										
								if(local_score>max_score){
									max_score<-local_score
								}
								if(length(screened_listed[[i]][[m]][[k]]$Peaks[,1])>max_num_peaks){
									max_num_peaks<-length(screened_listed[[i]][[m]][[k]]$Peaks[,1])
								}
								if(is_sample & (local_score>cut_score)){
									for(d in 1:length(screened_listed[[i]][[m]][[k]][[1]][,1])){
										centro_sample[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]]<-c(
											centro_sample[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]],
											profileList[[2]][screened_listed[[i]][[m]][[k]][[1]][d,2],2]	
										)
									}
									next; # no need to check blank
								}
								if(is_blank & (local_score>cut_score)){
									for(d in 1:length(screened_listed[[i]][[m]][[k]][[1]][,1])){
										centro_blank[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]]<-c(
											centro_blank[[ screened_listed[[i]][[m]][[k]][[1]][d,1] ]],
											profileList[[2]][screened_listed[[i]][[m]][[k]][[1]][d,2],2]	
										)
									}
								}									
							}
						}
						if(is_sample){
							if(max_score>cut_score){
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
							if(max_score>cut_score){
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
				num_samples_all[[i]]<-num_samples
				num_blanks_all[[i]]<-num_blanks
				max_score_sample_all[[i]]<-max_score_sample
				max_score_blank_all[[i]]<-max_score_blank
				num_peaks_sample_all[[i]]<-num_peaks_sample
				num_peaks_blank_all[[i]]<-num_peaks_blank	
			}
		}
		results_table_1<-data.frame(
			IDed,
			named,
			adducted,
			num_samples_all,
			round(max_score_sample_all,digits=2),
			num_peaks_sample_all,
			num_blanks_all,
			round(max_score_blank_all,digits=2),
			num_peaks_blank_all,
			mean_int_ratio,
			stringsAsFactors=FALSE
		)
		names(results_table_1)<-c(
			"ID",
			"compound",
			"adduct",
			"Sample matches",
			"Max. sample score",
			"Max. sample peaks",
			"Blank matches",
			"Max. blank score",
			"Max. blank peaks",
			"Mean ratio sample/blank"
		)
		
		
		
		
		results<-list()
		results[[1]]<-results_table_1
		return(results)
	}
