
if(  
	#(logfile$workflow[2]=="yes" && logfile$summary[5,2]=="FALSE")  || 
	#(logfile$Tasks_to_redo[3]=="TRUE") 
	FALSE
){

    ############################################################################
	# for IS compounds #########################################################
	# positive ionization ######################################################
	if(TRUE){
	
		if(any(objects(envir=as.environment(".GlobalEnv"))=="peaklist")){rm(peaklist,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="peaklist")){rm(peaklist)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="profileList_pos")){rm(profileList_pos,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="profileList_pos")){rm(profileList_pos)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="pattern_pos_IS")){rm(pattern_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="pattern_pos_IS")){rm(pattern_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternRT_pos_IS")){rm(patternRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternRT_pos_IS")){rm(patternRT_pos_IS)}
		if(any(objects(envir=as.environment(".GlobalEnv"))=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS,envir=as.environment(".GlobalEnv"))}
		if(any(objects()=="patternDelRT_pos_IS")){rm(patternDelRT_pos_IS)}
		load(file=file.path(as.character(logfile[[1]]),"results","profileList_pos"),envir=as.environment(".GlobalEnv"));	
		load(file=file.path(logfile[[1]],"results","pattern_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		load(file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"),envir=as.environment(".GlobalEnv"));
		
		peaks<-profileList_pos[[7]];
		peaklist<-peaks[,c(14,16,15)];
		
		
ppm<-TRUE
mztol<-4
cutint<-2E3
int_tol<-.25
RT_tol_inside<-0.3

		
#pattern_pos_IS<-c(pattern_pos_IS,pattern_pos_IS,pattern_pos_IS)
#length(pattern_pos_IS)
#patternRT_pos_IS<-rep(patternRT_pos_IS,3)
#patternDelRT_pos_IS<-rep(patternDelRT_pos_IS,3)
		
#system.time({		

		# screen centroids
		count_nonmax<-0
		for(i in 1:length(pattern_pos_IS)){
			count_nonmax<-(count_nonmax+
				length(pattern_pos_IS[[i]][,1])
			)
		}
		centro_mass<-rep(0,count_nonmax)
		centro_ID<-rep(0,count_nonmax)
		centro_number<-rep(0,count_nonmax)
		centro_RT<-rep(0,count_nonmax)
		centro_dRT<-rep(0,count_nonmax)
		at_ID<-1
		screen_list<-as.list(rep("FALSE",length(pattern_pos_IS)))
		for(i in 1:length(pattern_pos_IS)){
			n<-length(pattern_pos_IS[[i]][,1])
			centro_mass[at_ID:(at_ID+n-1)]<-pattern_pos_IS[[i]][,1]
			centro_ID[at_ID:(at_ID+n-1)]<-i
			centro_number[at_ID:(at_ID+n-1)]<-(1:n)
			centro_RT[at_ID:(at_ID+n-1)]<-patternRT_pos_IS[i]
			centro_dRT[at_ID:(at_ID+n-1)]<-patternDelRT_pos_IS[i]
			screen_list[[i]]<-as.list(rep("FALSE",n))
			at_ID<-(at_ID+n)
		}
		getit <- search_peak( ### adapt mz tolerances
			peaklist, 
			centro_mass, 
			dmz=5, 
			ppm=TRUE, 
			RT=centro_RT, 
			dRT=centro_dRT)	
		for(i in 1:length(getit)){ # transfer to a fist list of compoundadduct x centroids
			screen_list[[centro_ID[i]]][[centro_number[i]]]<-getit[i]
		}
		# resort to a full result list: pattern x sample x centroids x matches ( = peak index in profileList_pos)
		res_IS_pos_screen<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				res_IS_pos_screen[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profile peaks = k
							if(profileList_pos[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_pos[[7]][profs[k],1]:profileList_pos[[7]][profs[k],2]){ # over their sample peaks
								delmass<-abs(profileList_pos[[2]][m,1]-pattern_pos_IS[[i]][j,1])		
								if(!ppm){
									if(delmass>mztol){next}
								}else{
									if(delmass*1E6/pattern_pos_IS[[i]][j,1]>mztol){next}
								}
								if(length(res_IS_pos_screen[[i]])<profileList_pos[[2]][m,6][[1]] ){
									res_IS_pos_screen[[i]][[profileList_pos[[2]][m,6][[1]] ]]<-list() 		# sample level
								}
								if(length(res_IS_pos_screen[[i]][[profileList_pos[[2]][m,6][[1]]]])<j){
									res_IS_pos_screen[[i]][[ profileList_pos[[2]][m,6][[1]]]][[j]]<-list()	# centroid level
								}									
								len<-length(res_IS_pos_screen[[i]][[profileList_pos[[2]][m,6][[1]]]][[j]]) # peak level
								res_IS_pos_screen[[i]][[profileList_pos[[2]][m,6][[1]]]][[j]][[len+1]]<-m
							}							
						}
					}
				}
			}
		}

		
system.time({	
		
		# calculate combinations over centroids & peaks per sample per compound_adduct
		many<-0
		doubled<-0
		max_score<-rep(0,length(pattern_pos_IS))
		max_peaks<-rep(0,length(pattern_pos_IS))
		res_IS_pos_screen_logical<-res_IS_pos_screen # store logical results - peak matched?
		for(i in 1:length(res_IS_pos_screen)){ # i - on compound_adduct
			if(length(res_IS_pos_screen[[i]])>0){
				for(m in 1:length(res_IS_pos_screen[[i]])){ # m - sample
					if(length(res_IS_pos_screen[[i]][[m]])>0){
						# retrieve ALL unique two-combinations of peaks ##########################################
						combis<-list();
						at_combis<-1;
						if(length(res_IS_pos_screen[[i]][[m]])>1){
						for(j in 1:(length(res_IS_pos_screen[[i]][[m]])-1)){ # j - on centroid
							if(length(res_IS_pos_screen[[i]][[m]][[j]])>0){
								for(k in 1:length(res_IS_pos_screen[[i]][[m]][[j]])){ # k - on peak 
									if(length(res_IS_pos_screen[[i]][[m]][[j]][[k]])>0){
										rescaled_intens_lower<-(pattern_pos_IS[[i]][,2]*(
											( profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]-(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]) )
											/pattern_pos_IS[[i]][j,1])
										)									
										rescaled_intens_upper<-(pattern_pos_IS[[i]][,2]*(
											( profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]+(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]) )
											/pattern_pos_IS[[i]][j,1])
										)									
										for(b in (j+1):length(res_IS_pos_screen[[i]][[m]])){ # b - on centroid
											if(length(res_IS_pos_screen[[i]][[m]][[b]])>0){			
												for(d in 1:length(res_IS_pos_screen[[i]][[m]][[b]])){
													if(res_IS_pos_screen[[i]][[m]][[j]][[k]]!=res_IS_pos_screen[[i]][[m]][[b]][[d]]){ # must be two different peaks!
														if(
															RT_tol_inside<abs(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),3]-profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),3])
														){next} # within small RT window?
														if(
															(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]+(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]))<rescaled_intens_lower[b]
														){next}
														if(
															(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]-(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]))>rescaled_intens_upper[b]
														){next}													
															combis[[at_combis]]<-list()
															combis[[at_combis]][[1]]<-c(j,k)
															combis[[at_combis]][[2]]<-c(b,d)
															at_combis<-(at_combis+1)
															many<-(many+1)
													}else{
														doubled<-(doubled+1)
													}
												}
											}
										}									
									}
								}
							}
						}
						}else{"finish me"} # singleton
						# recombine these 2-tupels ###############################################################
						#if(length(combis)>0){print(length(combis))}
						#if(length(combis)>2){stop()}
						fin_combis<-list()
						if(length(combis)>1){
							do_combis<-TRUE
							while(do_combis){
								new_combis<-list()
								for(b in 1:length(combis)){
									for(d in (b+1):length(combis)){
								


								
									}
								}
								
								if(length(new_combis)<2){do_combis<-FALSE}
							}
						}else{fin_combis<-combis}
						# annotate results #######################################################################
						for(b in 1:length(fin_combis)){
						
						
						
						}
						
						
						
					}
				}
			}
		}
		
				
})		
		
		
system.time({	
		
		# calculate combinations over centroids&peaks per sample per compound_adduct
		many<-0
		max_score<-rep(0,length(pattern_pos_IS))
		max_peaks<-rep(0,length(pattern_pos_IS))
		res_IS_pos_screen_logical<-res_IS_pos_screen # store logical results - peak matched?
		for(i in 1:length(res_IS_pos_screen)){ # i - on compound_adduct
			if(length(res_IS_pos_screen[[i]])>0){
				for(m in 1:length(res_IS_pos_screen[[i]])){ # m - sample
					if(length(res_IS_pos_screen[[i]][[m]])>0){
						for(j in 1:length(res_IS_pos_screen[[i]][[m]])){ # j - on centroid
							if(length(res_IS_pos_screen[[i]][[m]][[j]])>0){
								for(k in 1:length(res_IS_pos_screen[[i]][[m]][[j]])){ # k - on peak 
									if(length(res_IS_pos_screen[[i]][[m]][[j]][[k]])>0){
									#if(res_IS_pos_screen_logical[[i]][[m]][[j]][[k]]!=0){ # peak positively screened before?
										rescaled_intens_lower<-(pattern_pos_IS[[i]][,2]*(
											( profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]-(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]) )
											/pattern_pos_IS[[i]][j,1])
										)									
										rescaled_intens_upper<-(pattern_pos_IS[[i]][,2]*(
											( profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]+(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),2]) )
											/pattern_pos_IS[[i]][j,1])
										)									
										found_intens<-rep(FALSE,length(pattern_pos_IS[[i]][,1]))
										found_intens[k]<-TRUE
										# scale to that peak and screen over all others - requires an inner loop
										for(b in 1:length(res_IS_pos_screen[[i]][[m]])){ # b - on centroid
										if(length(res_IS_pos_screen[[i]][[m]][[b]])>0){			
										if(b!=j){
											for(d in 1:length(res_IS_pos_screen[[i]][[m]][[b]])){
												# check lower intensity bound
												if(
													(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]+(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]))<rescaled_intens_lower[[b]]
												){next}
												if(
													(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]-(int_tol*profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),2]))>rescaled_intens_upper[[b]]
												){next}
												# check adapted RT window
												if(
													RT_tol_inside<abs(profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[b]][[d]]),3]-profileList_pos[[2]][(res_IS_pos_screen[[i]][[m]][[j]][[k]]),3])
												){next}
												# matched!
												found_intens[b]<-TRUE
												many<-(many+1)
											}
										}
										}
										}
										if(sum(found_intens)>max_peaks[i]){max_peaks[i]<-sum(found_intens)}
										must_intens<-(rescaled_intens_lower>=cutint)
										aim_intens<-sum(rescaled_intens_lower[must_intens])
										if(any(must_intens)){
											local_score<-(sum(rescaled_intens_lower[found_intens & must_intens])/aim_intens)
											if(local_score>max_score[i]){max_score[i]<-local_score}
										}else{ # found below threshold!
											if(max_score[i]==0){
												aim_intens<-sum(rescaled_intens_lower[rescaled_intens_lower>=min(rescaled_intens_lower[found_intens])])
												local_score<--(sum(rescaled_intens_lower[found_intens])/aim_intens)
												if(local_score<max_score[i]){
													max_score[i]<-local_score
												}
											}
										}
										#res_IS_pos_screen_logical[[i]][[m]][[b]][[d]]<-0
										}
									#}
								}
							}
						}
					}
				}
			}
		}
		
				
})		
	
	
			
		# define the result table #################################################
		IS_pos_ID<-sample(c(1,2,3),10, replace =TRUE)
		IS_pos_conc<-sample(c(0.5,0.10,300),10, replace =TRUE)
		IS_pos_score<-sample(c(1,0.8,0.5),10, replace =TRUE)
		IS_pos_flag<-rep('<img src="ID_1.png" height="22"></img>',10)	
	
		
		
		
for(i in 1:1000){		
png(file = file.path(logfile$project_folder,"results","screening",paste("ID_",as.character(i),sep=""),fsep = "\\"), 
	width = 300, height = 80, units = "px", pointsize = 12, bg = "transparent")
	par(mar=c(0,0,0,0))
	plot(1:10,pch=19,cex=.3,col="red")
	rect(1, 5, 3, 7, col = "white")
dev.off()
}		
		
		
		IS_screening_pos<- data.frame(
			ID = IS_pos_ID,
			Concentration = IS_pos_conc,
			Score = IS_pos_score,
			flag = IS_pos_flag#c(
				#'<img src="ID_1.png" height="22"></img>',
				#'<img src="project22/results/screening/ID_1.png" height="32"></img>',
				#'<img src="http://upload.wikimedia.org/wikipedia/commons/thumb/f/fa/Flag_of_the_People%27s_Republic_of_China.svg/200px-Flag_of_the_People%27s_Republic_of_China.svg.png" height="52"></img>'
			#)
		)
		save(
			IS_screening_pos,
			file=file.path(logfile$project_folder,"results","screening","IS_screening_pos",fsep = "\\"),
			precheck=FALSE
		)

		
		
		

		
	
	
	}
	
	
    ############################################################################
	# then on targets pos, before switching to negative list - upload takes too long
	
	
}else{
if(FALSE){
	if(logfile$workflow[2]=="no"){
		logfile$summary[5,2]<<-"FALSE";
		logfile$summary[5,2]<-"FALSE";
		path=file.path(logfile[[1]],"pics","recal_none")
			png(filename = path, bg = "white")
			plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
			dev.off()
		    exprrec<-list(src=path)
			output$recal_pic<-renderImage(exprrec, deleteFile = FALSE);		
			output$peakhist_pic<-renderImage(exprrec, deleteFile = FALSE);
			output$peakmzRT_pic<-renderImage(exprrec, deleteFile = FALSE);	
	}
	logfile$Tasks_to_redo[3]<-"FALSE";
	logfile$Tasks_to_redo[3]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[5,2]<-"skipped"
	summa[5,2]<<-"skipped"
	output$summar<<-renderTable(summa);
    cat("Mass recalibration skipped \n");
    output$dowhat<<-renderText("Recalibration skipped ... wait")
}
}



