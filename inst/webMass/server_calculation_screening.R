
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
RT_tol_inside<-15

		
#pattern_pos_IS<-c(pattern_pos_IS,pattern_pos_IS,pattern_pos_IS)
#length(pattern_pos_IS)
#patternRT_pos_IS<-rep(patternRT_pos_IS,3)
#patternDelRT_pos_IS<-rep(patternDelRT_pos_IS,3)
		
system.time({		

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
		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos)
		IS_pos_screen_listed<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				IS_pos_screen_listed[[i]]<-list() # m-level		
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
								if(length(IS_pos_screen_listed[[i]])<profileList_pos[[2]][m,6][[1]] ){
									IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-matrix(ncol=2,nrow=0)	# sample level
								}
								IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]]<-rbind(
									IS_pos_screen_listed[[i]][[profileList_pos[[2]][m,6][[1]]]],c(j,m)
								)
							}							
						}
					}
				}
			}
		}
		
		
		
		
})
		
system.time({	

		recomb_score<-function(cent_peak_mat,pattern,profileList,LOD,RT_tol_inside){
			
			#######################################################################
			if(!is.matrix(cent_peak_mat)){stop("cent_peak_mat must be a matrix")}
			if(!is.matrix(pattern)){stop("pattern must be a matrix")}
			if((!is.numeric(LOD))||(length(LOD)>1)){stop("LOD must be numeric")}
			#######################################################################		
			# internal function to check plausibility #############################
			check_plaus<-function(cent_peak_mat,profileList,RT_tol_inside){
				# if only one row with (centroid,peak) left = all is plausible
				# - plausibility of this one (centroid,peak) must be given at this stage, can be inherited further!
				# unique centroids combined?
				if(any(duplicated(cent_peak_mat[,1]))){return(FALSE)}
				# unique peaks combined?
				if(any(duplicated(cent_peak_mat[,2]))){return(FALSE)}
				# all within small RT window ?
				rangeRT<-range(profileList[[2]][cent_peak_mat[,2],3])
				if((rangeRT[2]-rangeRT[1])>RT_tol_inside){return(FALSE)}
				# does intensity pattern match?
				if(length(cent_peak_mat[,1])>1){ # ... which needs more than one peak
	# FINISH			
	# use intercept = 0?
				}
				# else all is plausible
				return(TRUE)
			}
			#######################################################################		
			results<-list(0)
			at_results<-1
			checked<-TRUE
			check_nodes<-list()
			check_nodes[[1]]<-cent_peak_mat # initialize with full set
			check_nodes_index<-list()
			check_nodes_index[[1]]<-length(cent_peak_mat[,1])
			while(checked){
				new_nodes<-list()
				new_nodes_index<-list()
				at_new_nodes<-1
				checked<-FALSE
				for(k in 1:length(check_nodes)){
					if( check_plaus(check_nodes[[k]],profileList,RT_tol_inside) ){
						results[[at_results]]<-check_nodes[[k]]
						at_results<-(at_results+1)
					}else{
						# maker smaller combinations by omission of one (centroid,peak)
						if(check_nodes_index[[k]]>0){ 
							# nothing to inherit - 
							# - this combination is either part of another larger one or
							# - has been build from low-combining that one
							len<-(length(check_nodes[[k]][,1]):1)
							len<-(len[1:check_nodes_index[[k]]])
							len<-rev(len)
							for(z in 1:check_nodes_index[[k]]){
								new_nodes[[at_new_nodes]]<-check_nodes[[k]][-(len[z]),,drop=FALSE]
								new_nodes_index[[at_new_nodes]]<-(check_nodes_index[[k]]-z)
								at_new_nodes<-(at_new_nodes+1)
							}
							checked<-TRUE
						}
					}
				}
				check_nodes<-new_nodes
				check_nodes_index<-new_nodes_index
			}
			return(results)
			#######################################################################
		
		}
		
		many<-0
		many_unamb<-0
		res_IS_pos_screen<-list()  # default: no match at all
		j<-200
		for(i in j:length(IS_pos_screen_listed)){ # i - on compound_adduct
			if(length(IS_pos_screen_listed[[i]])>0){	
				for(m in 1:length(IS_pos_screen_listed[[i]])){ # m - sample
					if(length(IS_pos_screen_listed[[i]][[m]])>0){
						combination_matches<-recomb_score(
							cent_peak_mat=IS_pos_screen_listed[[i]][[m]],
							pattern=pattern_pos_IS[[i]],
							profileList=profileList_pos,
							LOD=cutint,
							RT_tol_inside=RT_tol_inside
						)
						res_IS_pos_screen[[i]]<-combination_matches
						if(length(res_IS_pos_screen[[i]])>0){
							plot(
								log10(pattern_pos_IS[[i]][res_IS_pos_screen[[i]][[1]][,1],2]),
								log10(profileList_pos[[2]][res_IS_pos_screen[[i]][[1]][,2],2]),
								pch=19,cex=.8,main=paste(i,"-",m))
								Sys.sleep(.01)
						}
						if(length(combination_matches)>1){many_unamb<-(many_unamb+1)}
						many<-(many+1)
					}
				}
			}
		}
	


	
})		
		
		

k<-433
these_peaks<-res_IS_pos_screen[[k]][[1]][,2]
pattern_pos_IS[[k]]
res_IS_pos_screen[[k]][[1]][,1]
cbind(res_IS_pos_screen[[k]][[1]][,1],profileList_pos[[2]][these_peaks,c(1,2)])
		
		
		
		
system.time({	



		# resort to a full result list: pattern x sample x (centroids,matches) ( = peak index in profileList_pos)
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



