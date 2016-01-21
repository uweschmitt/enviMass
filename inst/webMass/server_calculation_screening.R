
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
		# resort to a full result list: pattern x sample x centroids x matches ( = peak index in profileList_pos)
		res_IS_pos_screen<-list()  # default: no match at all
		for(i in 1:length(screen_list)){ # over compound x adduct = i
			if(any(is.na(screen_list[[i]]==FALSE))){
				res_IS_pos_screen[[i]]<-list() # m-level		
				for(j in 1:length(screen_list[[i]])){ # over its centroids = j
					if(screen_list[[i]][[j]]!="FALSE"){ 
						profs<-as.numeric(strsplit(screen_list[[i]][[j]]," / ")[[1]])
						for(k in 1:length(profs)){ # over their matched profiles = k
							if(profileList_pos[[7]][profs[k],4]!=profs[k]){cat("\n debug me: profile ID mismatch");stop();} # just a check
							for(m in profileList_pos[[7]][profs[k],1]:profileList_pos[[7]][profs[k],2]){ # over their samples
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
		# calculate combinations over centroids&peaks per sample per compoundadduct
		many<-0
		for(i in 1:length(res_IS_pos_screen)){
			if(!is.null(res_IS_pos_screen[[i]])){
				for(j in 1:length(res_IS_pos_screen[[i]])){
					if(!is.null(res_IS_pos_screen[[i]][[j]])){
						many<-(many+1)
			
			
			
					}
				}
			}
		}
				
})		
		

		
	
	
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



