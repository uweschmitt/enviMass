

measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
if(  
	any(duplicated(measurements$tag3[measurements$tag3!="FALSE"]))
){

	# DEVEL - TO BE DELETED !!! ####################################
	if(TRUE){
		measurements$tag3[c(2,3)]<-"group_1"
		measurements$tag3[c(9,11,12)]<-"group_2"	
		ppm<-TRUE
		mz_tol<-3
		rt_tol<-30
		with_test<-TRUE
	}
	
	replic<-(measurements$tag3[measurements$tag3!="FALSE"])
	replic<-replic[duplicated(replic)]
	replic<-unique(replic)
	for(i in 1:length(replic)){
	
			cat(paste("\n    Replicate intersection",replic[i],":"));
			IDs<-measurements$ID[measurements$tag3==replic[i]]
			if(any(duplicated(IDs))){stop("replicates: non-unique IDs found!")} # should not happen anyway
			# initialize intersection rectangles with first peaklist ################################
			load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[1])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
			peaklist<-peaklist[,c(12,13,14)]
			if(ppm){
				low_mass<-(peaklist[,1]-(mz_tol/1E6*peaklist[,1]))
				high_mass<-(peaklist[,1]+(mz_tol/1E6*peaklist[,1]))
			}else{
				low_mass<-(peaklist[,1]-mz_tol)
				high_mass<-(peaklist[,1]+mz_tol)	
			}
			low_rt<-(peaklist[,3]-rt_tol)
			high_rt<-(peaklist[,3]+rt_tol)			
			ID<-seq(1,length(peaklist[,1]),1)
			intersection<-cbind(low_mass,high_mass,low_rt,high_rt,ID)
			# test 
			if(with_test){
				intersection<-rbind(
					intersection,
					c(20,20,5,15,-2),
					c(.3,.3001,39,41,-3),
					c(800.001,800.002,1400,1430,-4))
			}
			rm(peaklist)
			for(j in 2:length(IDs)){
				# load next peaklist
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				peaklist<-peaklist[,c(12,13,14)]
				if(ppm){
					low_mass<-(peaklist[,1]-(mz_tol/1E6*peaklist[,1]))
					high_mass<-(peaklist[,1]+(mz_tol/1E6*peaklist[,1]))
				}else{
					low_mass<-(peaklist[,1]-mz_tol)
					high_mass<-(peaklist[,1]+mz_tol)	
				}
				low_rt<-(peaklist[,3]-rt_tol)
				high_rt<-(peaklist[,3]+rt_tol)		
				query<-cbind(low_mass,high_mass,low_rt,high_rt)			
				if(with_test){
					query<-rbind(
						query,
						c(20,20,5,15),
						c(.3,.3001,39,41),
						c(800.001,800.002,1400,1430))
				}
				# build a boxtree with next peaklist
				tree <- .Call("boxtree", 
					as.matrix(query),
					PACKAGE="nontarget"
				)	
				colnames(tree)<-c("LOSON","MIDSON","HISON","level","disc")		
				# query
				search_bounds<-rep(0,4)
				new_intersection<-matrix(nrow=10000,ncol=(4+j),-1)
				len<-10000
				at_new<-1
				for(k in 1:length(intersection[,1])){
					search_bounds[1]<-intersection[k,1]
					search_bounds[2]<-intersection[k,2]				
					search_bounds[3]<-intersection[k,3]
					search_bounds[4]<-intersection[k,4]					
					found <- .Call("search_boxtree", 
						query,
						tree,
						as.numeric(search_bounds),
						as.integer(1), # return full findings
						PACKAGE="nontarget"
					)
					if(length(found)>0){
						for(m in 1:length(found)){
							new_intersection[at_new,1]<-max(c(intersection[k,1],query[found[m],1]))
							new_intersection[at_new,2]<-min(c(intersection[k,2],query[found[m],2]))							
							new_intersection[at_new,3]<-max(c(intersection[k,3],query[found[m],3]))							
							new_intersection[at_new,4]<-min(c(intersection[k,4],query[found[m],4]))
							new_intersection[at_new,(5:(5+j-2))]<-intersection[k,(5:(5+j-2))]
							new_intersection[at_new,(4+j)]<-found[m]
							at_new<-(at_new+1)
							if(at_new>len){ # extend new_intersection
								new_intersection<-rbind(
									new_intersection,
									matrix(nrow=10000,ncol=(4+j),-1)
								)
								len<-(len+10000)
							}
						}	
					}				
				}
				intersection<-new_intersection[new_intersection[,5]!=-1,,drop=FALSE] # omit empty entries
				if(at_new==1){break} # no further intersections detected - abort
			}
			if(with_test){
				if(!any(intersection[,5]==-2)){stop("intersection test_1 failed")}
				if(!any(intersection[,5]==-3)){stop("intersection test_1 failed")}
				if(!any(intersection[,5]==-4)){stop("intersection test_1 failed")}	
				intersection<-intersection[intersection[,5]>0,]				
			}
			# clean peaklists
			for(j in 1:length(IDs)){
				load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				keep<-rep(0,length(peaklist[,1]))
				if(length(intersection[,1])>0){ # any intersections at all?
					keep_those<-unique(intersection[,(4+j)])
					keep[keep_those]<-1
				}
				cat(paste("\n Keep",sum(keep==1),"of",length(keep),"peaks"))
				peaklist[,colnames(peaklist)=="keep"]<-keep
				save(peaklist,file=file.path(logfile[[1]],"peaklist",as.character(IDs[j])))
				rm(peaklist)
			}
			
			
	}
	logfile$summary[11,2]<<-"TRUE";
    logfile$summary[11,2]<-"TRUE";
	#logfile$Tasks_to_redo[3]<-"FALSE";
	#logfile$Tasks_to_redo[3]<<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[11,2]<-"done"
	summa[11,2]<<-"done"
	output$summa_html<<-renderText(summary_html(summa));
    cat("Replicate intersection done \n");
    output$dowhat<<-renderText("Replicate intersection done ... wait")
			
}else{

	logfile$summary[11,2]<<-"FALSE";
    logfile$summary[11,2]<-"FALSE";
    save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
    summa[11,2]<-"skipped"
	summa[11,2]<<-"skipped"
	output$summa_html<<-renderText(summary_html(summa));
    cat("Replicate intersection skipped \n");
    output$dowhat<<-renderText("Replicate intersection skipped ... wait")

}




