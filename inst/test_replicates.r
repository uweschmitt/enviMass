########################################################################
# On blind subtraction #################################################
cat("\n Checking replicate intersection ...")
if(
	(logfile$workflow[names(logfile$workflow)=="replicates"])=="yes"
){
	cat(" included ...")
	####################################################################
	ppm<-logfile$parameters$replicate_ppm
	mz_tol<-as.numeric(logfile$parameters$replicate_dmz)
	rt_tol<-as.numeric(logfile$parameters$replicate_delRT)
	int_tol<-10^(as.numeric(logfile$parameters$replicate_IS_dInt))
	####################################################################
	measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
	measurements<-measurements[measurements[,"include"]=="TRUE",]
	if(any(measurements$tag3!="FALSE")){
		cat(" files available ...")
		for(i in 1:length(measurements[,"ID"])){
			cat("\n - ")
			if( (measurements$Type[i]=="sample") & (measurements$tag3[i]!="FALSE") ){
				load(file=file.path(logfile[[1]],"peaklist",as.character(measurements$ID[i])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
				peaks_sam<-peaklist;rm(peaklist)
				# on peaks NOT present in all replicates ###############
				if(any(peaks_sam[,colnames(peaks_sam)=="keep"]==0)){
					those<-which(peaks_sam[,colnames(peaks_sam)=="keep"]==0)
					if(length(those)>sam){
						those<-sample(those,sam,replace=FALSE)
					}
					found<-vector(mode="list", length=length(those))
					for(k in 1:length(those)){
						MASS<-peaks_sam[those[k],colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_sam[those[k],colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_sam[those[k],colnames(peaks_sam)=="sum_int_corr"]
						if(ppm){
							mz_low<-(MASS-(mz_tol*MASS/1E6))
							mz_up<-(MASS+(mz_tol*MASS/1E6))
						}else{
							mz_low<-(MASS-(mz_tol/1000))
							mz_up<-(MASS+(mz_tol/1000))							
						}
						RT_low<-(RT-rt_tol)
						RT_up<-(RT+rt_tol)
						int_low<-(INT-int_tol)
						int_up<-(INT+int_tol)						
						found[[k]]<-matrix(nrow=1,ncol=7,
							c(mz_low,mz_up,int_low,int_up,RT_low,RT_up,1)
						)
					}
					IDs<-measurements[measurements$tag3==measurements$tag3[i],"ID"]
					IDs<-IDs[IDs!=measurements$ID[i]]
					for(k in 1:length(IDs)){
						cat("*");
						load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[k])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
						peaks_repl<-peaklist;rm(peaklist)
						MASS<-peaks_repl[,colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_repl[,colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_repl[,colnames(peaks_sam)=="sum_int_corr"]						
						if(ppm){
							mz_low<-(MASS-(mz_tol*MASS/1E6))
							mz_up<-(MASS+(mz_tol*MASS/1E6))
						}else{
							mz_low<-(MASS-(mz_tol/1000))
							mz_up<-(MASS+(mz_tol/1000))							
						}
						RT_low<-(RT-rt_tol)
						RT_up<-(RT+rt_tol)
						int_low<-(INT-int_tol)
						int_up<-(INT+int_tol)			
						for(l in 1:length(found)){
							new_matrix<-matrix(ncol=7,nrow=1,0)
							for(m in 1:length(found[[l]][,1])){
								got<-which(
									(mz_low<=found[[l]][m,2]) &
									(mz_up>=found[[l]][m,1]) &
									(int_low<=found[[l]][m,4]) &
									(int_up>=found[[l]][m,3]) &
									(RT_low<=found[[l]][m,6]) &
									(RT_up>=found[[l]][m,5]) 
								)						
								if(length(got)>0){
									for(n in 1:length(got)){
										new_matrix<-rbind(
											new_matrix,
											matrix(
												ncol=7,nrow=1,
												c(
													max(c(mz_low[got[n]],found[[l]][m,1])),	
													min(c(mz_up[got[n]],found[[l]][m,2])),
													max(c(int_low[got[n]],found[[l]][m,3])),	
													min(c(int_up[got[n]],found[[l]][m,4])),											
													max(c(RT_low[got[n]],found[[l]][m,5])),	
													min(c(RT_up[got[n]],found[[l]][m,6])),
													(found[[l]][m,7]+1)
												)
											)
										)
									}
								}else{
									new_matrix<-rbind(
										new_matrix,
										found[[l]][m,,drop=FALSE]
									)
								}
							}	
							found[[l]]<-new_matrix[-1,,drop=FALSE]
						}
					}				
					for(k in 1:length(found)){
						if(any(found[[k]][,7]>=(length(IDs)+1))){
							cat("WARNING: missing replicate peaks detected!")				
						}
					}
				}
				# on peaks NOT present in all replicates ###############		
				if(any(peaks_sam[,colnames(peaks_sam)=="keep"]==1)){
					those<-which(peaks_sam[,colnames(peaks_sam)=="keep"]==1)
					if(length(those)>sam){
						those<-sample(those,sam,replace=FALSE)
					}
					found<-vector(mode="list", length=length(those))
					for(k in 1:length(those)){
						MASS<-peaks_sam[those[k],colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_sam[those[k],colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_sam[those[k],colnames(peaks_sam)=="sum_int_corr"]
						if(ppm){
							mz_low<-(MASS-(mz_tol*MASS/1E6))
							mz_up<-(MASS+(mz_tol*MASS/1E6))
						}else{
							mz_low<-(MASS-(mz_tol/1000))
							mz_up<-(MASS+(mz_tol/1000))							
						}
						RT_low<-(RT-rt_tol)
						RT_up<-(RT+rt_tol)
						int_low<-(INT-int_tol)
						int_up<-(INT+int_tol)						
						found[[k]]<-matrix(nrow=1,ncol=7,
							c(mz_low,mz_up,int_low,int_up,RT_low,RT_up,1)
						)
					}
					IDs<-measurements[measurements$tag3==measurements$tag3[i],"ID"]
					IDs<-IDs[IDs!=measurements$ID[i]]
					for(k in 1:length(IDs)){
						cat(".");
						load(file=file.path(logfile[[1]],"peaklist",as.character(IDs[k])),envir=as.environment(".GlobalEnv"),verbose=FALSE);
						peaks_repl<-peaklist;rm(peaklist)
						MASS<-peaks_repl[,colnames(peaks_sam)=="m/z_corr"]
						RT<-peaks_repl[,colnames(peaks_sam)=="RT_corr"]
						INT<-peaks_repl[,colnames(peaks_sam)=="sum_int_corr"]						
						if(ppm){
							mz_low<-(MASS-(mz_tol*MASS/1E6))
							mz_up<-(MASS+(mz_tol*MASS/1E6))
						}else{
							mz_low<-(MASS-(mz_tol/1000))
							mz_up<-(MASS+(mz_tol/1000))							
						}
						RT_low<-(RT-rt_tol)
						RT_up<-(RT+rt_tol)
						int_low<-(INT-int_tol)
						int_up<-(INT+int_tol)			
						for(l in 1:length(found)){
							new_matrix<-matrix(ncol=7,nrow=1,0)
							for(m in 1:length(found[[l]][,1])){
								got<-which(
									(mz_low<=found[[l]][m,2]) &
									(mz_up>=found[[l]][m,1]) &
									(int_low<=found[[l]][m,4]) &
									(int_up>=found[[l]][m,3]) &
									(RT_low<=found[[l]][m,6]) &
									(RT_up>=found[[l]][m,5]) 
								)						
								if(length(got)>0){
									for(n in 1:length(got)){
										new_matrix<-rbind(
											new_matrix,
											matrix(
												ncol=7,nrow=1,
												c(
													max(c(mz_low[got[n]],found[[l]][m,1])),	
													min(c(mz_up[got[n]],found[[l]][m,2])),
													max(c(int_low[got[n]],found[[l]][m,3])),	
													min(c(int_up[got[n]],found[[l]][m,4])),											
													max(c(RT_low[got[n]],found[[l]][m,5])),	
													min(c(RT_up[got[n]],found[[l]][m,6])),
													(found[[l]][m,7]+1)
												)
											)
										)
									}
								}else{
									new_matrix<-rbind(
										new_matrix,
										found[[l]][m,,drop=FALSE]
									)
								}
							}	
							found[[l]]<-new_matrix[-1,,drop=FALSE]
						}
					}				
					for(k in 1:length(found)){
						if(all(found[[k]][,7]<(length(IDs)+1))){
							cat("WARNING: non-replicated replicate peaks detected!")				
						}
					}
				}
			}
		}
	}
	####################################################################
}
cat("\n Checking replicate intersection ... done.")
########################################################################


