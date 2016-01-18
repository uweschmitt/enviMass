if(

	( logfile$summary[4,2]=="FALSE" &  logfile$workflow[2]=="yes") || # run: recal
	( logfile$summary[4,2]=="FALSE" &  logfile$workflow[4]=="yes") || # run: norm
	( logfile$summary[4,2]=="FALSE" &  logfile$workflow[3]=="yes") || # run: allign
	( logfile$Tasks_to_redo[8]=="TRUE") 

){
            measurements<-read.csv(file=file.path(logfile[[1]],"dataframes","measurements"),colClasses = "character");
            ######################################################################
			# on IS ##############################################################
 			intstand<-read.table(file=file.path(logfile[[1]],"dataframes","IS.txt"),header=TRUE,sep="\t",colClasses = "character");
			if(length(intstand$Formula[intstand$Formula!="-"])>0){ # any IS formulas available?

				##################################################################
				if(any(measurements[,4]=="positive") & any(intstand[,7]=="positive") ){
					pattern_pos_IS<-list(0);
					counter<-c(1);
					finform<-c();
					fineidis<-c();
					finadd<-c();
					mz_pos_IS<-c();
					RT_pos_IS<-c();
					patternRT_pos_IS<-c();
					patternDelRT_pos_IS<-c();
					formeln<-intstand[,3];
					eidis<-intstand[,1];
					rets<-intstand[,4];
					tolrets<-intstand[,5];
					rec<-intstand[,8];
					tag1<-intstand[,12];
					tag2<-intstand[,13];
					tag3<-intstand[,14];
					for(i in 1:length(formeln)){
						if(intstand[i,7]=="positive"){
							cat(paste("\nDoing IS positive ",as.character(i)," of ",as.character(length(formeln)),":",sep=""));
							if( (intstand[i,10]=="TRUE") & (intstand[i,6]=="FALSE") ){
								warning(paste("IS pattern: restrict to adduct with no adduct selected? skipped IS with ID:",eidis[i],"\n",sep=""))
								next;
							}
							takeall<-TRUE # use all adducts for recal?
							if(intstand[i,10]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
								if(intstand[i,6]=="FALSE"){ # no main adduct chosen?
									with_add<-c(logfile[[7]])
									takeall<-TRUE
								}else{ # with main adduct?
									with_add<-c(intstand[i,6],logfile[[7]])
									with_add<-unique(with_add)
									takeall<-FALSE
								}
							}else{ # use main adduct only
								with_add<-c(intstand[i,6])
								takeall<-TRUE # there is just one entry anyway
							}
							for(j in 1:length(with_add)){
								formelone<-formeln[i];
								formelone<-enviPat::multiform(formelone,as.numeric(adducts[with_add[j]==adducts[,1],9]))
								if(as.character(adducts[with_add[j]==adducts[,1],7])!="FALSE"){
									formelone<-enviPat::mergeform(formelone,as.character(adducts[with_add[j]==adducts[,1],7]))
								}
								if(as.character(adducts[with_add[j]==adducts[,1],8])!="FALSE"){
									if(
										(enviPat::check_ded(formelone, as.character(adducts[with_add[j]==adducts[,1],8])))==FALSE
									){
										formelone<-enviPat::subform(formelone,as.character(adducts[with_add[j]==adducts[,1],8]))
									}else{
										next; # just skip that deduct!
									}
								}
								finform<-c(finform,formelone);
								fineidis<-c(fineidis,paste(eidis[i],"_",with_add[j],"_",tag1[i],"_",tag2[i],"_",tag3[i],sep=""));
								pattern<-enviPat::isopattern(
									isotopes,
									formelone,
									threshold=0.5,
									plotit=FALSE,
									charge=as.numeric(strsplit(as.character(adducts[with_add[j]==adducts[,1],3]),"+",fixed=TRUE)[[1]][1]),
									emass=0.00054858,
									algo=1
								)
								checked<-enviPat::check_chemform(isotopes, formelone)
								res<-try(enviPat::getR(checked, resolution_list[names(resolution_list) == logfile$parameters[22]][[1]], nknots = 7, spar = 0.1, plotit = FALSE), silent=TRUE)
								if(grepl("Error",res)){
									cat("\n Mass out of range of Resolution data - set to range!");
									if(checked[[3]]<=min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])){
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]
									}else{
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==max(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]									
									}
								}
								pattern<-enviPat::envelope(
									pattern,
									ppm=FALSE,
									dmz=0.0001,
									frac=1/4,
									env="Gaussian",
									resolution=res,
									plotit=FALSE
								)
								pattern<-enviPat::vdetect(
									pattern,
									detect="centroid",
									plotit=FALSE
								)
								pattern_pos_IS[[counter]]<-pattern[[1]]; 								
								patternRT_pos_IS<-c(patternRT_pos_IS,as.numeric(as.character(rets[i]))*60);
								if(tolrets[i]!="FALSE"){
									patternDelRT_pos_IS<-c(patternDelRT_pos_IS,as.numeric(as.character(tolrets[i]))*60)
								}else{
									patternDelRT_pos_IS<-c(patternDelRT_pos_IS,isolate(input$screen_IS_delRT))
								}								
								if(rec[i]=="TRUE"){
									if(takeall){
										mz_pos_IS<-c(mz_pos_IS,pattern_pos_IS[[counter]][pattern_pos_IS[[counter]][,2]==max(pattern_pos_IS[[counter]][,2]),1][counter=1]);
										RT_pos_IS<-c(RT_pos_IS,
											as.numeric(as.character(rets[i]))*60
										);
									}else{
										if(j==1){
											mz_pos_IS<-c(mz_pos_IS,pattern_pos_IS[[counter]][pattern_pos_IS[[counter]][,2]==max(pattern_pos_IS[[counter]][,2]),1][counter=1]);
											RT_pos_IS<-c(RT_pos_IS,
												as.numeric(as.character(rets[i]))*60
											);
										}
									}
								}
								counter<-c(counter+1);
							}
						}
					}
					if(length(mz_pos_IS)>0){
						intmass_pos_IS<-data.frame(mz_pos_IS,RT_pos_IS)
						names(intmass_pos_IS)<-c("m/z","RT")
						save(intmass_pos_IS,file=file.path(logfile[[1]],"results","intmass_pos_IS"))                    
					}else{
						cat("\n no masses z=pos for recalibration formed!")
					}
					names(pattern_pos_IS)<-fineidis    
					save(pattern_pos_IS,file=file.path(logfile[[1]],"results","pattern_pos_IS"))
					save(patternRT_pos_IS,file=file.path(logfile[[1]],"results","patternRT_pos_IS"))
					save(patternDelRT_pos_IS,file=file.path(logfile[[1]],"results","patternDelRT_pos_IS"))
				}
				###################################################################
				
				##################################################################
				if(any(measurements[,4]=="negative") & any(intstand[,7]=="negative")){
					pattern_neg_IS<-list(0);
					counter<-c(1);
					finform<-c();
					fineidis<-c();
					finadd<-c();
					mz_neg_IS<-c();
					RT_neg_IS<-c();
					patternRT_neg_IS<-c();
					patternDelRT_neg_IS<-c();
					formeln<-intstand[,3];
					eidis<-intstand[,1];
					rets<-intstand[,4];
					tolrets<-intstand[,5];
					rec<-intstand[,8];
					tag1<-intstand[,12];
					tag2<-intstand[,13];
					tag3<-intstand[,14];
					for(i in 1:length(formeln)){
						if(intstand[i,7]=="negative"){
							cat(paste("\nDoing IS negative ",as.character(i)," of ",as.character(length(formeln)),":",sep=""));
							if( (intstand[i,10]=="TRUE") & (intstand[i,6]=="FALSE") ){
								warning(paste("IS pattern: restrict to adduct with no adduct selected? skipped IS with ID:",eidis[i],"\n",sep=""))
								next;
							}
							takeall<-TRUE # use all adducts for recal?
							if(intstand[i,10]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
								if(intstand[i,6]=="FALSE"){ # no main adduct chosen?
									with_add<-c(logfile[[8]])
									takeall<-TRUE
								}else{ # with main adduct?
									with_add<-c(intstand[i,6],logfile[[8]])
									with_add<-unique(with_add)
									takeall<-FALSE
								}
							}else{ # use main adduct only
								with_add<-c(intstand[i,6])
								takeall<-TRUE # there is just one entry anyway
							}
							for(j in 1:length(with_add)){
								formelone<-formeln[i];
								formelone<-enviPat::multiform(formelone,as.numeric(adducts[with_add[j]==adducts[,1],9]))
								if(as.character(adducts[with_add[j]==adducts[,1],7])!="FALSE"){
									formelone<-enviPat::mergeform(formelone,as.character(adducts[with_add[j]==adducts[,1],7]))
								}
								if(as.character(adducts[with_add[j]==adducts[,1],8])!="FALSE"){
									if(
										(enviPat::check_ded(formelone, as.character(adducts[with_add[j]==adducts[,1],8])))==FALSE
									){
										formelone<-enviPat::subform(formelone,as.character(adducts[with_add[j]==adducts[,1],8]))
									}else{
										next; # just skip that deduct!
									}
								}
								finform<-c(finform,formelone);
								fineidis<-c(fineidis,paste(eidis[i],"_",with_add[j],"_",tag1[i],"_",tag2[i],"_",tag3[i],sep=""));
								pattern<-enviPat::isopattern(
									isotopes,
									formelone,
									threshold=0.5,
									plotit=FALSE,
									charge=as.numeric(strsplit(as.character(adducts[with_add[j]==adducts[,1],3]),"+",fixed=TRUE)[[1]][1]),
									emass=0.00054858,
									algo=1
								)
								checked<-enviPat::check_chemform(isotopes, formelone)
								res<-try(enviPat::getR(checked, resolution_list[names(resolution_list) == logfile$parameters[22]][[1]], nknots = 7, spar = 0.1, plotit = FALSE), silent=TRUE)
								if(grepl("Error",res)){
									cat("\n Mass out of range of Resolution data - set to range!");
									if(checked[[3]]<=min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])){
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]
									}else{
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==max(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]									
									}
								}
								pattern<-enviPat::envelope(
									pattern,
									ppm=FALSE,
									dmz=0.0001,
									frac=1/4,
									env="Gaussian",
									resolution=res,
									plotit=FALSE
								)
								pattern<-enviPat::vdetect(
									pattern,
									detect="centroid",
									plotit=FALSE
								)
								pattern_neg_IS[[counter]]<-pattern[[1]];                    								
								patternRT_neg_IS<-c(patternRT_neg_IS,as.numeric(as.character(rets[i]))*60);
								if(tolrets[i]!="FALSE"){
									patternDelRT_neg_IS<-c(patternDelRT_neg_IS,as.numeric(as.character(tolrets[i]))*60)
								}else{
									patternDelRT_neg_IS<-c(patternDelRT_neg_IS,isolate(input$screen_IS_delRT))
								}
								if(rec[i]=="TRUE"){
									if(takeall){
										mz_neg_IS<-c(mz_neg_IS,pattern_neg_IS[[counter]][pattern_neg_IS[[counter]][,2]==max(pattern_neg_IS[[counter]][,2]),1][counter=1]);
										RT_neg_IS<-c(RT_neg_IS,
											as.numeric(as.character(rets[i]))*60
										);
									}else{
										if(j==1){
											mz_neg_IS<-c(mz_neg_IS,pattern_neg_IS[[counter]][pattern_neg_IS[[counter]][,2]==max(pattern_neg_IS[[counter]][,2]),1][counter=1]);
											RT_neg_IS<-c(RT_neg_IS,
												as.numeric(as.character(rets[i]))*60
											);
										}
									}
								}
								counter<-c(counter+1);
							}
						}
					}
					if(length(mz_neg_IS)>0){
						intmass_neg_IS<-data.frame(mz_neg_IS,RT_neg_IS)
						names(intmass_neg_IS)<-c("m/z","RT")
						save(intmass_neg_IS,file=file.path(logfile[[1]],"results","intmass_neg_IS"))                    
					}else{
						cat("\n no masses z=neg for recalibration formed!")
					}
					names(pattern_neg_IS)<-fineidis    
					save(pattern_neg_IS,file=file.path(logfile[[1]],"results","pattern_neg_IS"))
					save(patternRT_neg_IS,file=file.path(logfile[[1]],"results","patternRT_neg_IS"))
					save(patternDelRT_neg_IS,file=file.path(logfile[[1]],"results","patternDelRT_neg_IS"))
				}
				###################################################################

			}
            ######################################################################
				
				
				
				
            ######################################################################
			# on targets #########################################################
			targets<-read.table(file=file.path(logfile[[1]],"dataframes","targets.txt"),header=TRUE,sep="\t",colClasses = "character")
			if(length(targets$Formula[targets$Formula!="-"])>0){ # any target formulas available?

				##################################################################
				if(any(measurements[,4]=="positive") & any(targets[,8]=="positive") ){
					pattern_pos_target<-list(0);
					counter<-c(1);
					finform<-c();
					fineidis<-c();
					finadd<-c();
					mz_pos_target<-c();
					RT_pos_target<-c();
					patternRT_pos_target<-c();
					patternDelRT_pos_target<-c();
					formeln<-targets[,3];
					eidis<-targets[,1];
					rets<-targets[,4];
					tolrets<-targets[,5];
					rec<-targets[,9];
					tag1<-targets[,13];
					tag2<-targets[,14];
					tag3<-targets[,15];
					for(i in 1:length(formeln)){
						if(targets[i,8]=="positive"){
							cat(paste("\nDoing target positive ",as.character(i)," of ",as.character(length(formeln)),":",sep=""));
							if( (targets[i,11]=="TRUE") & (targets[i,7]=="FALSE") ){
								warning(paste("target pattern: restrict to adduct with no adduct selected? skipped target with ID:",eidis[i],"\n",sep=""))
								next;
							}
							takeall<-TRUE # use all adducts for recal?
							if(targets[i,11]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
								if(targets[i,7]=="FALSE"){ # no main adduct chosen?
									with_add<-c(logfile[[7]])
									takeall<-TRUE
								}else{ # with main adduct?
									with_add<-c(targets[i,7],logfile[[7]])
									with_add<-unique(with_add)
									takeall<-FALSE
								}
							}else{ # use main adduct only
								with_add<-c(targets[i,7])
								takeall<-TRUE # there is just one entry anyway
							}
							for(j in 1:length(with_add)){
								formelone<-formeln[i];
								formelone<-enviPat::multiform(formelone,as.numeric(adducts[with_add[j]==adducts[,1],9]))
								if(as.character(adducts[with_add[j]==adducts[,1],7])!="FALSE"){
									formelone<-enviPat::mergeform(formelone,as.character(adducts[with_add[j]==adducts[,1],7]))
								}
								if(as.character(adducts[with_add[j]==adducts[,1],8])!="FALSE"){
									if(
										(enviPat::check_ded(formelone, as.character(adducts[with_add[j]==adducts[,1],8])))==FALSE
									){
										formelone<-enviPat::subform(formelone,as.character(adducts[with_add[j]==adducts[,1],8]))
									}else{
										next; # just skip that deduct!
									}
								}
								finform<-c(finform,formelone);
								fineidis<-c(fineidis,paste(eidis[i],"_",with_add[j],"_",tag1[i],"_",tag2[i],"_",tag3[i],sep=""));
								pattern<-enviPat::isopattern(
									isotopes,
									formelone,
									threshold=0.5,
									plotit=FALSE,
									charge=as.numeric(strsplit(as.character(adducts[with_add[j]==adducts[,1],3]),"+",fixed=TRUE)[[1]][1]),
									emass=0.00054858,
									algo=1
								)
								checked<-enviPat::check_chemform(isotopes, formelone)
								res<-try(enviPat::getR(checked, resolution_list[names(resolution_list) == logfile$parameters[22]][[1]], nknots = 7, spar = 0.1, plotit = FALSE), silent=TRUE)
								if(grepl("Error",res)){
									cat("\n Mass out of range of Resolution data - set to range!");
									if(checked[[3]]<=min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])){
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]
									}else{
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==max(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]									
									}
								}
								pattern<-enviPat::envelope(
									pattern,
									ppm=FALSE,
									dmz=0.0001,
									frac=1/4,
									env="Gaussian",
									resolution=res,
									plotit=FALSE
								)
								pattern<-enviPat::vdetect(
									pattern,
									detect="centroid",
									plotit=FALSE
								)
								pattern_pos_target[[counter]]<-pattern[[1]];                    
								patternRT_pos_target<-c(patternRT_pos_target,as.numeric(as.character(rets[i]))*60);
								if(tolrets[i]!="FALSE"){
									patternDelRT_pos_target<-c(patternDelRT_pos_target,as.numeric(as.character(tolrets[i]))*60)
								}else{
									patternDelRT_pos_target<-c(patternDelRT_pos_target,isolate(input$screen_target_delRT))
								}
								if(rec[i]=="TRUE"){
									if(takeall){
										mz_pos_target<-c(mz_pos_target,pattern_pos_target[[counter]][pattern_pos_target[[counter]][,2]==max(pattern_pos_target[[counter]][,2]),1][counter=1]);
										RT_pos_target<-c(RT_pos_target,
											as.numeric(as.character(rets[i]))*60
										);
									}else{
										if(j==1){
											mz_pos_target<-c(mz_pos_target,pattern_pos_target[[counter]][pattern_pos_target[[counter]][,2]==max(pattern_pos_target[[counter]][,2]),1][counter=1]);
											RT_pos_target<-c(RT_pos_target,
												as.numeric(as.character(rets[i]))*60
											);
										}
									}
								}
								counter<-c(counter+1);
							}
						}
					}
					if(length(mz_pos_target)>0){
						intmass_pos_target<-data.frame(mz_pos_target,RT_pos_target)
						names(intmass_pos_target)<-c("m/z","RT")
						save(intmass_pos_target,file=file.path(logfile[[1]],"results","intmass_pos_target"))                    
					}else{
						cat("\n no masses z=pos for recalibration formed!")
					}
					names(pattern_pos_target)<-fineidis    
					save(pattern_pos_target,file=file.path(logfile[[1]],"results","pattern_pos_target"))
					save(patternRT_pos_target,file=file.path(logfile[[1]],"results","patternRT_pos_target"))
					save(patternDelRT_pos_target,file=file.path(logfile[[1]],"results","patternDelRT_pos_target"))
				}
				###################################################################
				
				##################################################################
				if(any(measurements[,4]=="negative")  & any(targets[,8]=="negative")  ){
					pattern_neg_target<-list(0);
					counter<-c(1);
					finform<-c();
					fineidis<-c();
					finadd<-c();
					mz_neg_target<-c();
					RT_neg_target<-c();
					patternRT_neg_target<-c();
					patternDelRT_neg_target<-c();
					formeln<-targets[,3];
					eidis<-targets[,1];
					rets<-targets[,4];
					tolrets<-targets[,5];
					rec<-targets[,9];
					tag1<-targets[,13];
					tag2<-targets[,14];
					tag3<-targets[,15];
					for(i in 1:length(formeln)){
						if(targets[i,8]=="negative"){
							cat(paste("\nDoing target negative ",as.character(i)," of ",as.character(length(formeln)),":",sep=""));
							if( (targets[i,11]=="TRUE") & (targets[i,7]=="FALSE") ){
								warning(paste("target pattern: restrict to adduct with no adduct selected? skipped target with ID:",eidis[i],"\n",sep=""))
								next;
							}
							takeall<-TRUE # use all adducts for recal?
							if(targets[i,11]=="FALSE"){ # use main adduct & all adducts - main adduct MUST come first -> for mz_ and RT_ below
								if(targets[i,7]=="FALSE"){ # no main adduct chosen?
									with_add<-c(logfile[[8]])
									takeall<-TRUE
								}else{ # with main adduct?
									with_add<-c(targets[i,7],logfile[[8]])
									with_add<-unique(with_add)
									takeall<-FALSE
								}
							}else{ # use main adduct only
								with_add<-c(targets[i,7])
								takeall<-TRUE # there is just one entry anyway
							}
							for(j in 1:length(with_add)){
								formelone<-formeln[i];
								formelone<-enviPat::multiform(formelone,as.numeric(adducts[with_add[j]==adducts[,1],9]))
								if(as.character(adducts[with_add[j]==adducts[,1],7])!="FALSE"){
									formelone<-enviPat::mergeform(formelone,as.character(adducts[with_add[j]==adducts[,1],7]))
								}
								if(as.character(adducts[with_add[j]==adducts[,1],8])!="FALSE"){
									if(
										(enviPat::check_ded(formelone, as.character(adducts[with_add[j]==adducts[,1],8])))==FALSE
									){
										formelone<-enviPat::subform(formelone,as.character(adducts[with_add[j]==adducts[,1],8]))
									}else{
										next; # just skip that deduct!
									}
								}
								finform<-c(finform,formelone);
								fineidis<-c(fineidis,paste(eidis[i],"_",with_add[j],"_",tag1[i],"_",tag2[i],"_",tag3[i],sep=""));
								pattern<-enviPat::isopattern(
									isotopes,
									formelone,
									threshold=0.5,
									plotit=FALSE,
									charge=as.numeric(strsplit(as.character(adducts[with_add[j]==adducts[,1],3]),"+",fixed=TRUE)[[1]][1]),
									emass=0.00054858,
									algo=1
								)
								checked<-enviPat::check_chemform(isotopes, formelone)
								res<-try(enviPat::getR(checked, resolution_list[names(resolution_list) == logfile$parameters[22]][[1]], nknots = 7, spar = 0.1, plotit = FALSE), silent=TRUE)
								if(grepl("Error",res)){
									cat("\n Mass out of range of Resolution data - set to range!");
									if(checked[[3]]<=min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])){
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==min(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]
									}else{
										use_this<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1]==max(resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][,1])
										res<-resolution_list[names(resolution_list) == logfile$parameters[22]][[1]][use_this,2][1]									
									}
								}
								pattern<-enviPat::envelope(
									pattern,
									ppm=FALSE,
									dmz=0.0001,
									frac=1/4,
									env="Gaussian",
									resolution=res,
									plotit=FALSE
								)
								pattern<-enviPat::vdetect(
									pattern,
									detect="centroid",
									plotit=FALSE
								)
								pattern_neg_target[[counter]]<-pattern[[1]];                    
								patternRT_neg_target<-c(patternRT_neg_target,as.numeric(as.character(rets[i]))*60);
								if(tolrets[i]!="FALSE"){
									patternDelRT_neg_target<-c(patternDelRT_neg_target,as.numeric(as.character(tolrets[i]))*60)
								}else{
									patternDelRT_neg_target<-c(patternDelRT_neg_target,isolate(input$screen_target_delRT))
								}								
								if(rec[i]=="TRUE"){
									if(takeall){
										mz_neg_target<-c(mz_neg_target,pattern_neg_target[[counter]][pattern_neg_target[[counter]][,2]==max(pattern_neg_target[[counter]][,2]),1][counter=1]);
										RT_neg_target<-c(RT_neg_target,
											as.numeric(as.character(rets[i]))*60
										);
									}else{
										if(j==1){
											mz_neg_target<-c(mz_neg_target,pattern_neg_target[[counter]][pattern_neg_target[[counter]][,2]==max(pattern_neg_target[[counter]][,2]),1][counter=1]);
											RT_neg_target<-c(RT_neg_target,
												as.numeric(as.character(rets[i]))*60
											);
										}
									}
								}
								counter<-c(counter+1);
							}
						}
					}
					if(length(mz_neg_target)>0){
						intmass_neg_target<-data.frame(mz_neg_target,RT_neg_target)
						names(intmass_neg_target)<-c("m/z","RT")
						save(intmass_neg_target,file=file.path(logfile[[1]],"results","intmass_neg_target"))                    
					}else{
						cat("\n no masses z=neg for recalibration formed!")
					}
					names(pattern_neg_target)<-fineidis    
					save(pattern_neg_target,file=file.path(logfile[[1]],"results","pattern_neg_target"))
					save(patternRT_neg_target,file=file.path(logfile[[1]],"results","patternRT_neg_target"))
					save(patternDelRT_neg_target,file=file.path(logfile[[1]],"results","patternDelRT_neg_target"))
				}
				###################################################################
	
			}
			####################################################################	
				
	############################################################################
	# (3) logfile entry 
	logfile$summary[4,2]<<-"TRUE";
	logfile$summary[4,2]<-"TRUE";
	logfile$Tasks_to_redo[8]<-"FALSE";	
	logfile$Tasks_to_redo[8]<<-"FALSE";	
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	summa[4,2]<-"done"
	summa[4,2]<<-"done"
	output$summar<<-renderTable(summa);
	cat("\nIsotope pattern calculations completed \n");
	output$dowhat<<-renderText("Isotope pattern calculations done ... wait")

}else{

	logfile$Tasks_to_redo[8]<-"FALSE";
	logfile$Tasks_to_redo[8]<<-"FALSE";
	save(logfile,file=file.path(as.character(logfile[[1]]),"logfile.emp"));
	summa[4,2]<-"skipped"
	summa[4,2]<<-"skipped"
	output$summar<<-renderTable(summa);
	cat("\nIsotope pattern calculations skipped \n");
	output$dowhat<<-renderText("Isotope pattern calculations skipped ... wait")

}




