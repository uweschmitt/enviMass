#' @title
#'
#' @export
#'
#' @description 
#'
#' @param profileList List.
#' @param prof_ID Integer.
#' @param links_profiles List.
#' @param min_peaks Integer. Minimum number of peaks co-occuring between profiles to check their correlation
#' @param skip_peaks Logical. TRUE: omit profile links which are below min_peaks AND min_cor. Else (FALSE), keep those which are below min_peaks OR above min_cor.
#' @param min_cor Numeric. Minimum spearman correlation threshold between profile pairs to retain.
#' @param with_test Logical. Run a few tests along?
#' @param only_direct Logical. Check only profiles which are directly linked to the one with prof_ID	
#' @param del_RT Numeric or logical FALSE. Maximum permitted RT deviation between profiles. Otherwise, set to FALSE
#' @param omit_profiles Logical. Logical FALSE or Integer vector of length of profiles. For the latter, entries other than 0 omit a profile (e.g., if it was grouped before).
#' 
#' @details enviMass workflow function. 
#' 

get_isotopol<-function(
		profileList,
		prof_ID,
		links_profiles,
		min_peaks=3,
		skip_peaks=FALSE,
		min_cor=.9,
		with_test=FALSE,
		only_direct=FALSE,
		del_RT=30,
		omit_profiles=FALSE
	){

	prof_isot_IDs<-c()
	################################################################
	if(omit_profiles[1]!="FALSE"){
		if(length(omit_profiles)!=length(profileList[["index_prof"]][,"links"])){stop("omit_profiles not equal to number of profiles")}
	}
	if(with_test){if(prof_ID>length(profileList[["index_prof"]][,"profile_ID"])){stop("\n Debug get_isotopol _1!")}}
	if(with_test){if(profileList[["index_prof"]][prof_ID,"profile_ID"]!=prof_ID){stop("\n Debug get_isotopol _2!")}}
	if(profileList[["index_prof"]][prof_ID,"links"]==0){ # no links?
		return(prof_ID) # only main profile
	}
	################################################################
	# collect direct isotopologue links of main peak ###############
	in_link<-profileList[["index_prof"]][prof_ID,"links"][[1]]	
	inter_links<-c()
	not_inter_links<-c()
	if(length(links_profiles[[in_link]]$isot[,1])>0){
		if(skip_peaks){
			those<-which(
				(links_profiles[[in_link]]$isot[,"ref"]>=min_cor) &
				(links_profiles[[in_link]]$isot[,"link counts"]>=min_peaks)
			)
			if(length(those)>0){
				those<-links_profiles[[in_link]]$isot[those,"linked profile"]
				if(omit_profiles[1]!="FALSE"){
					those<-those[omit_profiles[those]==0]
				}
			}
			not_those<-which(
				(links_profiles[[in_link]]$isot[,"ref"]<min_cor) |
				(links_profiles[[in_link]]$isot[,"link counts"]<min_peaks)
			)
			if(length(not_those)>0){
				not_those<-links_profiles[[in_link]]$isot[not_those,"linked profile"]
			}
		}else{
			those<-which(
				(links_profiles[[in_link]]$isot[,"ref"]>=min_cor) |
				(links_profiles[[in_link]]$isot[,"link counts"]<min_peaks)
			)		
			if(length(those)>0){
				those<-links_profiles[[in_link]]$isot[those,"linked profile"]
				if(omit_profiles[1]!="FALSE"){
					those<-those[omit_profiles[those]==0]
				}
			}
			not_those<-which(
				(links_profiles[[in_link]]$isot[,"ref"]<min_cor) &
				(links_profiles[[in_link]]$isot[,"link counts"]>min_peaks)
			)		
			if(length(not_those)>0){
				not_those<-links_profiles[[in_link]]$isot[not_those,"linked profile"]
			}			
		}
		if(length(those)==0){
			return(prof_ID) # only main profile
		}else{
			inter_links<-c(inter_links,those)
			not_inter_links<-c(not_inter_links,not_those)		
			if(only_direct){return(c(prof_ID,inter_links))} # omit below indirectly related links
		}
	}else{
		return(prof_ID) # only main profile
	}
	prof_isot_IDs<-c(prof_ID) # claimed all links for main profile
	# all below can be skipped if main profile has not enough peaks.
	################################################################
	# collect all indirect links ###################################
	if(skip_peaks & (profileList[["index_prof"]][prof_ID,"number_peaks_total"][[1]]<min_peaks)){
		return(prof_isot_IDs)
	} # all below can be skipped if main profile has not enough peaks.
	while(length(inter_links)>0){
		in_link<-profileList[["index_prof"]][inter_links[1],"links"][[1]]
		if(length(links_profiles[[in_link]]$isot[,1])>0){
			those<-links_profiles[[in_link]]$isot[,"linked profile"]
			# remove: finished or intermediate links ############
			those<-those[is.na(match(those,c(inter_links,prof_isot_IDs)))]
			# remove: rejected links ############################
			if(length(those)>0){
				those<-those[is.na(match(those,c(not_inter_links)))]
			}
			#
			if(length(those)>0 & omit_profiles[1]!="FALSE"){
				those<-those[omit_profiles[those]==0]		
			}
			# check RT difference to main profile ###############
			if(length(those)>0 & del_RT!=FALSE){
				those<-those[abs(profileList[["index_prof"]][prof_ID,"mean_RT"]-profileList[["index_prof"]][those,"mean_RT"])<=del_RT]
			}
			# correlate with main profile #######################
			if(length(those)>0){ # main profile has enough peaks?
				this<-profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"sampleIDs"]
				del1<-profileList[["index_prof"]][prof_ID,"number_peaks_total"][[1]]
				keep<-rep(FALSE,length(those))
				for(k in 1:length(those)){
					if(skip_peaks){ # del1 checked above				 	
						del2<-profileList[["index_prof"]][those[k],"number_peaks_total"][[1]]
						if(del2<min_peaks){next}						
						that<-profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"sampleIDs"]
						matched<-match(this,that)							
						if(sum(!is.na(matched))<min_peaks){next}
						int1<-((profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"intensity"])[!is.na(matched)])
						int2<-((profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"intensity"])[matched[!is.na(matched)]])					
						correl<-cor(int1,int2)
						if(correl>=min_cor){
							keep[k]<-TRUE
						}
					}else{
						del2<-profileList[["index_prof"]][those[k],"number_peaks_total"][[1]]
						if(del2<min_peaks){
							keep[k]<-TRUE
							next
						}								
						that<-profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"sampleIDs"]
						matched<-match(this,that)							
						if(sum(!is.na(matched))<min_peaks){
							keep[k]<-TRUE
							next
						}								
						int1<-((profileList[["peaks"]][profileList[["index_prof"]][prof_ID,"start_ID"]:profileList[["index_prof"]][prof_ID,"end_ID"],"intensity"])[!is.na(matched)])
						int2<-((profileList[["peaks"]][profileList[["index_prof"]][those[k],"start_ID"]:profileList[["index_prof"]][those[k],"end_ID"],"intensity"])[matched[!is.na(matched)]])					
						correl<-cor(int1,int2)
						if(correl>=min_cor){
							keep[k]<-TRUE
						}
					}
				}		
				# keep remaining links
				if(length(those[keep])>0){
					inter_links<-c(inter_links,those[keep])
				}				
				# save rejected links
				if(length(those[!keep])>0){
					not_inter_links<-c(not_inter_links,those[!keep])
				}			
			}
		}
		prof_isot_IDs<-c(prof_isot_IDs,inter_links[1])
		inter_links<-inter_links[-1]
	}
	################################################################	
	if(with_test){
		if(any(profileList[["index_prof"]][prof_isot_IDs,"profile_ID"]!=prof_isot_IDs)){
			stop("\n Debug get_isotopol _3!")
		}
	}
	################################################################
	return(prof_isot_IDs)
	
}
