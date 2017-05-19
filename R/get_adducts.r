#' @title
#'
#' @export
#'
#' @description Given any measurements dataframe, extracts the ID of the newest file by date&time
#'
#' @param measurements Dataframe listing available enviMass files
#' 
#' @details enviMass workflow function. 
#' 

get_adducts<-function(
		profileList,
		prof_ID,
		links_profiles,
		min_peaks=3,
		skip_peaks=FALSE,
		min_cor=.9,
		with_test=FALSE,
		omit_profiles=FALSE
	){

	################################################################
	if(omit_profiles[1]!="FALSE"){
		if(length(omit_profiles)!=length(profileList[["index_prof"]][,"links"])){stop("omit_profiles not equal to number of profiles")}
	}
	if(with_test){if(prof_ID>length(profileList[["index_prof"]][,"profile_ID"])){stop("\n Debug get_isotopol _1!")}}
	if(with_test){if(profileList[["index_prof"]][prof_ID,"profile_ID"]!=prof_ID){stop("\n Debug get_isotopol _2!")}}
	if(profileList[["index_prof"]][prof_ID,"links"]==0){ # no links?
		return(c()) # only main profile
	}
	################################################################
	# collect direct isotopologue links of main peak ###############
	in_link<-profileList[["index_prof"]][prof_ID,"links"][[1]]	
	if(length(links_profiles[[in_link]]$adduc[,1])>0){
		if(skip_peaks){
			those<-which(
				(links_profiles[[in_link]]$adduc[,"ref"]>=min_cor) &
				(links_profiles[[in_link]]$adduc[,"link counts"]>=min_peaks)
			)
			if(length(those)>0){
				those<-links_profiles[[in_link]]$adduc[those,"linked profile"]
				if(omit_profiles[1]!="FALSE"){
					those<-those[omit_profiles[those]==0]
				}
			}			
		}else{
			those<-which(
				(links_profiles[[in_link]]$adduc[,"ref"]>=min_cor) |
				(links_profiles[[in_link]]$adduc[,"link counts"]<min_peaks)
			)		
			if(length(those)>0){
				those<-links_profiles[[in_link]]$adduc[those,"linked profile"]
				if(omit_profiles[1]!="FALSE"){
					those<-those[omit_profiles[those]==0]
				}
			}
		}
		if(length(those)==0){
			prof_adduc_IDs<-c() # only main profile
		}else{
			prof_adduc_IDs<-those
		}
	}else{
		return(c()) 	
	}
	 # claimed all links for main profile
	################################################################	
	if(with_test){
		if(any(profileList[["index_prof"]][prof_adduc_IDs,"profile_ID"]!=prof_adduc_IDs)){
			stop("\n Debug get_isotopol _3!")
		}
	}
	################################################################
	return(prof_adduc_IDs)
	
}
