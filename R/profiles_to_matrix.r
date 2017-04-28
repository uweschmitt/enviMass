#' @title Converts profileList peak table to matrix
#'
#'
#' @description Converts profileList peak table to matrix of peak intensities. 
#'
#' @param profileList A profile list.
#' @param n_profiles Integer. How many of the most intense profiles to include?
#' @param only_sample. Logical. TRUE = should peaks of sample files (and not, e.g., blank files) be included only?
#' @param n_latest. NULL or Integer. If integer, number of latest file peaks to include.
#' @param normalize. Logical. TRUE = should intensities of each profile be divided by their maximum to range in [0,1]?
#'
#' @return Matrix with intensities.
#' 
#' @details The function writes the sample peak intensities of the most intense (= maximum peak intensity per profile) n_profiles 
#' into a non-sparse matrix. Columns contain profiles, rows sample peaks. Missing peak intensities are set to 0.
#' Rows are ordered by date and time, columns by decreasing amximum profile intensity.  
#' 

profiles_to_matrix<-function(
	profileList,
	n_profiles=NULL,
    only_sample=FALSE,
    n_latest=NULL,
    normalize=FALSE
){

	############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
    len<-dim(profileList[["index_prof"]])[1]
    if(!is.logical(normalize)){stop("normalize must be logical")}
	############################################################################

    ############################################################################    
    # get sample IDs ###########################################################
    ord<-order(as.POSIXct(profileList[["datetime"]]),decreasing=TRUE)
    type<-profileList[["type"]][ord]
    file_ID<-profileList[["sampleID"]][ord]
    if(!is.null(n_latest)[1]){
        type<-type[1:n_latest]
		file_ID<-file_ID[1:n_latest]
    }
	keep<-rep(TRUE,len)
	profile_IDs<-profileList[["index_prof"]][,"profile_ID"]
    if(only_sample){
        file_ID<-file_ID[type=="sample"]
		# filter out profiles which do not contain sample peaks
		for(i in 1:len){
			if(
				all(is.na(match(
					profileList[["peaks"]][
							profileList[["index_prof"]][i,"start_ID"]:profileList[["index_prof"]][i,"end_ID"]
						,"sampleIDs"],
					file_ID
				)))
			){
				keep[i]<-FALSE
			}
		}	
	}
    ############################################################################ 

    ############################################################################    
    # get maximum intensity of each profile ####################################
    max_int<-rep(0,len)
    for(i in 1:len){
        if(!keep[i]){next}
        max_int[i]<-
            max(profileList[["peaks"]][
                    profileList[["index_prof"]][i,"start_ID"]:profileList[["index_prof"]][i,"end_ID"]
                ,"intensity"]
            )
    }
	profile_IDs<-profile_IDs[keep]
	max_int<-max_int[keep]
    max_int_ord<-order(max_int,decreasing=TRUE)
	profile_IDs<-profile_IDs[max_int_ord]
    ############################################################################    

    ############################################################################    
    # write to matrix ##########################################################	
	if(is.null(n_profiles)){
		n_profiles<-length(profile_IDs)
	}else{
		if(n_profiles>length(profile_IDs)){
			n_profiles<-length(profile_IDs)
		}
	}
    mat<-matrix(nrow=length(file_ID),ncol=n_profiles,0)
    rownames(mat)<-as.character(file_ID)
	colnames(mat)<-as.character(profile_IDs[1:n_profiles])
    sub_peaks_ind<-(
        !is.na(
            match(
                profileList[["peaks"]][,"sampleIDs"],
                as.numeric(rownames(mat))
            )
        ) &
        !is.na(
            match(
                profileList[["peaks"]][,"profileIDs"],
                as.numeric(colnames(mat))
            )
        )
    )
	mat[
        cbind(
            match(as.character(profileList[["peaks"]][sub_peaks_ind,"sampleIDs"]),rownames(mat)),
            match(as.character(profileList[["peaks"]][sub_peaks_ind,"profileIDs"]),colnames(mat))
        )
    ]<-profileList[["peaks"]][sub_peaks_ind,"intensity"]
    ############################################################################
	
    ############################################################################
	# checks ###################################################################
	if(any(is.na(mat))){ 
		stop("Sth went wrong in function profiles_to_matrix - debug_1!")
	}	
	if(any(apply(mat,2,max)==0)){ 
		stop("Sth went wrong in function profiles_to_matrix - debug_2!")
	}
    ############################################################################ 	

    ############################################################################
    # (0,1)-normalize ##########################################################
    if(normalize){
        mat<-sweep(mat,2,apply(mat,2,max),"/")
    }
    ############################################################################

    ############################################################################
    return(mat)
	
}
