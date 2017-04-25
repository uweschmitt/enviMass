#' @title Converts profileList peak table to matrix
#'
#'
#' @description Converts profileList peak table to matrix of peak intensities. 
#'
#' @param profileList A profile list.
#' @param n_latest. NULL or Integer. If integer, number of latest file peaks to include.
#' @param n_profiles Integer. How many of the most intense profiles to include?
#' @param incl_profiles. NULL or Logical vector of length equal to that of the number of profiles, with FALSE indicating profiles to exclude.
#' @param only_sample. Logical. TRUE = should peaks of sample files (and not, e.g., blank files) be included only?
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
	n_profiles=1000,
    n_latest=NULL,
    incl_profiles=NULL,
    only_sample=TRUE,
    normalize=TRUE
){

    ############################################################################
    if(!profileList[[1]][[3]]){stop("profileList not profiled; aborted.")}
    len<-dim(profileList[["index_prof"]])[1]
    if(!is.null(incl_profiles[1])){
        if(length(incl_profiles)!=len){
            stop("If not set to NULL, incl_profiles must be logical and of length equal to that of the number of profiles.")
        }
    }else{
        incl_profiles<-rep(TRUE,len)
    }
    if(!is.logical(normalize)){
        stop("normalize must be logical")
    }
	############################################################################

    ############################################################################    
    # get maximum intensity of each profile ####################################
    max_int<-rep(0,len)
    for(i in 1:len){
        if(!incl_profiles[i]){next}
        max_int[profileList_pos[["index_prof"]][i,"profile_ID"]]<-
            max(profileList_pos[["peaks"]][
                    profileList_pos[["index_prof"]][i,"start_ID"]:profileList_pos[["index_prof"]][i,"end_ID"]
                ,"intensity"]
            )
    }
    max_int_ord<-order(max_int,decreasing=TRUE)
    ############################################################################    

    ############################################################################    
    # get sample IDs ###########################################################
    ord<-order(as.POSIXct(profileList[["datetime"]]),decreasing=TRUE)
    if(!is.null(n_latest)[1]){
        ord<-ord[1:n_latest]
    }
    type<-profileList[["type"]][ord]
    file_ID<-profileList[["sampleID"]][ord]
    if(only_sample){
        file_ID<-file_ID[type=="sample"]
    }
    ############################################################################ 

    ############################################################################    
    # write to matrix ##########################################################
    mat<-matrix(nrow=length(file_ID),ncol=n_profiles,0)
    rownames(mat)<-as.character(file_ID)
    # which sample profiles to include -> most intense & not excluded?
    max_int_ord_nb<-max_int_ord[incl_profiles[max_int_ord]]
    max_int_ord_nb<-max_int_ord_nb[1:n_profiles]
    colnames(mat)<-as.character(max_int_ord_nb)
    sub_peaks_ind<-(
        !is.na(
            match(
                profileList_pos[["peaks"]][,"sampleIDs"],
                as.numeric(rownames(mat))
            )
        ) &
        !is.na(
            match(
                profileList_pos[["peaks"]][,"profileIDs"],
                as.numeric(colnames(mat))
            )
        )
    )
    mat[
        cbind(
            as.character(profileList_pos[["peaks"]][sub_peaks_ind,"sampleIDs"]),
            as.character(profileList_pos[["peaks"]][sub_peaks_ind,"profileIDs"])
        )
    ]<-profileList_pos[["peaks"]][sub_peaks_ind,"intensity"]
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
