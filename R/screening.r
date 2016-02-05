#' @title Screens a peak- and blanklist for matches with compounds
#'
#' @export
#'
#' @description \code{screening} screens a peak- and blanklist for matches with the isotope patterns of compounds
#'
#' @param peaklist Matrix or data.frame of sample peaks. 3 columns: m/z, intensity, RT.
#' @param blanklist Logical FALSE or matrix or data.frame of blind/blank peaks. 3 columns: m/z, intensity, RT.
#' @param pattern List of isotope patterns (i.e., centroids) of compounds to be screened for.
#' @param RT Numeric vector. Retention time to be matched; units equal to those of the input files. Must match length of \code{pattern}.
#' @param dmz Numeric. +/- m/z tolerance (precision)
#' @param ppm Logical. \code{dmass} given in ppm?
#' @param dRT Numeric. RT tolerance for matches with \code{peaklist}; units equal to those of the input files
#' @param dRTwithin Numeric. More restrictive RT tolerance relative to the most intense match.
#' @param dRTblank Numeric. RT tolerance window for matches with \code{blanklist}; units equal to those of the input files
#' @param dInt Numeric. Intensity tolerance, given as a fraction of a peak intensity.
#' @param Intcut Numeric. Centroids below this intensities are not screened for.
#' @param w1 0<=Numeric<=1. Weight for fraction of expected vs. matched peaks within mass tolerances.
#' @param w2 0<=Numeric<=1. Weight for fraction of expected vs. matched peaks within mass & intensity tolerances.
#' @param w3 0<=Numeric<=1. Weight for fraction of expected vs. matched peaks in \code{blanklist}.
#'
#' @return List of matches with length equal to \code{pattern}
#' 
#' @details  enviMass workflow function. Each entry in the list \code{pattern} refers to the centroid masses and intensities of a specific compound.
#' Hence, different adducts of the same compound result in different entries in this list. \code{w1,w2,w3} must sum to 1.
#' First, and for each entry in \code{pattern}, a match with the highest centroid is screened for. Having detected a match,
#' the intensities of the remaining theretical centroids from \code{pattern} are scaled accordingly. Those remaining above
#' \code{Intcut} are then screened for.
#'

screening<-function (	peaklist, 
						blanklist=FALSE, 
						pattern, 
						RT = 1000,
						dmz = 10, 
						ppm = TRUE, 
						dRT = 500, 
						dRTwithin = 50, 
						dRTblank = 200, 
						dInt = 0.2, 
						Intcut = 5000, 
						w1 = 0.8, 
						w2 = 0.2, 
						w3 = 0.0				
					){

	############################################################################
	# check arguments ##########################################################
	if(any(ls()=="logfile")){stop("illegal logfile detected #1 in screening.r!")}
    options(digits = 10)
	cat("\n Check inputs ...")
    if (is.list(pattern) == FALSE) {
        stop("Argument pattern not a list")
    }
    if (length(pattern) != length(RT) & length(RT)!=1) {
        stop("Length RT must equal length of pattern list!")
    }
    if (ppm != TRUE & ppm != FALSE) {
        stop("Wrong ppm setting!")
    }
    if (sum(w1, w2, w3) != 1) {
        stop("w1, w2, w3 must sum to 1!")
    }
    if (length(blanklist) == 1 & w3 != 0) {
        stop("No blanklist but w3 not equal to 0?")
    }
    if (length(dRT) != 1 & length(dRT) != length(pattern)) {
        stop("Invalid length of vector dRT")
    }
    if (length(dRTwithin) != 1 & length(dRTwithin) != length(pattern)) {
        stop("Invalid length of vector dRTwithin")
    }
    if (length(dInt) != 1 & length(dInt) != length(pattern)) {
        stop("Invalid length of vector dInt")
    }
    if (length(dRTblank) != 1 & length(dRTblank) != length(pattern)) {
        stop("Invalid length of vector dRTblank")
    }
    if (length(Intcut) != 1 & length(Intcut) != length(pattern)) {
        stop("Invalid length of vector Intcut")
    }
    if (length(dmz) != 1 & length(dmz) != length(pattern)) {
        stop("Invalid length of vector dmz")
    }
    if (length(dmz) == 1) {
        dmz <- rep(as.numeric(dmz), length(pattern))
    }
    if (length(dRT) == 1) {
        dRT <- rep(as.numeric(dRT), length(pattern))
    }
    if (length(RT) == 1) {
        RT <- rep(as.numeric(RT), length(pattern))
		cat("RT concatenated\n ")
    }
    if (length(dRTwithin) == 1) {
        dRTwithin <- rep(as.numeric(dRTwithin), length(pattern))
    }
    if (length(dInt) == 1) {
        dInt <- rep(as.numeric(dInt), length(pattern))
    }
    if (length(Intcut) == 1) {
        Intcut <- rep(as.numeric(Intcut), length(pattern))
    }
    cat(" done.")
	############################################################################
    # Find most intense peak ###################################################
	cat("\n Screen most intense ...")
    #for (i in 1:length(pattern)) {
    #    pattern[[i]] <- pattern[[i]][order(pattern[[i]][, 2], 
    #        decreasing = TRUE), ,drop=FALSE]
    #}
    results <- list(0)
    result <- data.frame( "-", "-", "-", "-", "-", "-", "-", 
        "-", "-", "-")
    names(result) <- c("max_hit", "all_hits", "dm/z", 
        "dRT", "dInt", "score_1", "score_2", "score_3", "score_sum", "below_Intcut")
    for (i in 1:length(pattern)) {
        results[[i]] <- result
    }
    names(results) <- names(pattern)
    mon <- rep()
    for (i in 1:length(pattern)) {
        mon <- c(mon, pattern[[i]][1,1])
    }	
    getit <- search_peak(peaklist, mon, dmz, ppm, RT, dRT)	
    cat(" done.")
	############################################################################
    # Full pattern screening ###################################################	
    cat("\n Full pattern screen ...")
	from <- c(1:length(pattern))
    getpeaks <- seq(1, length(peaklist[, 1]), 1)
    for ( i in 1:length(getit) ) {
        if ( getit[i] != FALSE ) {
            hihit <- as.numeric(strsplit(getit[i], " / ")[[1]])
            pat <- pattern[[from[i]]][, c(1, 2),drop=FALSE]
            for (j in 1:length(hihit)) {
                pat2 <- as.data.frame(pat)
                pat2[, 2] <- pat2[, 2]/max(pat2[, 2])
                pat2[, 2] <- pat2[, 2] * peaklist[hihit[j], 2]
                if (any(pat2[, 2] >= Intcut[from[i]])) {
					pat3 <- pat2[pat2[, 2] < Intcut[from[i]],]
					pat2 <- pat2[pat2[, 2] >= Intcut[from[i]],]
					many <- length(pat2[, 1])
                } else {
					pat2 <- pat2[1, ]
					pat3 <- pat2[rep(FALSE,length(pat2[, 1])), ]
					many <- 0
                }
                scoreitsample <- list(0)
                scoreitblank <- list(0)
                for (m in 1:length(pat2[, 1])) {
					if (ppm == TRUE) {
						mztol <- c(dmz[from[i]] * pat2[m, 1]/1e+06 * 2)
					} else {
						mztol <- c(dmz[from[i]] * 2)
					}
					these <- getpeaks[
						peaklist[, 1] >= (pat2[m,1] - mztol) & 
						peaklist[, 1] <= (pat2[m, 1] + mztol) & 
						peaklist[, 3] >= (peaklist[hihit[j],3] - dRTwithin[from[i]]) & 
						peaklist[, 3] <= (peaklist[hihit[j], 3] + dRTwithin[from[i]])]
					scoreitsample[[m]] <- these;
					if (length(blanklist) != 1) {
						these <- getpeaks[
							blanklist[,1] >= (pat2[m,1] - mztol) & 
							blanklist[, 1] <= (pat2[m,1] + mztol) & 
							blanklist[, 3] >= (RT[from[i]] - dRTblank[from[i]]) & 
							blanklist[, 3] <= (RT[from[i]] + dRTblank[from[i]]) ]
						scoreitblank[[m]] <- these;
					}else{
						scoreitblank[[m]] <- numeric(0)
					}
                }
                scoreitsample[[1]] <- hihit[j]
				get_other_peaks <- list(0)		
				found<-FALSE	
				if(length(pat3[,1])>0){
					for (m in 1:length(pat3[, 1])) {
						if (ppm == TRUE) {
							mztol <- c(dmz[from[i]] * pat3[m, 1]/1e+06 * 2)
						} else {
							mztol <- c(dmz[from[i]] * 2)
						}
						these <- getpeaks[
							peaklist[, 1] >= (pat3[m,1] - mztol) & 
							peaklist[, 1] <= (pat3[m, 1] + mztol) & 
							peaklist[, 3] >= (peaklist[hihit[j],3] - dRTwithin[from[i]]) & 
							peaklist[, 3] <= (peaklist[hihit[j], 3] + dRTwithin[from[i]])]
						get_other_peaks[[m]] <- these;
						if(length(these)>0){
							found<-TRUE	
						}
					}				
				}
				if(found){
					below_cut<-""
					for(m in 1:length(get_other_peaks)){
						if(length(get_other_peaks[[m]])>0){
							for(n in 1:length(get_other_peaks[[m]])){
								below_cut<-paste(below_cut,get_other_peaks[[m]][n],"/",sep="")
							}
						}
					}
				}else{
					below_cut<-"-"
				}
                if (length(scoreitsample) > 1) {
                  for (m in 2:length(scoreitsample)) {
                    if (length(scoreitsample[[m]]) > 1) {
                      scoreitsample[[m]] <- scoreitsample[[m]][abs(peaklist[hihit[j], 
                        3] - peaklist[scoreitsample[[m]], 3]) == 
                        min(abs(peaklist[hihit[j], 3] - peaklist[scoreitsample[[m]], 
                          3]))]
                    }
                    if (length(scoreitsample[[m]]) > 1) {
                      scoreitsample[[m]] <- scoreitsample[[m]][abs(pat2[m, 
                        1] - peaklist[scoreitsample[[m]], 1]) == 
                        min(abs(pat2[m, 1] - peaklist[scoreitsample[[m]], 
                          1]))]
                    }
                    if (length(scoreitsample[[m]]) > 1) {
                      scoreitsample[[m]] <- scoreitsample[[m]][abs(pat2[m, 
                        2] - peaklist[scoreitsample[[m]], 2]) == 
                        min(abs(pat2[m, 2] - peaklist[scoreitsample[[m]], 
                          2]))]
                    }
                    if (length(scoreitsample[[m]]) > 1) {
                      scoreitsample[[m]] <- sample(scoreitsample[[m]], 
                        1)
                    }
                  }
                }
                many_1 <- c(0)
                many_2 <- c(0)
                many_3 <- c(0)
                deltamz <- ""
                deltaRT <- ""
                deltaInt <- ""
                all_hits <- ""
                for (m in 1:length(scoreitsample)) {
					if (length(scoreitsample[[m]]) > 0) {
						many_1 <- c(many_1 + 1)
						all_hits <- paste(all_hits, scoreitsample[[m]], 
						  "/", sep = "")
						deltamz <- paste(deltamz, as.character(round(((pat2[m, 
						  1] - peaklist[scoreitsample[[m]], 1])/peaklist[scoreitsample[[m]], 
						  1] * 1e+06), digits = 2)), "/", sep = "")
						deltaRT <- paste(deltaRT, as.character(round((RT[from[i]] - 
						  peaklist[scoreitsample[[m]], 3]), digits = 2)), 
						  "/", sep = "")
						delInt <- round(((pat2[m, 2] - peaklist[scoreitsample[[m]], 
						  2])/peaklist[scoreitsample[[m]],2]), digits = 2)
						if( abs(delInt) <= (2*dInt[from[i]]) ) {
						  many_2 <- c(many_2 + 1)
						}
						deltaInt <- paste(deltaInt, as.character(delInt), 
						  "/", sep = "")
						if (length(scoreitblank[[m]]) == 0) {
						  many_3 <- c(many_3 + 1)
						}
					}else{
						all_hits <- paste(all_hits, "-/", sep="");
						deltamz <- paste(deltamz, "-/", sep="");
						deltaRT <- paste(deltaRT, "-/", sep="");
						deltaInt <- paste(deltaInt, "-/", sep="");
					}
                }				
                score1 <- paste(many_1, " of ", many, sep = "")
                score2 <- paste(many_2, " of ", many, sep = "")
                score3 <- paste(many_3, " of ", many, sep = "")
                sumscore <- as.character(round(c((many_1/many * 
                  w1) + (many_2/many * w2) + (many_3/many * w3)), 
                  digits = 3))
                result2 <- data.frame(as.character(hihit[j]), 
                  all_hits, deltamz, deltaRT, deltaInt, score1, 
                  score2, score3, sumscore, below_cut)
                names(result2) <- c("max_hit", "all_hits", 
                  "dm/z", "dRT", "dInt", "score_1", "score_2", 
                  "score_3", "score_sum", "below_Intcut")
                results[[from[i]]] <- rbind(results[[from[i]]],result2)
            }
        }
    }
    cat(" done.")
	############################################################################
    cat("\n Generate outputs ...")
    for (i in 1:length(results)) {
        if (length(results[[i]][, 1]) == 1) {
            results[[i]] <- "No most intense hit found: no isotopic pattern fit conducted."
        }else {
            results[[i]] <- results[[i]][-1, ]
			results[[i]]<-results[[i]][order(as.numeric(as.character(results[[i]][,9])),decreasing=TRUE),];				
        }
    }
    cat(" done.\n")
	############################################################################
    return(results)
	
}


