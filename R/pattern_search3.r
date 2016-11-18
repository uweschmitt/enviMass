#' @title 
#'
#' @description 
#'
#' @param
#' 
#' @details 
#' 


pattern.search3<-function(
	peaklist, 
	quantiz, 
	mztol = 2, 
	ppm = TRUE, 
	inttol = 0.5, 
    rttol = 0.3, 
	use_isotopes = c("13C", "37Cl", "15N", "81Br", "34S", "18O"), 
	use_charges = c(1, 2), 
	use_marker = TRUE, 
    quick = TRUE, 
	isotopes, 
	exclude=FALSE
){
    size_deltamass <- quantiz[[5]][[1]]
    size_mass <- quantiz[[5]][[2]]
    size_intens <- quantiz[[5]][[3]]
    adductmass_LB <- quantiz[[5]][[7]]
    adductmass_UB <- quantiz[[5]][[8]]
    max_d_delmz <- quantiz[[5]][[9]]
    max_d_ratio <- quantiz[[5]][[10]]
    max_d_mass <- quantiz[[5]][[11]]
    isotope_key <- quantiz[[3]]
    charge_key <- quantiz[[4]]
    size_deltamass <- (size_deltamass + max_d_delmz)
    size_mass <- (size_mass + max_d_mass)
    size_intens <- (size_intens + max_d_ratio)
    mass_slots <- quantiz[[7]]
    cat("\n(1) Check inputs ...")
    if (mztol < 0) {
        warning("mztol should be >=0!")
    }
    if (inttol > 1 || inttol < 0) {
        stop("inttol must be >0 and <=1")
    }
    if (!is.data.frame(peaklist)) {
        stop("peaklist must be a data.frame")
    }
    if (length(peaklist[1, ]) > 3) {
        stop("peaklist with > 3 columns not allowed")
    }
    if (!length(peaklist[, 1]) > 1) {
        stop("peaklist with one entry - doesn`t make sense ...")
    }
    if (!is.numeric(peaklist[, 1]) || !is.numeric(peaklist[, 
        2]) || !is.numeric(peaklist[, 3])) {
        stop("peaklist columns not numeric")
    }
    if (ppm == "TRUE") {
        ppm2 <- 1
    }else {
        ppm2 <- 0
    }
    if (use_isotopes[1] != "FALSE") {
        if (any(is.na(match(use_isotopes, isotope_key)))) {
            paste("invalid use_isotopes, available only: ", sep = "")
            print(isotope_key)
            stop()
        }
    }
    if (use_charges[1] != "FALSE") {
        use_charges2 <- abs(use_charges)
        if (any(is.na(match(use_charges2, charge_key)))) {
            paste("invalid use_charges, available only: ", sep = "")
            print(charge_key)
            stop()
        }
    }
    cat(" done.")
    cat("\n(2) Build peaklist kd-tree, screen slots, query quantized data: \n")
    pBar <- txtProgressBar(min = 0, max = length(peaklist[, 1]), 
        style = 3)
    inter <- as.numeric(interactive())
    peakTree <- .Call("kdtree4", as.matrix(peaklist[, 1:3]), 
        as.integer(inter), pBar, PACKAGE = "nontarget")
    close(pBar)
    peakTree <- peakTree[, 1:4, drop = FALSE]
    cat("\n screen ... ")
    mass_slots <- quantiz[[7]]
    int_slots <- (10^quantiz[[8]])
    pBar <- txtProgressBar(min = 0, max = length(peaklist[, 1]), 
        style = 3)
    relat <- .Call("peak_search", as.matrix(peaklist[, 1:3]), 
        as.matrix(peakTree), as.matrix(mass_slots), as.matrix(int_slots), 
        as.numeric(mztol), as.numeric(ppm2), as.numeric(inttol), 
        as.numeric(rttol), as.integer(inter), pBar, PACKAGE = "nontarget")
    close(pBar)
	#########################################################################################################	
	# DEBUG?
	if(any(duplicated(relat))){
		cat("\n WARNING: duplicated relat entries - debug me!")
		relat<-unique(relat)
	}
	#########################################################################################################
	# NEW
	if(exclude[1]!=FALSE){
		# exclude<-exclude[order(exclude[,1],exclude[,2],decreasing=FALSE),] # already ordered in do_EIC_correlation.r
		those<-(relat[,1]>relat[,2])
		relat[those,]<-relat[those,c(2,1)]
		relat<-relat[order(relat[,1],relat[,2],decreasing=FALSE),]
		found<-enviMass:::rows_compare(relat,exclude,row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=FALSE)
		cat("Exclusion: ")
		cat(paste(as.character(round(sum(found)/length(found)*100,digits=2)),"% of potential potential pairs removed by missing EIC correlation.",sep=""))
		relat<-relat[!found,]
	}else{
		cat("\n exclusion skipped.")
	}
	#########################################################################################################	
    if (length(relat) < 1) {
        return("\n No matches found \n ")
    }
    done <- matrix(ncol = length(charge_key), nrow = length(isotope_key), 
        FALSE)
    colnames(done) <- charge_key
    rownames(done) <- isotope_key
    search_bounds <- rep(0, 6)
    marker_bounds <- matrix(ncol = 2, nrow = 3, 0)
    marker_bounds[2, 1] <- min(peaklist[, 2])
    marker_bounds[2, 2] <- max(peaklist[, 2])
    bound_int <- log10((1 + inttol)/(1 - inttol))
    from_peak <- c()
    to_peak <- c()
    isotope <- c()
    charge <- c()
    retr_1 <- 0
    pBar <- txtProgressBar(min = 0, max = length(relat[, 1]), 
        style = 3)
    for (j in 1:length(relat[, 1])) {
        done[, ] <- FALSE
        got <- FALSE
        del_mass <- (peaklist[relat[j, 2], 1] - peaklist[relat[j, 
            1], 1])
        if (ppm == TRUE) {
            search_bounds[1] <- (del_mass - (2 * mztol * peaklist[relat[j, 
                1], 1]/1e+06))
            search_bounds[2] <- (del_mass + (2 * mztol * peaklist[relat[j, 
                1], 1]/1e+06))
        }
        else {
            search_bounds[1] <- (del_mass - (2 * mztol))
            search_bounds[2] <- (del_mass + (2 * mztol))
        }
        search_bounds[3] <- (peaklist[relat[j, 1], 1] + adductmass_LB)
        search_bounds[4] <- (peaklist[relat[j, 1], 1] + adductmass_UB)
        log_int <- log10(peaklist[relat[j, 1], 2]/peaklist[relat[j, 
            2], 2])
        search_bounds[5] <- (log_int - bound_int)
        search_bounds[6] <- (log_int + bound_int)
        for (i in 1:length(quantiz[[6]])) {
            do <- TRUE
            if (use_isotopes[1] != FALSE) {
                if (!any(use_isotopes == isotope_key[as.numeric(strsplit(names(quantiz[[6]])[i], 
                  "_")[[1]][1])])) {
                  do <- FALSE
                }
            }
            if (use_charges[1] != FALSE) {
                if (!any(use_charges2 == as.numeric(strsplit(names(quantiz[[6]])[i], 
                  "_")[[1]][2]))) {
                  do <- FALSE
                }
            }
            if (do == FALSE) {
                next
            }
            if (strsplit(names(quantiz[[6]])[i], "_")[[1]][3] == 
                "wo") {
                found <- .Call("search_boxtree", quantiz[[6]][[i]][, 
                  1:6], quantiz[[6]][[i]][, 16:20], as.numeric(search_bounds), 
                  as.integer(0), PACKAGE = "nontarget")
                retr_1 <- c(retr_1 + 1)
                if (found == -2) {
                  done[as.numeric(strsplit(names(quantiz[[6]])[i], 
                    "_")[[1]][1]), as.numeric(strsplit(names(quantiz[[6]])[i], 
                    "_")[[1]][2])] <- TRUE
                  from_peak <- c(from_peak, relat[j, 1])
                  to_peak <- c(to_peak, relat[j, 2])
                  isotope <- c(isotope, as.numeric(strsplit(names(quantiz[[6]])[i], 
                    "_")[[1]][1]))
                  charge <- c(charge, as.numeric(strsplit(names(quantiz[[6]])[i], 
                    "_")[[1]][2]))
                  got <- TRUE
                }
            }
            if (got & quick) 
                break
            if (got) 
                next
            if (strsplit(names(quantiz[[6]])[i], "_")[[1]][3] == 
                "w") {
                if (done[as.numeric(strsplit(names(quantiz[[6]])[i], 
                  "_")[[1]][1]), as.numeric(strsplit(names(quantiz[[6]])[i], 
                  "_")[[1]][2])] == FALSE) {
                  if (use_marker != "TRUE") {
                    found <- .Call("search_boxtree", quantiz[[6]][[i]][, 
                      1:6], quantiz[[6]][[i]][, 16:20], as.numeric(search_bounds), 
                      as.integer(0), PACKAGE = "nontarget")
                    if (found == -2) {
                      done[as.numeric(strsplit(names(quantiz[[6]])[i], 
                        "_")[[1]][1]), as.numeric(strsplit(names(quantiz[[6]])[i], 
                        "_")[[1]][2])] <- TRUE
                      from_peak <- c(from_peak, relat[j, 1])
                      to_peak <- c(to_peak, relat[j, 2])
                      isotope <- c(isotope, as.numeric(strsplit(names(quantiz[[6]])[i], 
                        "_")[[1]][1]))
                      charge <- c(charge, as.numeric(strsplit(names(quantiz[[6]])[i], 
                        "_")[[1]][2]))
                      got <- TRUE
                    }
                    retr_1 <- c(retr_1 + 1)
                  }
                  else {
                    found <- .Call("search_boxtree", quantiz[[6]][[i]][, 
                      1:6], quantiz[[6]][[i]][, 16:20], as.numeric(search_bounds), 
                      as.integer(1), PACKAGE = "nontarget")
                    retr_1 <- c(retr_1 + 1)
                    if (length(found) > 0) {
                      for (k in 1:length(found)) {
                        marker_delmass <- c((peaklist[relat[j, 
                          2], 1] - quantiz[[6]][[i]][found[k], 
                          8]), (peaklist[relat[j, 2], 1] - quantiz[[6]][[i]][found[k], 
                          7]))
                        if (ppm == TRUE) {
                          marker_bounds[1, 1] <- (min(marker_delmass) - 
                            (2 * mztol * peaklist[relat[j, 1], 
                              1]/1e+06))
                          marker_bounds[1, 2] <- (max(marker_delmass) + 
                            (2 * mztol * peaklist[relat[j, 1], 
                              1]/1e+06))
                        }
                        else {
                          marker_bounds[1, 1] <- (min(marker_delmass) - 
                            (2 * mztol))
                          marker_bounds[1, 2] <- (max(marker_delmass) + 
                            (2 * mztol))
                        }
                        max_int <- max(peaklist[relat[j, ], 2])
                        marker_bounds[2, 1] <- (max_int * (1 - 
                          inttol))
                        marker_bounds[3, 1] <- (peaklist[relat[j, 
                          1], 3] - rttol)
                        marker_bounds[3, 2] <- (peaklist[relat[j, 
                          1], 3] + rttol)
                        found_m <- .Call("search_kdtree", as.matrix(peaklist[, 
                          1:3]), as.matrix(peakTree), as.matrix(marker_bounds), 
                          PACKAGE = "nontarget")
                        if (length(found_m) > 1) {
                          done[as.numeric(strsplit(names(quantiz[[6]])[i], 
                            "_")[[1]][1]), as.numeric(strsplit(names(quantiz[[6]])[i], 
                            "_")[[1]][2])] <- TRUE
                          from_peak <- c(from_peak, relat[j, 
                            1])
                          to_peak <- c(to_peak, relat[j, 2])
                          isotope <- c(isotope, as.numeric(strsplit(names(quantiz[[6]])[i], 
                            "_")[[1]][1]))
                          charge <- c(charge, as.numeric(strsplit(names(quantiz[[6]])[i], 
                            "_")[[1]][2]))
                          got <- TRUE
                          break
                        }
                      }
                    }
                  }
                }
            }
            if (inter == 1) {
                setTxtProgressBar(pBar, j, title = NULL, label = NULL)
            }
            if (got & quick) 
                break
        }
    }
    close(pBar)
    if (length(from_peak) == 0) {
        stop("\n No matches found \n")
    }
	relat<-cbind(from_peak,to_peak)
	return(relat)
	
	############################################################################################################
	
	if(FALSE){
    use <- c(rep(1, length(to_peak)), rep(2, length(to_peak)))
    isotope <- c(isotope, isotope)
    charge <- c(charge, charge)
    to_peak2 <- c(to_peak, from_peak)
    from_peak2 <- c(from_peak, to_peak)
    use <- use[order(from_peak2, decreasing = FALSE)]
    isotope <- isotope[order(from_peak2, decreasing = FALSE)]
    charge <- charge[order(from_peak2, decreasing = FALSE)]
    to_peak2 <- to_peak2[order(from_peak2, decreasing = FALSE)]
    from_peak2 <- from_peak2[order(from_peak2, decreasing = FALSE)]
    groups <- .Call("metagroup", as.integer(from_peak2), as.integer(to_peak2), 
        PACKAGE = "nontarget")
    to_peak <- to_peak2[use == 1]
    from_peak <- from_peak2[use == 1]
    charge <- charge[use == 1]
    isotope <- isotope[use == 1]
    groups <- groups[use == 1]
    cat(paste("\n  ", length(to_peak), " of ", length(relat[, 
        1]), " candidate linkages accepted.", "\n", sep = ""))
    cat("(3) Create output ...")
    pattern <- list(0)
    alls <- length(peaklist[, 1])
    ID <- seq(1:alls)
    getit1 <- rep("0", alls)
    getit2 <- rep("none", alls)
    getit3 <- rep("0", alls)
    getit4 <- rep("0", alls)
    getit5 <- rep("0", alls)
    getit6 <- rep("0", alls)
    for (i in 1:length(from_peak)) {
        getit1[from_peak[i]] <- paste(getit1[from_peak[i]], as.character(to_peak[i]), 
            sep = "/")
        getit2[from_peak[i]] <- paste(getit2[from_peak[i]], quantiz[[3]][isotope[i]], 
            sep = "/")
        getit3[from_peak[i]] <- paste(getit3[from_peak[i]], "large", 
            sep = "/")
        getit4[from_peak[i]] <- paste(getit4[from_peak[i]], as.character(charge[i]), 
            sep = "/")
        getit5[from_peak[i]] <- paste(getit5[from_peak[i]], as.character(groups[i]), 
            sep = "/")
        getit5[to_peak[i]] <- paste(getit5[to_peak[i]], as.character(groups[i]), 
            sep = "/")
    }
    for (i in 1:alls) {
        if (getit1[i] != "0") {
            getit1[i] <- substr(getit1[i], 3, nchar(getit1[i]))
        }
        if (getit2[i] != "none") {
            getit2[i] <- substr(getit2[i], 6, nchar(getit2[i]))
        }
        if (getit3[i] != "0") {
            getit3[i] <- substr(getit3[i], 3, nchar(getit3[i]))
        }
        if (getit4[i] != "0") {
            getit4[i] <- substr(getit4[i], 3, nchar(getit4[i]))
        }
        if (getit5[i] != "0") {
            getit5[i] <- substr(getit5[i], 3, nchar(getit5[i]))
            those <- strsplit(getit5[i], "/", fixed = TRUE)[[1]]
            those <- unique(those)
            getit5[i] <- those[1]
            if (length(those) > 1) {
                for (j in 2:length(those)) {
                  getit5[i] <- paste(getit5[i], those[j], sep = "/")
                }
            }
        }
        if (getit6[i] != "0") {
            getit6[i] <- substr(getit6[i], 3, nchar(getit6[i]))
        }
    }
    grouped_peaks <- data.frame(peaklist, ID, getit5, getit6, 
        getit1, getit2, getit3, getit4, stringsAsFactors = FALSE)
    names(grouped_peaks) <- c(names(peaklist), "peak ID", "group ID", 
        "interaction level", "to ID", "isotope(s)", "mass tolerance", 
        "charge level")
    pattern[[1]] <- grouped_peaks
    parameters <- data.frame(-rttol, rttol, mztol, 0, ppm, inttol, 
        0, 0, adductmass_LB, adductmass_UB, size_deltamass, size_mass, 
        size_intens, stringsAsFactors = FALSE)
    names(parameters) <- c("rttol", "rttol", "mztol", "mzfrac", 
        "ppm", "inttol", "cutint", "deter", "adductmass_LB", 
        "adductmass_UB", "size_deltamass", "size_mass", "size_intens")
    pattern[[2]] <- parameters
    groupID <- unique(groups)
    groupID <- groupID[order(groupID)]
    peakIDs <- c()
    charge_group <- c()
    charge_count <- unique(charge)
    counted <- rep(0, length(charge_count))
    for (i in 1:length(groupID)) {
        those <- c(from_peak[groups == groupID[i]], to_peak[groups == 
            groupID[i]])
        those <- unique(those)
        get1 <- as.character(those[1])
        for (j in 2:length(those)) {
            get1 <- paste(get1, ",", those[j], sep = "")
        }
        peakIDs <- c(peakIDs, get1)
        charges <- charge[groups == groupID[i]]
        charges <- unique(charges)
        counted[charge_count == charges[1]] <- (counted[charge_count == 
            charges[1]] + 1)
        charge_level <- as.character(charges[1])
        if (length(charges) > 1) {
            for (j in 2:length(charges)) {
                counted[charge_count == charges[1]] <- (counted[charge_count == 
                  charges[j]] + 1)
                charge_level <- paste(charges[j], "/", charge_level, 
                  sep = "")
            }
        }
        charge_group <- c(charge_group, charge_level)
    }
    for (i in 1:length(groupID)) {
        groupID[i] <- paste("/", as.character(groupID[i]), "/", 
            sep = "")
    }
    grouping <- data.frame(groupID, peakIDs, charge_group, stringsAsFactors = FALSE)
    names(grouping) <- c("group ID", "peak IDs", "charge")
    pattern[[3]] <- grouping
    pattern[[4]] <- "no information"
    charge_count <- data.frame(charge_count, counted, stringsAsFactors = FALSE)
    names(charge_count) <- c("Charge level", "Counts")
    pattern[[5]] <- charge_count
    pattern[[6]] <- numeric(0)
    pattern[[7]] <- "no information"
    pattern[[8]] <- "no information"
    if (quick) {
        pattern[[9]] <- "no information using quick"
    }
    else {
        isos <- rep(isotope_key, length(charge_key))
        chrgs <- rep(charge_key, length(isotope_key))
        chrgs <- chrgs[order(chrgs)]
        incr_count <- rep(0, length(chrgs))
        group_count <- rep(0, length(chrgs))
        element <- c()
        for (i in 1:length(isos)) {
            element <- c(element, as.character(isotopes[isotopes[, 
                2] == isos[i], 1][1]))
        }
        for (i in 1:length(getit2)) {
            if (getit2[i] != "none") {
                get1 <- strsplit(getit2[i], "/")[[1]]
                get2 <- strsplit(getit4[i], "/")[[1]]
                for (j in 1:length(get1)) {
                  incr_count[isos == get1[j] & chrgs == get2[j]] <- (incr_count[isos == 
                    get1[j] & chrgs == get2[j]] + 1)
                }
            }
        }
        for (i in 1:length(peakIDs)) {
            peaks <- as.numeric(strsplit(peakIDs[i], ",")[[1]])
            get3 <- c()
            get4 <- c()
            for (j in 1:length(peaks)) {
                if (getit2[peaks[j]] != "none") {
                  get3 <- c(get3, strsplit(getit2[peaks[j]], 
                    "/")[[1]])
                  get4 <- c(get4, as.numeric(strsplit(getit4[peaks[j]], 
                    "/")[[1]]))
                }
            }
            get1 <- unique(data.frame(get3, get4))
            for (j in 1:length(get1[, 1])) {
                group_count[isos == get1[j, 1] & chrgs == get1[j, 
                  2]] <- (group_count[isos == get1[j, 1] & chrgs == 
                  get1[j, 2]] + 1)
            }
        }
        counts <- data.frame(isos, chrgs, incr_count, group_count, 
            element, stringsAsFactors = FALSE)
        names(counts) <- c("isotope", "charge", "peak counts", 
            "group counts", "element")
        counts <- counts[(counts[, 3] != 0 | counts[, 4] != 0), 
            ]
        pattern[[9]] <- counts
    }
    if (quick) {
        pattern[[10]] <- "no information"
    }
    else {
        pattern[[10]] <- as.character(unique(counts[, 5]))
    }
    pattern[[11]] <- use_charges
    pattern[[12]] <- "no information"
    names(pattern) <- c("Patterns", "Parameters", "Peaks in pattern groups", 
        "Atom counts", "Count of pattern groups", "Removals by rules", 
        "Number of peaks with pattern group overlapping", "Number of peaks per within-group interaction levels", 
        "Counts of isotopes", "Elements", "Charges", "Rule settings")
    cat(paste(" queries: ", retr_1, " - done.\n", sep = ""))
    return(pattern)
	}
}

