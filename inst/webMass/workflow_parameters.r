########################
cat("\n Reading external parameters ...")
logfile$parameters$external<-list()
########################
# add custom parameters below, using the format:
# logfile$parameters$external$x<-"y";
# with:
# x = new variable name
# y = new variable values, set as character string; must be converted to aspired format in scripts, e.g.,
# as.numeric(logfile$parameters$external$x), as.logical(), ...
########################

# Parameters for EIC correlation 
logfile$parameters$external$EICor_delRT<-5		# [s] RT window for candidate peak pairs
logfile$parameters$external$EICor_minpeaks<-10 	# min. number of data points per EIC & shared in EIC pair to check for correlation - otherwise set to 1000
logfile$parameters$external$EICor_mincor<-.95 	# minimum correlation

# Parameters for Isotopologue grouping
logfile$parameters$external$isotop_mztol<-2.5
logfile$parameters$external$isotop_ppm<-TRUE
logfile$parameters$external$isotop_inttol<-0.5
logfile$parameters$external$isotop_rttol<-5
logfile$parameters$external$isotop_use_charges<-c(1,2)

# Parameters for adduct grouping
logfile$parameters$external$adducts_rttol<-5
logfile$parameters$external$adducts_mztol<-2.5
logfile$parameters$external$adducts_ppm<-TRUE
logfile$parameters$external$adducts_pos<-c("M+H","M+Na","M+K","M+NH4")
logfile$parameters$external$adducts_neg<-c("M-H","M-","2M-H")

# Parameters for homologue grouping
logfile$parameters$external$homol_units<-c("CH2","CH4O")
logfile$parameters$external$homol_charges<-c(1,2)
logfile$parameters$external$homol_minmz<-10
logfile$parameters$external$homol_maxmz<-120
logfile$parameters$external$homol_minrt<-10
logfile$parameters$external$homol_maxrt<-60
logfile$parameters$external$homol_ppm<-TRUE
logfile$parameters$external$homol_mztol<-2.5
logfile$parameters$external$homol_rttol<-20
logfile$parameters$external$homol_minlength<-6
logfile$parameters$external$homol_vec_size<-1E8

# Parameters for componentization




########################
cat("done.")
########################