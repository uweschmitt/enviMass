##############################################################################################
# check upper atom count bounds ##############################################################
##############################################################################################
data(chemforms)


	masses<-centro[[1]][,1]
	intensities<-centro[[1]][,2]
	int_cut<-(max(intensities)*0.1)
	inttol=0.2
	use_C=FALSE
	charges=c(1,2)
	ppm=TRUE
	dmz=c(20,20,20,20,3,3,0)
	elements=c("C","H","N","O","Cl","Br","P","S")
	must_peak=FALSE

	
##############################################################################################
for(j in 1:20){ # different settings
	for(i in 1:length(chemforms)){
		use_chem<-chemforms[i]
		counts<-check_chemform(isotopes,use_chem,get_sorted=FALSE,get_list=TRUE)[[1]]
		at_charge<-sample(1:4,1)
		#at_charge<-3
		pattern<-isopattern(
		  isotopes,
		  use_chem,
		  threshold=0.1,
		  plot_it=FALSE,
		  charge=at_charge,
		  emass=0.00054858,
		  algo=2
		)
		res<-sample(seq(1E4,5E5,1E4),1)
		profiles<-envelope(
			pattern,
			ppm=FALSE,
			dmz=0.0001,   
			frac=1/10,
			env="Gaussian",
			resolution=1E5,
			plot_it=FALSE
		)
		centro<-vdetect(
		  profiles,
		  detect="centroid",
		  plot_it=FALSE
		)
		###########################################################################################
		# estimate bounds #########################################################################
		bounds<-enviMass:::atoms(
			masses=centro[[1]][,1],
			intensities=centro[[1]][,2],
			elements=names(counts),
			dmz=rep(100,length(names(counts))),
			ppm=TRUE,
			charges=at_charge,
			isotopes,
			int_cut=(max(centro[[1]][,2])*0.1),
			inttol=0.2,
			use_C=TRUE,
			must_peak=FALSE
		)	
		for(j in 1:length(counts)){
			here<-which(colnames(bounds)==names(counts)[j])
			if(counts[j]>max(bounds[,here])){
				stop("")
			}
		}
		
	}
}




