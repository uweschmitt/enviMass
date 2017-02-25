#' @title 
#'
#' @description 
#'
#' @param masses Masses
#' @param intensities 
#' @param elements 
#' @param dmz 
#' @param ppm 
#' @param charges
#' @param isotopes
#' @param int_cut
#' @param inttol
#' @param use_C
#' @param must_peak Must a non-monoisotopic peak be present to have an atom count >0 for an element? Non applicable to single-isotopic elements.
#' 
#' @details enviMass workflow function
#' 

atoms <-
function(
	masses,
	intensities,
	elements=c("C","H","N","O","Cl","Br","P"),
	dmz=c(20,20,20,20,3,3,0),
	ppm=TRUE,
	charges=c(1,2),
	isotopes,
	int_cut,
	inttol=0.2,
	use_C=FALSE,
	must_peak=FALSE
){

    ############################################################################
    # check inputs #############################################################
	if(length(masses)!=length(intensities)){stop("\n Vector masses must be of same length as intensities.")}
	if(length(elements)!=length(dmz)){stop("\n Vector dmz must be of same length as elements vector.")}	
	if(!is.logical(ppm)){stop("\n Argument ppm must be logical.")}
	if(!is.logical(use_C)){stop("\n Argument use_C must be logical.")}
	if(any(is.na(match(elements,isotopes[,1])))){stop("\n elements not found in isotopes table.")}
	if(!is.numeric(dmz)){stop("Argument dmz must be numeric.")}
	if(!is.numeric(charges)){stop("Argument charges must be numeric.")}	
	if(!is.numeric(int_cut)){stop("Argument int_cut must be numeric.")}		
	if(int_cut==0){stop("Argument int_cut must be >0")}		
	int_cut_up<-(int_cut+(int_cut*inttol))
	ord<-order(masses,decreasing=FALSE)
	masses<-masses[ord]
	intensities<-intensities[ord]
	max_intensities<-(intensities+(intensities*inttol))
	min_intensities<-(intensities-(intensities*inttol))	
    ############################################################################	
	bounds<-matrix(ncol=length(elements),nrow=length(charges),0) # default: zero atoms
	colnames(bounds)<-elements
	rownames(bounds)<-as.character(charges)
    ############################################################################
	# restrict by pattern ######################################################
	for(i in 1:length(elements)){
		at_isotopes<-which(as.character(isotopes[,1])==elements[i])
		w_C<-as.numeric(unique(isotopes[at_isotopes,5]))
		# single-isotopic element = estimate by monoisotopic mass alone ########	
		if( length(at_isotopes)==1 & !must_peak ){
			for(n in 1:length(charges)){ # make entry for all charges
				iso_mass<-as.numeric(isotopes[at_isotopes,3])
				if(!use_C || (w_C==0)){
					bounds[n,i]<-floor(masses[1]/(iso_mass[1]/charges[n]))
				}else{	
					iso_mass<-(iso_mass+(1/w_C*12))
					bounds[n,i]<-floor(masses[1]/(iso_mass[1]/charges[n]))
				}
			}
			next; # no pattern to be evaluated 
		}
		# non-single isotopes? check non-monoisotopic ones ######################
		at_isotopes<-at_isotopes[order(as.numeric(isotopes[at_isotopes,3]),decreasing=FALSE)]
		# set mass bounds 
		if(ppm){
			min_masses<-(masses-(dmz[i]*masses/1E6))
			max_masses<-(masses+(dmz[i]*masses/1E6))
		}else{
			min_masses<-(masses-(dmz[i]))
			max_masses<-(masses+(dmz[i]))		
		}		
		for(n in 1:length(charges)){	
			max_over_isotopes<-c() 		# get maximum atom count for each isotope
			peaks_over_isotopes<-c()	# indicate if an isotopic peak was found
			for(j in 2:length(at_isotopes)){ # only non-monoisotopic isotopes
				if(as.numeric(isotopes[at_isotopes[j],4])==0){next} # must omit zero-abundance isotopes
				iso_del_mass<-((as.numeric(isotopes[at_isotopes[j],3])-as.numeric(isotopes[at_isotopes[1],3]))/charges[n])
				iso_rat_abund<-(as.numeric(isotopes[at_isotopes[1],4])/as.numeric(isotopes[at_isotopes[j],4]))
				for_mass<-which( (masses>=(min_masses[1]+iso_del_mass)) & (masses<=(max_masses[1]+iso_del_mass)) )
				if(length(for_mass)==0){ # no match on first M+X position
					# get estimate from int_cut_up:					
					n_atom<-floor(int_cut_up/min_intensities[1]*iso_rat_abund)					
					max_over_isotopes<-c(max_over_isotopes,n_atom)					
					peaks_over_isotopes<-c(peaks_over_isotopes,FALSE)
				}else{ # match : solve full pattern iteratively, starting from maximum estimate set by M+X position
					# get most intense mass
					for_mass<-for_mass[intensities[for_mass]==max(intensities[for_mass])][1] 
					# derive initial upper bound
					n_atom<-floor(max_intensities[for_mass]/min_intensities[1]*iso_rat_abund)
					if(n_atom>0){
						# iterate isotope pattern
						for(w in n_atom:0){
							found_peak<-FALSE
							if(w==0){
								n_atom_iter<-0
								break;
							}
							# derive pattern, stepwise decrease by one atom
							keep_doing<-TRUE
							keep_w<-TRUE # default
							at_step<-1
							from_int<-min_intensities[1]
							from_mass_min<-min_masses[1]
							from_mass_max<-max_masses[1]							
							while(keep_doing){
								from_mass_min<-(from_mass_min+iso_del_mass)
								from_mass_max<-(from_mass_max+iso_del_mass)							
								from_int<-(
									(from_int)*
									(1/iso_rat_abund)/
									(at_step)*
									(w-at_step+1)
								)
								if(any(
									(masses>=from_mass_min) &
									(masses<=from_mass_max)	&
									(max_intensities>=from_int) 
								)){								
									found_peak<-TRUE
									if(from_int<=int_cut_up){ # reached end = pattern below int_cut
										keep_doing<-FALSE
									}else{
										at_step<-(at_step+1)
									}
								}else{
									if(from_int>int_cut_up){ # pattern still above int_cut?
										keep_w<-FALSE
									}				
									keep_doing<-FALSE
								}
							}
							if(keep_w){
								n_atom_iter<-w
								break;
							}
						}
						max_over_isotopes<-c(max_over_isotopes,n_atom_iter)
						peaks_over_isotopes<-c(peaks_over_isotopes,found_peak)
					}else{ # initial upper bound == 0
						# get estimate from int_cut_up:
						n_atom<-floor(int_cut_up/min_intensities[1]*iso_rat_abund)
						max_over_isotopes<-c(max_over_isotopes,n_atom)					
						peaks_over_isotopes<-c(peaks_over_isotopes,FALSE)				
					}
				} # if mass at M+X found
			} # over all isotopes
			# compare bounds over isotopes; possibly restrict by mass & use_C -> C-peak!
			if(must_peak & (!any(peaks_over_isotopes))){next}
			end_count<-min(max_over_isotopes) # = minimum over maxima
			if(!use_C || (w_C==0)){
				iso_mass<-(as.numeric(isotopes[at_isotopes[1],3])*end_count/charges[n])
				if(iso_mass<=masses[1]){ # count larger than monoisotopic mass?
					bounds[n,i]<-end_count
				}else{ # shrink count to lowest possible mass
					bounds[n,i]<-floor(masses[1]/(as.numeric(isotopes[at_isotopes[1],3])/charges[n]))
				}
			}else{
				iso_mass<-((as.numeric(isotopes[at_isotopes[1],3])+(1/w_C*12))*end_count/charges[n])
				if(iso_mass<=masses[1]){ # count larger than monoisotopic mass?
					bounds[n,i]<-end_count
				}else{ # shrink count to lowest possible mass
					iso_mass<-((as.numeric(isotopes[at_isotopes[1],3])+(1/w_C*12))*1/charges[n])
					bounds[n,i]<-floor(masses[1]/iso_mass)
				}				
			}
		} # over all charges	
	} # over all elements
	############################################################################
	# use_C -> maximum C-number -> maximum atoms for other elements
	if(any(elements=="C") & use_C){
		for(n in 1:length(charges)){
			max_C<-bounds[n,colnames(bounds)=="C"]
			for(i in 1:length(bounds[1,])){	
				if(colnames(bounds)[i]=="C"){next}
				if(bounds[n,i]==0){next}
				at_isotopes<-which(as.character(isotopes[,1])==colnames(bounds)[i])
				w_C<-as.numeric(unique(isotopes[at_isotopes,5]))
				if(w_C==0){next}
				max_count_with_C<-(max_C*w_C)
				if(max_count_with_C<bounds[n,i]){
					bounds[n,i]<-max_count_with_C
				}
			}
		}
	}
	############################################################################
	return(bounds)
	
}




















