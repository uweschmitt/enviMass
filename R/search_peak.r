#' @title Match with peaklist
#'
#' @export
#'
#' @description \code{search_peak} searchs for a peak with a certain mass and RT in a peaklist.
#'
#' @param peaklist Matrix or data.frame of peaks. 3 columns: m/z, intensity, RT.
#' @param mz Numeric vector. Mass(es) to be matched.
#' @param dmz Numeric. +/- m/z tolerance (precision)
#' @param ppm Logical. \code{dmass} given in ppm?
#' @param RT Numeric vector. Retention times to be matched; units equal to those of the input files
#' @param dRT Numeric. RT tolerance window; units equal to those of the input files
#' @param onlymax Logical. Return only the matched peak of highest intensity. 
#' @param intratio. Logical FALSE or numeric. Only matching peaks >= this intensity ratio are returned.  
#' @param int. Logical FALSE or numeric vector. Intensities of the peaks to be matched; required when using intratio.
#' @param get_matches. Logical. Return all matches (TRUE, default) or only indicate if any have been found.
#'
#' @return Vector of characters. For each code{mz,RT} input, the vector entry specifies the match indices in \code{peaklist}.
#' 
#' @details  enviMass workflow function
#' 
#' @seealso \code{screening}



search_peak<-function(peaklist,mz,dmz=5,ppm=TRUE,RT,dRT,onlymax=FALSE,int_ratio=FALSE,int=FALSE,get_matches=TRUE){

  ##############################################################################
  # redefine order of peaklist & mz, keep ord1 and ord2 to resort
  peaks<-peaklist[order(peaklist[,1],decreasing=TRUE),];
  ord1<-order(peaklist[,1],decreasing=TRUE);
  mass<-mz[order(mz,decreasing=TRUE)];
  ord2<-order(mz,decreasing=TRUE);
  ##############################################################################
  if(length(dmz)==1){dmz<-rep(dmz,length(mz))}
  if(length(dRT)==1){dRT<-rep(dRT,length(mz))}
  result<-rep("FALSE",length(mz));
  leng<-length(peaks[,1]);
  k<-c(1);
  options(digits=10)
  ##############################################################################
  for(i in 1:length(mass)){
    deletes<-c();
    # calculate bounds per target mass #########################################
    if(ppm==TRUE){
    	target_low<-as.numeric(mass[i])-(as.numeric(dmz[i])*(as.numeric(mass[i])/1e6)*2);
    	target_up<-as.numeric(mass[i])+(as.numeric(dmz[i])*(as.numeric(mass[i])/1e6)*2);
    }else{
    	target_low<-as.numeric(mass[i])-((as.numeric(dmz[i])/1000)*2);
    	target_up<-as.numeric(mass[i])+((as.numeric(dmz[i])/1000)*2);
    }
    # reset k ##################################################################
    # decrease k ...    
    while((k-1)>0 ){
      if(target_up >= as.numeric(peaks[k-1,1])){    
        k<-c(k-1);
      }else{
        break;
      }
    };    
    # or increase k:
    while( (k+1)<=as.numeric(leng) ){
      if(target_up < as.numeric(peaks[k+1,1])){
        k<-c(k+1);
      }else{
        break
      }
    };
   	n<-k;   	
   	# increase n:
   	while( (n+1)<=as.numeric(leng) ){
   	  if(target_low <= as.numeric(peaks[n+1,1])){
		    n<-c(n+1);
		  }else{
		    break;
		  }
    };
   	for(f in k:n){
    	if(   (as.numeric(peaks[f,1]) >= as.numeric(target_low)) &&
    		    ( as.numeric(peaks[f,1]) <= as.numeric(target_up))  &&
    		    ( as.numeric(peaks[f,3]) >= (as.numeric(RT[ord2[i]])-as.numeric(dRT[ord2[i]]))) &&
    		    ( as.numeric(peaks[f,3]) <= (as.numeric(RT[ord2[i]])+as.numeric(dRT[ord2[i]])))){
          deletes<-c(deletes,f);
        }
    }
    # save results #############################################################
    if(length(deletes)>0){
		if(is.numeric(int_ratio)){
			deletes<-deletes[
				int_ratio>=(int[ord2[i]]/peaklist[ord1[deletes],2])
			]
		}
	}
	if(length(deletes)>0){
		if(!get_matches){
			result[ord2[i]]<-"TRUE"
			next;
		}
		if(onlymax){ # only save matched peak with highest intensity
			if(length(deletes)==1){
				result[ord2[i]]<-as.character(ord1[deletes]);
			}else{
				result[ord2[i]]<-as.character(
					(ord1[deletes])[which.max(peaklist[ord1[deletes],2])]
				)
			}
		}else{
			result[ord2[i]]<-as.character(ord1[deletes[1]]);
			if(length(deletes)>1){
				for(j in 2:length(deletes)){
					result[ord2[i]]<-paste(result[ord2[i]],"/",as.character(ord1[deletes[j]]))
				}
			}
		}	
	}
  }
  ##############################################################################
  return(result);

}


