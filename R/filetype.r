#' @title Checks file for consistency with .raw or .mzXML files
#'
#' @description Checks file for consistency with .raw or .mzXML files, returns file extension.
#'
#' @param filename Character string. Filename.
#' @param check Logical. Only return whether extension matches .mzXML, .raw, .RAW or .Raw?
#' 
#' @details enviMass workflow function
#' 

filetype<-function(filename,check){
	
	################################################
	if(!is.character(filename)){stop("filename not a character string.")}
	if(!is.logical(check)){stop("check not logical")}
	################################################
	done<-FALSE
	if(nchar(filename)>=4){
		if(
			substr(filename,nchar(filename)-3,nchar(filename))==".RAW" ||
			substr(filename,nchar(filename)-3,nchar(filename))==".raw" ||
			substr(filename,nchar(filename)-3,nchar(filename))==".Raw"		
		){
			f_type<-".RAW"
			done<-TRUE
		}
	}
	if(nchar(filename)>=6 & !done){
		if(
			substr(filename,nchar(filename)-5,nchar(filename))==".mzXML" 	
		){
			f_type<-".mzXML"
			done<-TRUE
		}
	}
	################################################
	if(check){
		return(done)
	}else{
		return(f_type)
	}
}
