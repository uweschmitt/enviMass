#' @title Thermo .raw file conversion and centroidization with ProteoWizard (PW)
#'
#' @export
#'
#' @description \code{PWfile} calls PW msconvert
#'
#' @param infile Path to input file
#' @param folderout Path to output folder
#' @param msconvert_path Path to PW msconvert executable (including \\msconvert).
#' @param notintern Ignore
#' @param use_format Output format
#' 
#' @details  enviMass workflow function. PW (not just msconvert) needs to be installed seperately, including the vendor library.
#' 


PWfile <-
function(infile,folderout,msconvert_path,notintern=FALSE,use_format="mzXML"){

      ##########################################################################
      # checks & setups ########################################################
      if(nchar(Sys.which("msconvert")[[1]])==0){
        cat("msconvert not in system path - ok if msconvert_path correct")
      }
      if(
          sum(substr(infile,nchar(infile)-3,nchar(infile))!=".RAW",substr(infile,nchar(infile)-3,nchar(infile))!=".raw")==1
      ){cat("running .RAW file conversion.")}	  
      ##########################################################################
      # convert ################################################################
      there2<-paste(" -o ",shQuote(folderout),sep="")
	  filtered0<-paste(shQuote("--"),use_format,sep="")
	  filtered1<-paste(shQuote("--32"),sep="")
	  filtered2<-paste(shQuote("--zlib"),sep="")
      filtered3<-paste(" --filter ",shQuote("peakPicking true 1-2"),sep="")
      filtered4<-paste(" --filter ",shQuote("msLevel 1"),sep="")
      system(
              paste(
                shQuote(msconvert_path),
                shQuote(infile),
				filtered1,
				filtered2,
				filtered0,
                filtered3,
                filtered4,
                there2
              )
      ,intern=notintern)
      ##########################################################################

}
