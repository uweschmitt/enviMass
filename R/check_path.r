#' @title Check enviMass folder path
#'
#' @description \code{check_path} checks if a path to an enviMass folder is valid and accessible / writable.
#'
#' @param project_path String of path
#' 
#' @details enviMass workflow function
#' 

project_path<-"F:/MS/Projects/Package_enviMass"

check_path<-function(project_path){
	
	say_path<-"Project path ok"
	##################################################################
	skip_rest<-FALSE
	if(!file.exists(project_path)){
		say_path<-paste("File or folder ",project_path," does not exist!")
		skip_rest<-TRUE
	}
	if(file.access(project_path, mode = 0)!=0 & !skip_rest){ # slightly different than above
		say_path<-paste("File or folder ",project_path," does not exist!")
	}
	if(file.access(project_path, mode = 2)!=0 & !skip_rest){
		say_path<-paste("Not allowed to write into ",project_path,". Revise user permissions?")
	}	
	if(file.access(project_path, mode = 4)!=0 & !skip_rest){
		say_path<-paste("Not allowed to read from ",project_path,". Revise user permissions?")
	}	
	# try to write into folder #######################################
	if(FALSE){ # windows fuck yourself
		a<-try({
			dir.create(
				file.path(project_path,"delete_this_trial_folder")
			)
		},silent=TRUE)
		if((a==FALSE) || (class(a)=="try-error" )){
			say_path<-paste("Cannot create ",project_path,". Revise user permissions, path validity, ... ?")	
		}else{
			file.remove(file.path(project_path,"delete_this_trial_folder"))
		}
	}
	##################################################################
	return(say_path);
	
}

