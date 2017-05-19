#' @title Compare matrix rows
#'
#' @description Compare if rows in one matrix can be found in another
#'
#' @param a Matrix with n columns and x rows
#' @param b Matrix with n columns and y rows
#' @param row_order Logical. If TRUE, the column order is irrelevant; otherwise it is. Does a sort in each row.
#' @param column_order Logical. Default TRUE. 
#' @param get_index Logical. If TRUE, the index of the matching row of a in b is returned. This is 0 if no match is found and refers to the first matach if duplcaies are present in b. 
#' 
#' @details Rows from matrix a are checked for their presence in matrix b. This is done column-wise if row_order=FALSE, i.e., the first match is sought in the first column, the second in the second, etc. Argument column_order ensures that rows are sorted so that columns contain only increasing numbers - a must for the algorithm to work. Can be skipped if matrices have been sorted accordingly beforehand.
#' 

rows_compare<-function(a,b,row_order=FALSE,column_order_a=TRUE,column_order_b=TRUE,get_index=FALSE){

	#################################################
	if(dim(a)[2]!=dim(b)[2]){
		stop("\n Matrices a and b have different number of columns?")
	}
	#################################################
    len1<-dim(a)[1]
    if(row_order){
        a<-t(apply(t(a), 2, sort))
        b<-t(apply(t(b), 2, sort))
    }
    if(column_order_a){	
		ord1<-do.call(order,as.data.frame(a))
		a<-a[ord1,]
	}
	len2<-dim(b)[1]
    if(column_order_b){	
		ord2<-do.call(order,as.data.frame(b))
		b<-b[ord2,] 
	}	
	#################################################
	if(FALSE){ # R version of below C code
		found<-rep(0,len1)  
		dims<-dim(a)[2]
		do_dims<-c(1:dim(a)[2])	
		at<-1
		for(i in 1:len1){
			for(m in do_dims){
				while(b[at,m]<a[i,m]){
					at<-(at+1)      
					if(at>len2){break}              
				}
				if(at>len2){break}
				if(b[at,m]>a[i,m]){break}
				if(m==dims){found[i]<-at}
			}
			if(at>len2){break}
		}
	}
	#################################################	
	if(TRUE){
		results<-rep(0,len1)
		found <- .Call("compare",
			as.matrix(a),
			as.matrix(b),	
			as.matrix(results),				
			PACKAGE="enviMass" 
		)
	}
	#################################################
	if(!get_index){
		if(!column_order_a & !column_order_b){	
			return(found>0)		
		}	
		if(column_order_a & !column_order_b){	
			found<-found[order(ord1)]
			return(found>0)
		}	
		if(column_order_a & column_order_b){	
			found<-found[order(ord1)]
			found<-order(ord2)[found]
			return(found>0)
		}
	}else{
		if(!column_order_a & !column_order_b){	
			return(found)		
		}	
		if(column_order_a & !column_order_b){	
			found<-found[order(ord1)]
			return(found)
		}	
		if(column_order_a & column_order_b){	
			found<-found[order(ord1)]
			found<-ord2[found]
			return(found)
		}		
	}
}




