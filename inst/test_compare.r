########################################################################
# compare function #####################################################
cat("\n Testing compare function ...")


# retrieve simulated data, all ordered
for(i in 1:20){

	ncols<-sample(seq(20),1)
	nrows<-sample(seq(1E4),1)
	if(ncols==1){next}
	# example data sets:
	a <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- rbind(a,b) # example of b containing a
	
	ord1<-do.call(order,as.data.frame(a))
	a<-a[ord1,]
	
	ord2<-do.call(order,as.data.frame(b))
	b<-b[ord2,]	
		
	found<-rows_compare(a,b,row_order=FALSE,column_order_a=FALSE,column_order_b=FALSE,get_index=TRUE)
	
	if( any(found==0) ){
		cat("WARNING: issue_1 in function compare.r")
	}else{
		for(j in 1:length(found)){
			if(!all(b[found[j],]==a[j,])){
				cat("\nWARNING: issue_2 in function compare.r")			
			}
		}
	
	}

}


# retrieve simulated data, search matrix a unordered
for(i in 1:20){

	ncols<-sample(seq(20),1)
	nrows<-sample(seq(1E4),1)
	if(ncols==1){next}
	# example data sets:
	a <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- rbind(a,b) # example of b containing a
	
	ord2<-do.call(order,as.data.frame(b))
	b<-b[ord2,]	
		
	found<-rows_compare(a,b,row_order=FALSE,column_order_a=TRUE,column_order_b=FALSE,get_index=TRUE)
	
	if( any(found==0) ){
		cat("WARNING: issue_3 in function compare.r")
	}else{
		for(j in 1:length(found)){
			if(!all(b[found[j],]==a[j,])){
				cat("\nWARNING: issue_4 in function compare.r")			
			}
		}
	
	}

}


# retrieve simulated data, search matrix a and target matrix b unordered
for(i in 1:20){

	ncols<-sample(seq(20),1)
	nrows<-sample(seq(1E4),1)
	if(ncols==1){next}
	# example data sets:
	a <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- matrix(sample.int(1E6,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
	b <- rbind(a,b) # example of b containing a
	b <- b[sample(dim(b)[1],dim(b)[1],replace = F),] 	
		
	found<-rows_compare(a,b,row_order=FALSE,column_order_a=TRUE,column_order_b=TRUE,get_index=TRUE)
	
	if( any(found==0) ){
		cat("WARNING: issue_5 in function compare.r")
	}else{
		for(j in 1:length(found)){
			if(!all(b[found[j],]==a[j,])){
				cat("\nWARNING: issue_6 in function compare.r")			
			}
		}
	
	}

}


cat(" done.")
########################################################################



#The below function `compare` compares data.frames or matrices `a,b` to find row matches of `a` in `b`. It returns the first row position in `b` which matches (after some internal sorting required to speed thinks up). Rows in `a` which have no match in `b` will have a return value of `0`. Should handle numeric, character and factor column types and mixtures thereof (the latter for `data.frames` only). Check the example below the function definition.


    compare<-function(a,b){
    
    	#################################################
    	if(dim(a)[2]!=dim(b)[2]){
    		stop("\n Matrices a and b have different number of columns!")
    	}
    	if(!all(sapply(a, class)==sapply(b, class))){
    		stop("\n Matrices a and b have incomparable column data types!")	
    	}
    	#################################################
    	if(is.data.frame(a)){
    		i <- sapply(a, is.factor)
    		a[i] <- lapply(a[i], as.character)
    	}
    	if(is.data.frame(b)){
    		i <- sapply(b, is.factor)
    		b[i] <- lapply(b[i], as.character)
    	}
    	len1<-dim(a)[1]
    	len2<-dim(b)[1]
    	ord1<-do.call(order,as.data.frame(a))
    	a<-a[ord1,]
    	ord2<-do.call(order,as.data.frame(b))
    	b<-b[ord2,] 	
    	#################################################
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
    	#################################################
    	found<-found[order(ord1)]
    	found<-ord2[found]
    	return(found)
    	
    }
    
    ncols<-10
    nrows<-1E4
    if(ncols==1){next}
    # example data sets:
    a <- matrix(sample(LETTERS,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
    b <- matrix(sample(LETTERS,size = (ncols*nrows), replace = T), ncol = ncols, nrow = nrows)
    b <- rbind(a,b) # example of b containing a
    b <- b[sample(dim(b)[1],dim(b)[1],replace = F),] 
    found<-compare(a,b)
    
    a<-as.data.frame(a) # = conversion to factors
    b<-as.data.frame(b) # = conversion to factors
    found<-compare(a,b)




