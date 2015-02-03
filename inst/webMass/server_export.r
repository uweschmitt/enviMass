##############################################################################
# Export peaklist selection ##################################################
##############################################################################
exported_1<-reactive({
    input$expo_profiles

    if( (isolate(init$a)=="TRUE") & isolate(input$expo_profiles) ){




	
		return("Export finished")	
	}else{


	
	

		return("Export failed")
	}
})	
output$expo1<-renderText(paste(exported_1()))
##############################################################################


