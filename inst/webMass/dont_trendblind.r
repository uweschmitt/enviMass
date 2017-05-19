

			path=file.path(logfile[[1]],"pics","boxprofile_pos")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr4p<-list(src=path)
				output$boxprofile<-renderImage(expr4p, deleteFile = FALSE)
			path=file.path(logfile[[1]],"pics","boxprofile_neg")
				png(filename = path, bg = "white")
				plot.new();plot.window(xlim=c(0,1),ylim=c(0,1));text(0.5,0.5,"nothing selected \n or not available",cex=1)
				dev.off()
				expr4n<-list(src=path)
				output$boxprofile<-renderImage(expr4n, deleteFile = FALSE)			
