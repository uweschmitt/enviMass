
library(shiny)
library(plotly)

ui <- fluidPage(
	plotly:::plotlyOutput("plotlii", width = "100%", height = "700px")
)

server <- function(input, output) {
	output$plotlii <- renderPlotly({
		# the data:		
		obs1 <- matrix(ncol=3,
			c(0,2,3,1,2,2,6,6,4)
		)
		df1 <- setNames(data.frame(obs1), c("a", "b", "c"))
		df1$c<-(.5*df1$c) # half values to make them the centre point of the bars
		obs2 <- matrix(ncol=3,
			c(14,16,10,11,12,12,23,23,22)
		)
		df2 <- setNames(data.frame(obs2), c("a", "b", "c"))
		df2$c<-(.5*df2$c) # half values to make them the centre point of the bars
		
		use_lines<-data.frame(
			cbind(
				c(1,2,2,1),
				c(1,1,2,1),
				c(-1,-1,-1,-1)
			)
		)
		use_lines<- setNames(data.frame(use_lines), c("x", "y", "z"))
		
		# the plot:
		p <- plot_ly(type="scatter3d",mode="markers",showlegend=FALSE)%>%
		add_trace(p,
			x = ~a, y = ~b, z = ~c,
			data = df1[df1$a<2,],
			color=I("red"),
			size = I(1),
			name = "Group a",
			error_z=list(
				color="red",
				thickness=0,
				symmetric = TRUE, 
				type = "data" ,
				array = df1[df1$a<2,]$c
			)
		)%>%
		add_trace(p,
			x = ~a, y = ~b, z = ~c,
			data = df1[df1$a>=2,],
			color=I("gray"),
			size = I(1),
			name = "Group a",
			error_z=list(
				color="gray",
				symmetric = TRUE, 
				type = "data" ,
				array = df1[df1$a>=2,]$c
			)
		)%>%	
		add_lines(p, 
			x = ~x, y = ~y, z = ~z,
			data = use_lines,
			color=I("blue")
		)
	})

	
		
	
}

shinyApp(ui, server)















ui <- fluidPage(
HTML('<p><font>
	<span style="color:black"><b>&bull;</b></span> Picked peaks <span style="color:gray"><b>&bull;</b></span> Raw data centroids <span style="color:red"><b>&bull;</b></span> Peak centroids 							
</p></font>')
)
server <- function(input, output) {
}

shinyApp(ui, server)






		load(file.path(logfile[[1]],"MSlist","21"), envir=as.environment(".GlobalEnv"))
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			