	# busy indicator
	HTML("
			<head>
			<script type='text/javascript'>
				setInterval(function(){
				  if( ($('html').attr('class')=='shiny-busy') ){
					setTimeout(function() {
					if ($('html').attr('class')=='shiny-busy') {
						if($('#textit').html()!='Waiting...' ){
							$('div.busy1').show()
						}
						if($('#textit').html()=='Waiting...'){
							$('div.busy2').show()
						}
					}	
					},1000)	
				  } else {
					$('div.busy1').hide()
					$('div.busy2').hide()
				  }
				},100)	
			</script>
			</head>
			<body>
			<div class='busy1' style='
				position:fixed;
				top: 60%;
				left: 25%;
				margin-top: -100px;
				margin-left: -10px;
				background: rgba(260,260,260,1);
				text-align: center;
				padding-top: 20px;
				padding-left: 30px;
				padding-bottom: 40px;
				padding-right: 30px;
				border-radius: 6px;
				border: solid;
				border-width: 4px;
				font-size:1.6em;
				z-index:30;
				color: green'>
				<p> Calculation in progress ... please wait</p>
				<div id='logobusy1' class='shiny-image-output'></div>	
			</div>
			<div class='busy2' style='
				position:fixed;
				top: 60%;
				left: 25%;
				margin-top: -100px;
				margin-left: -10px;
				background: rgba(260,260,260,1);
				text-align: center;
				padding-top: 20px;
				padding-left: 30px;
				padding-bottom: 40px;
				padding-right: 30px;
				border-radius: 6px;
				border: solid;
				border-width: 4px;
				font-size:1.6em;
				z-index:30;
				color: grey'>
				<p> Loading project ... please wait</p>
				<div id='logobusy2' class='shiny-image-output'></div>	
			</div>
			</body>
			<font color='white'>
	")