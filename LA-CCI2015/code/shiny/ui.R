# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

# Define UI
shinyUI(fluidPage(
	# Readme
	fluidRow(
		column(12,
			 includeHTML("README.html")
		)
	),
	tags$hr(),
	# Inputs
	fluidRow(
		column(4,
			 tags$h3("Main parameters"),
			 tags$h4("(Slide to change)"),
			 tags$hr(),
			 # Number of offspring
			 sliderInput(
			 	inputId = "lambda", 
			 	label   = h5("Number of offspring (lambda)"),
			 	min     = 1,
			 	max     = 20,
			 	step    = 1,
			 	value   = 1
			 ),
			 # Standard deviation of mutation
			 sliderInput(
			 	inputId = "sigma", 
			 	label   = h5("Mutation St. Dev. (sigma)"),
			 	min     = 0.01,
			 	max     = 0.5,
			 	step    = 0.01,
			 	value   = 0.05
			 ),
			 # Problem dimension
			 sliderInput(
			 	inputId = "N", 
			 	label   = h5("Problem dimension (n)"),
			 	min     = 2,
			 	max     = 100,
			 	step    = 1,
			 	value   = 10
			 )),
		#
		column(4,
			 tags$h3("Stop criteria parameters"),
			 tags$h4("(0 = not used)"),
			 tags$hr(),
			 # Normalized distance from optimum
			 sliderInput(
			 	inputId = "Rnorm", 
			 	label   = h5("Standardized dist. optimum (R/sigma)"),
			 	min     = 0.25,
			 	max     = 10,
			 	step    = 0.25,
			 	value   = 1.5
			 ),
			 # Max iterations
			 sliderInput(
			 	inputId = "maxiter", 
			 	label   = h5("Maximum number of iterations"),
			 	min     = 0,
			 	max     = 100000,
			 	step    = 1000,
			 	value   = 0
			 ),
			 # Significance level for statistical stop criterion
			 sliderInput(
			 	inputId = "alpha", 
			 	label   = h5("Significance level (alpha)"),
			 	min     = 0.01,
			 	max     = 0.5,
			 	step    = 0.01,
			 	value   = 0.05
			 )),
		#
		column(4,
			 tags$h3("Info and Actions"),
			 tags$h4("(Updates with sliders)"),
			 tags$hr(),
			 tags$h4(tags$strong(textOutput("mstar"))),
			 tags$h4("Click below to run optimization"),
			 tags$h4("(It may take a while)"),
			 actionButton("run", "Click to RUN")
		),
		#
		column(12,
			 tags$hr()),
		#
		#Output
		column(12, 
			 plotOutput("viewplot",
			 	     width  = "100%",
			 	     height = "600px")
		)
	)
))