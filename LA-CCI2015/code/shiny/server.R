# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

# Load required functions
filenames <- dir(path = "../")
toload    <- filenames[grep(".R", filenames)]
toload    <- toload[-grep("example", toload)]
for (file in toload) source(paste0("../",file))

library(ggplot2)

# Define server logic
shinyServer(function(input, output) {
	# Calculate m*
	output$mstar <- renderText(
		paste("Required iterations without improvement:",
			calc_mstar(alpha  = input$alpha, 
				     Rnorm  = input$Rnorm,
				     n      = input$N,
				     lambda = input$lambda)))
	
	# Run optimization after user clicks RUN 
	observeEvent(input$run,{
		# Optimization run - setup
		probpars <- list(name  = "sphere", 
				     xmin = rep(-0.5, input$N), 
				     xmax = rep(0.5,  input$N))
		stopcrit <- list(names = c("stop_stat", "stop_maxiter"), 
				     Rnorm = input$Rnorm, 
				     sig.level = input$alpha,
				     maxiter = ifelse(input$maxiter==0,
				     		     Inf,
				     		     input$maxiter))
		
		# Run optimization
		out <- run_ES(input$lambda, input$sigma, probpars, stopcrit)
		
		output$viewplot <- renderPlot({
			ggplot(data = data.frame(t = seq_along(out$fhist),
							 f = out$fhist),
				 aes(x = t, y = f)) + 
				geom_line() + geom_point() +
				ggtitle("Algorithm progress") + 
				xlab("Iteration") + ylab ("Function value") +
				theme(axis.title = element_text(size = 16),
					axis.text  = element_text(size = 14),
					plot.title = element_text(size = 18)) + 
				annotate("text", 
					   x = 0.1*out$iter, y = out$fhist[1], 
					   label = paste("Final value:", signif(out$f,6), 
					   		  "\nTotal iters:", out$iter), 
					   hjust = 0) + 
				scale_x_log10()
		})
	})
})