# Examples

# Clear workspace
rm(list=ls())

# Load all relevant files
filenames <- dir()
toload    <- filenames[grep(".R", filenames)]
toload    <- toload[-grep("example", toload)]
for (file in toload) source(file)

# Example 1
lambda   <- 1
sigma    <- 0.05
probpars <- list(name  = "sphere", 
			 xmin = rep(-2, 10), 
			 xmax = rep(1, 10))
stopcrit <- list(names = "stop_stat", Rnorm = 1, sig.level = 0.05)
seed     <- 99
out 	   <- run_ES(lambda, sigma, probpars, stopcrit, seed)


# Plot progress
library(ggplot2)
p <- ggplot(data = data.frame(t = seq_along(out$fhist),
					f = out$fhist),
		aes(x = t, y = f)) + 
	geom_line() + geom_point() +
	ggtitle("Algorithm progress") + 
	xlab("Iteration") + ylab ("Function value") +
	theme(axis.title = element_text(size = 16),
		axis.text  = element_text(size = 14),
		plot.title = element_text(size = 18))

p + annotate("text", 
		 x = 0.8*out$iter, y = out$fhist[1], 
		 label = paste("Final value:", signif(out$f,6), 
		 		  "\nTotal iters:", out$iter), 
		 hjust = 0)
ggsave(filename = "../figures/Example1-lin-lin.pdf",
	 width = 24,
	 height = 16)

p + annotate("text", 
		 x = 0.1*out$iter, y = out$fhist[1], 
		 label = paste("Final value:", signif(out$f,6), 
		 		  "\nTotal iters:", out$iter), 
		 hjust = 0) + 
	scale_x_log10()
ggsave(filename = "../figures/Example1-log-lin.pdf",
	 width = 24,
	 height = 16)
