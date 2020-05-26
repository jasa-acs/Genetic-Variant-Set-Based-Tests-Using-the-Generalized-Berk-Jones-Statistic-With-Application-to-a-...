# This script reads the bounds files produced by gen_data_fig1.R, and then it
# sets up the parameters necessary to produce the plots in Supplementary Figure 3 of the 
# submitted manuscript.

# Make sure that all the jobs from gen_data_fig1.R have finished running and have
# written all their output to disk before attempting to run this script (will not 
# work if the data does not exist to be read).

# Output will be four files with names like d20_alpha1_rho3_pct50_pscale.eps. These correspond
# to each of the four panels of Supplementary Figure 3.



###################################################################
# Change data_dir and libs to the appropriate directory (make sure you have all the necessary data files created in previous step),
# then run the rest of the code to produce the plots as in Figure 1 of the manuscript.

# Where are the tables of data saved from gen_data_fig1.R?
data_dir <- '/users/user/desktop'
setwd(data_dir)

########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'


###################################################################
# No manual changes need to be made after this point
###################################################################
library(dplyr, lib=libs)

###################################################################
# Define the parameter matrix so we know which run corresponds to which parameters.
runs_to_use <- c(3, 5, 19, 21)
d_vec <- c(20, 50, 100, 200)
alpha_vec <- c(0.01, 0.05)
rho_vec <- c(0.3)
percent_cor <- c(0.25, 0.5, 0.75, 1)
param_mat <- expand.grid(alpha_vec, rho_vec, percent_cor, d_vec)
param_mat <- data.frame(param_mat)
colnames(param_mat) <- c('alpha', 'rho', 'percent', 'd')

###################################################################
# Loop over the four plots to reproduce
for (i in 1:length(runs_to_use)) {
	param_row <- runs_to_use[i]
	xaxis_length <- 19
	
	# Get the parameters for that run
	d <- param_mat$d[param_row]
	alpha <- param_mat$alpha[param_row]
	rho <- param_mat$rho[param_row]
	percent_cor <- param_mat$percent[param_row]
	
	# Open the correct file
	temp_fname <- paste('d', d, '_alpha', alpha*100, '_rho', rho*10, 
			'_pct', percent_cor*100, '_bounds.txt', sep='')
	bounds_data <- read.table(temp_fname, header=T)
	
	# go to p-value scale
	bounds_data  <- bounds_data %>% mutate(BJ = -log10(1-pchisq(BJ^2, df=1))) %>%
	    mutate(GBJ = -log10(1-pchisq(GBJ^2, df=1))) %>%
	    mutate(GHC = -log10(1-pchisq(GHC^2, df=1))) %>%
	    mutate(HC = -log10(1-pchisq(HC^2, df=1)))

	# Now we plot it
	temp_plotname <- paste('d', d, '_alpha', alpha*100, '_rho', rho*10, 
	                       '_pct', percent_cor*100, 'pscale.eps', sep='')
	postscript(temp_plotname, width=720, height=500)
	
	# Some options on the axes, etc.
	percent_cor <- 100*percent_cor
	plot_data <- bounds_data
	par(mar=c(6,6,6,2))
	xlab_tex = expression(paste('Ordered p-values, ', group("",p,"")[(j)]))
	main_tex = bquote(paste(.(d), ' SNPs, ', .(percent_cor), '% correlated, ',
		 alpha, '=', .(alpha), sep=''))
	ymax <- ceiling(plot_data$BJ[d])
	ymin <- floor(min(plot_data[d-49,]))
	xmax <- max(plot_data$Index)
	xmin <- max(plot_data$Index)-xaxis_length
	plot(1, type='n', xlim=c(xmin, xmax), ylim=c(ymin, ymax), 
		main=main_tex, cex.main=3, xaxt="n", yaxt="n", xlab='', ylab='')
	axis(side=1, at=seq(from=xmin+1, to=xmax, by=2), 
	     labels=as.character(seq(from=20, to=2,by=-2)), cex.axis=2.5)
	axis(side=2, at=seq(from=ymin, to=ymax, by=1), 
	     labels=as.character(seq(from=ymin, to=ymax, by=1)), cex.axis=2.5)
	title(ylab=expression(paste(-log[10], "(p-value)")), line=3.5, cex.lab=2.7)
	title(xlab=xlab_tex, line=5, cex.lab=2.7)
	lines(cbind(plot_data$Index, plot_data$BJ), type='o', col=4, lwd=3, cex=2, lty=1, pch=19)
	lines(cbind(plot_data$Index, plot_data$GBJ), type='o', col=4, lwd=3, cex=2, lty=2, pch=19)
	lines(cbind(plot_data$Index, plot_data$HC), type='o', col=2, lwd=3, cex=2, lty=1, pch=18)
	lines(cbind(plot_data$Index, plot_data$GHC), type='o', col=2, lwd=3, cex=2, lty=2, pch=18)
	legend('topleft', c('BJ', 'GBJ', 'HC', 'GHC'), lty=c(1,2,1,2), lwd=c(3,3,3,3),
		col=c(4,4,2,2), pch=c(19, 19, 18, 18), bty='n', cex=2)
	
	# Save
	dev.off()
}

