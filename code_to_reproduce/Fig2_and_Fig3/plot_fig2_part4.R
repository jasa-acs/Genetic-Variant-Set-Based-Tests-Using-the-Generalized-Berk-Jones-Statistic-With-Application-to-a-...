# The purpose of this script is to plot the figure in the right panel of Figure 3
# of the submitted manuscript.

# We read in the output produced by gen_sim_data_fig2.R and gen_omni_data_fig2.R,
# calculate the omnibus statistic p-values, and set up the plotting window with various 
# options. 
# For example, we change options like the title of the plot as well as the colors and textures
# of the lines.
# Output is the image file d100_3_3_3.eps.

# Make sure all the jobs from gen_sim_data_fig2.R and gen_omni_data_fig2.R have finished running,
# otherwise the necessary data will not exist and there will be errors.

########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'

# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

########################################################################################
# No manual changes need to be made after this point
########################################################################################

########################################################################################
# Load library
library(mvtnorm, lib=libs)

# param_row is equal to 1,2,6, or 8. These are the four settings we will plot.
param_row <- 8

# First build the new parameter matrix
rho1_vec <- c(0, 0.3)
rho2_vec <- c(0, 0.3)
rho3_vec <- c(0, 0.3)
d_vec <- c(100)
param_mat <- expand.grid(rho1_vec, rho3_vec, rho2_vec, d_vec)
colnames(param_mat) <- c('rho1', 'rho2', 'rho3', 'd')
########################################################################################


########################################################################################
# Read in simulation results
setwd( data_dir )
param_results <- read.table( paste('power_results_S', param_row, '_1.txt', sep=''), sep='\t', header=T)
for (i in 2:50) {
	temp_fname <- paste('power_results_S', param_row, '_', i, '.txt', sep='')
	temp_file <- tryCatch(read.table(temp_fname, header=T, sep='\t'), error=function(e) e, warning=function(w) w)
	if (length(class(temp_file)) > 1) {
		next
		cat('Missing S ', Snum, ' aID ', aID, '\n')
	}
	
	param_results <- rbind(param_results, temp_file)
}
########################################################################################


########################################################################################	
# Now summarize the results to get power for each ncausal
alpha <- 0.01
power_mat <- matrix(data=NA, nrow=10, ncol=10)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('rho1', 'rho2', 'rho3', 'd', 
							'ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni')
for (i in 1:10)
{
	# Cut down to 500 results for each run
	temp_data <- param_results[which(param_results$ncausal==i),]
	cat(i, '\n')
	
	# Report parameters
	power_mat$rho1[i] <- temp_data$rho1[1]
	power_mat$rho2[i] <- temp_data$rho2[1]
	power_mat$rho3[i] <- temp_data$rho3[1]
	power_mat$d[i] <- temp_data$d[1]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	
	# Read in the omnibus correlation matrix file
	d <- param_mat$d[param_row]
	rho1 <- param_mat$rho1[param_row]
	rho2 <- param_mat$rho2[param_row]
	rho3 <- param_mat$rho3[param_row]
	temp_omni_fname <- paste('d', d, '_rho1_', 10*rho1, '_rho2_', 10*rho2, '_rho3_', 
				10*rho3, '_ncausal_', i, '.txt', sep='')
	temp_omni_mat <- read.table(temp_omni_fname, header=T)
		
	# Get the omni test statistics
	omni_pvec <- rep(NA, nrow(temp_data))
	for (j in 1:nrow(temp_data)) {
		omni_stat <- temp_data $omni_stat[j] 
		omni_bound <- tryCatch(qnorm(1-omni_stat), warning=function(w) w, error=function(e) e)
		omni_Z_bounds <- rep(omni_bound, 4)
		omni_pvec[j] <- 1 - pmvnorm(lower=rep(-Inf, 4), upper=omni_Z_bounds,
							sigma=unname(as.matrix(temp_omni_mat)))[1]
	}

	# Record power
	power_mat$GBJ[i] <- sum(temp_data$GBJ_p<alpha) / nrow(temp_data)
	power_mat$GHC[i] <- sum(temp_data$GHC_p<alpha) / nrow(temp_data)
	power_mat$minP[i] <- sum(temp_data$minP_pvalue<alpha) / nrow(temp_data)
	power_mat$skat[i] <- sum(temp_data$skat_p<alpha) / nrow(temp_data)
	power_mat$omni[i]<- sum(omni_pvec<alpha) / nrow(temp_data)
}
power_mat


##################################################################
# Done reading data.
##################################################################



##################################################################
# Make the plot for simulation 8 (corresponds to right panel of figure 3)
xlab_tex = expression('Number of causal SNPs')
main_tex = expression(paste('100 SNPs, ' , rho[1], '=0.3, ', rho[2], '=0.3, ', rho[3], '=0.3'))

postscript("d100_3_3_3.eps", width=720, height=500)

plot_data <- data.frame(power_mat[1:10,5:10])
par(mar=c(5,6,6,2))
xmax <- 10
plot(1, type='n', xlim=c(1,xmax), ylim=c(0,1), 
		main=main_tex, cex.main=3, xaxt='n', yaxt='n', xlab='', ylab='')
axis(side=1, at=seq(from=1, to=xmax, by=2), labels=as.character(seq(from=1, to=xmax, by=2)), cex.axis=2.5)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=as.character(seq(from=0, to=1, by=0.2)), cex.axis=2.5)
title(ylab='Power', line=4, cex.lab=2.7)
title(xlab=xlab_tex, line=4, cex.lab=2.7)
ss_GBJ <- smooth.spline(plot_data$GBJ ~ plot_data$ncausal, spar=0.4)
ss_GHC <- smooth.spline(plot_data$GHC ~ plot_data$ncausal, spar=0.4)
ss_minP <- smooth.spline(plot_data$minP ~ plot_data$ncausal, spar=0.4)
ss_skat <- smooth.spline(plot_data$skat ~ plot_data$ncausal, spar=0.4)
ss_omni <- smooth.spline(plot_data$omni ~ plot_data$ncausal, spar=0.4)
lines(ss_GBJ, col=2, lwd=5, lty=1, cex=2)
lines(ss_GHC, col=3, lwd=5, lty=2, cex=2)
lines(ss_minP, col=4, lwd=5, lty=3, cex=2)
lines(ss_skat, col=1, lwd=5, lty=4, cex=2)
lines(ss_omni, col=6, lwd=5, lty=5, cex=2)
legend(x=7.5, y=0.45, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
		lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')
		
dev.off()
