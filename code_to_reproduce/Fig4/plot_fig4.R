# The purpose of this script is to plot the two figures in Figure 4 of
# the submitted manuscript.

# We read in the output produced by gen_data_fig4.R, use splines to generate smoothed power 
# curves, and set up the plotting windows with various options.
# For example, we change options like the title of the plot as well as the colors and textures
# of the lines.
# Outputs are the image files S4_d40_chr5_all.eps and S4_d40_chr5_lowrho.eps.

# Make sure all the jobs from gen_data_fig4.R have finished running,
# otherwise the necessary data will not exist and there will be errors.


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'
library(data.table, lib=libs)

# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

###################################################################
# No manual changes need to be made after this point
###################################################################


########################################################################################
# Read in all the data files from disk
setwd(data_dir)
results <- NULL
alpha <- 0.00001
for (i in 1:320)
{
	# There should be 320 files.
	temp_fname <- paste('chr5_power_S4_', i, '.txt', sep='')
	temp_file <- tryCatch(read.table(temp_fname, header=T, sep='\t'), warning=function(w) w, 
						error=function(e) e)
	
	if (class(temp_file)[1] == 'simpleWarning' | class(temp_file)[1] == 'simpleError') {
		cat("Missing ", i, '\n')
		next
	}
	
 	# Sometimes GBJ will return NA due to lack of numerical precision in R, but we can
 	# tell if the p-value would be low enough from previous runs which did not run into 
 	# numerical precision issues. Corresponds to rerunning and getting an upper bound on the p-value.
	all_genes <- unique(temp_file$start_row)
	GBJ_fixed <- 0
	for (j in 1:length(all_genes)) {
	    temp_gene <- all_genes[j]
	    gene_data <- temp_file[which(temp_file$start_row==temp_gene), ]
	    good_GBJ_data <- gene_data[which(gene_data$GBJ_err == 0), ]
	    GBJ_threshold <- good_GBJ_data$GBJ[which(good_GBJ_data$GBJ_p < alpha)]
	    
	    # If we have good runs from before that did not run into numerical precision errors.
	    if (length(GBJ_threshold) > 0) {
	        GBJ_threshold <- min(GBJ_threshold)
	        fixable_rows <- which(gene_data$GBJ_err != 0 & gene_data$GBJ > GBJ_threshold)
	        if (length(fixable_rows) > 0) {
	            gene_data$GBJ_p[fixable_rows] <- alpha * 0.1
	            temp_file[which(temp_file$start_row==temp_gene), ] <- gene_data
	            GBJ_fixed <- GBJ_fixed + length(fixable_rows)
	        }
	    }
	}
		
	# Sometimes GHC will return NA due to lack of numerical precision in R, but we can
 	# tell if the p-value would be low enough from previous runs which did not run into 
 	# numerical precision issues. Corresponds to rerunning and getting an upper bound on the p-value.
	GHC_fixed <- 0
	for (j in 1:length(all_genes)) {
	    temp_gene <- all_genes[j]
	    gene_data <- temp_file[which(temp_file$start_row==temp_gene), ]
	    good_GHC_data <- gene_data[which(gene_data$GHC_err == 0), ]
	    GHC_threshold <- good_GHC_data$GHC[which(good_GHC_data$GHC_p < alpha)]
	    
	    	   # If we have good runs from before that did not run into numerical precision errors.
	    if (length(GHC_threshold) > 0) {
	        GHC_threshold <- min(GHC_threshold)
	        fixable_rows <- which(gene_data$GHC_err != 0 & gene_data$GHC > GHC_threshold)
	        if (length(fixable_rows) > 0) {
	            gene_data$GHC_p[fixable_rows] <- alpha * 0.1
	            temp_file[which(temp_file$start_row==temp_gene), ] <- gene_data
	            GHC_fixed <- GHC_fixed + length(fixable_rows)
	        }
	    }
	}
	
	# Attach to results
	if (is.null(results)) {
	    results <- temp_file
	} else {
	    mylist <- list(results, temp_file)
	    results <- rbindlist(mylist, use.names=F)
	}
}



################################################################################
# First plot the left panel of Figure 4, all results
################################################################################

temp_results <- results
power_mat <- matrix(data=NA, nrow=10, ncol=7)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni', 'BJ')
for (i in 1:8)
{
	# Cut down to 2000 results for each run
	temp_data <- temp_results[which(temp_results$ncausal==i),]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	# Power
	power_mat$GBJ[i] <- sum(temp_data$GBJ_p<alpha) / nrow(temp_data)
	power_mat$GHC[i] <- sum(temp_data$GHC_p<alpha) / nrow(temp_data)
	power_mat$minP[i] <- sum(temp_data$minP_pvalue<alpha) / nrow(temp_data)
	power_mat$skat[i] <- sum(temp_data$skat_p<alpha) / nrow(temp_data)
	power_mat$omni[i]<- sum(temp_data$omni_p<alpha) / nrow(temp_data)
	power_mat$BJ[i] <- sum(temp_data$BJ_p<alpha) / nrow(temp_data)
}

# Plot
plot_data <- data.frame(power_mat[1:8,1:7])
xlab_tex = expression('Number of causal SNPs')
main_tex = expression(paste('40 SNPs on Chr5, All Data', sep=''))

postscript("S4_d40_chr5_all.eps", width=720, height=500)
par(mar=c(5,6,6,2))
xmax <- 8
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
ss_BJ <- smooth.spline(plot_data$BJ ~ plot_data$ncausal, spar=0.4)
lines(ss_GBJ, col=2, lwd=5, lty=1, cex=2)
lines(ss_GHC, col=3, lwd=5, lty=2, cex=2)
lines(ss_minP, col=4, lwd=5, lty=3, cex=2)
lines(ss_skat, col=1, lwd=5, lty=4, cex=2)
lines(ss_omni, col=6, lwd=5, lty=5, cex=2)
legend(x=5, y=0.5, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
		lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')	
dev.off()






################################################################################
# Next plot the right panel of Figure 4
################################################################################

# Cut to low rho3
low_rho3 <- results[which(results$median_rho3_abs < 0.1), ]
temp_results <- low_rho3
power_mat <- matrix(data=NA, nrow=10, ncol=7)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni', 'BJ')
for (i in 1:8)
{
	# Cut down to 2000 results for each run
	temp_data <- temp_results[which(temp_results$ncausal==i),]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	# Power
	power_mat$GBJ[i] <- sum(temp_data$GBJ_p<alpha) / nrow(temp_data)
	power_mat$GHC[i] <- sum(temp_data$GHC_p<alpha) / nrow(temp_data)
	power_mat$minP[i] <- sum(temp_data$minP_pvalue<alpha) / nrow(temp_data)
	power_mat$skat[i] <- sum(temp_data$skat_p<alpha) / nrow(temp_data)
	power_mat$omni[i]<- sum(temp_data$omni_p<alpha) / nrow(temp_data)
	power_mat$BJ[i] <- sum(temp_data$BJ_p<alpha) / nrow(temp_data)
}


# Plot
xlab_tex = expression('Number of causal SNPs')
main_tex = expression(paste('40 SNPs on Chr5, Median |', rho[3], '|<0.1', sep=''))

postscript("S4_d40_chr5_lowrho.eps", width=720, height=500)
plot_data <- data.frame(power_mat[1:8,1:7])
par(mar=c(5,6,6,2))
xmax <- 8
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
ss_BJ <- smooth.spline(plot_data$BJ ~ plot_data$ncausal, spar=0.4)
lines(ss_GBJ, col=2, lwd=5, lty=1, cex=2)
lines(ss_GHC, col=3, lwd=5, lty=2, cex=2)
lines(ss_minP, col=4, lwd=5, lty=3, cex=2)
lines(ss_skat, col=1, lwd=5, lty=4, cex=2)
lines(ss_omni, col=6, lwd=5, lty=5, cex=2)
legend(x=6.1, y=0.5, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
		lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')
dev.off()
