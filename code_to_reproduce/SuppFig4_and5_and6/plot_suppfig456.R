# The purpose of this script is to plot the figures in Supplementary Figures 4,5,6 of
# the submitted manuscript.

# We read in the output produced by gen_data_suppfig456.R, use splines to generate smoothed power 
# curves, and set up the plotting windows with various options.
# For example, we change options like the title of the plot as well as the colors and textures
# of the lines.

# Make sure all the jobs from gen_data_suppfig456.R have finished running,
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
results1k <- NULL
results2k <- NULL
results4k <- NULL
alpha <- 3.34 * 10^(-6)
for (i in 1:40)
{
	# There should be 120 files.
    temp_fname1k <- paste('chr5_power_S1_', i, '.txt', sep='')
    temp_fname2k <- paste('chr5_power_S2_', i, '.txt', sep='')
    temp_fname4k <- paste('chr5_power_S4_', i, '.txt', sep='')

	temp_file1k <- tryCatch(read.table(temp_fname1k, header=T, sep='\t'), warning=function(w) w, 
						error=function(e) e)
	temp_file2k <- tryCatch(read.table(temp_fname2k, header=T, sep='\t'), warning=function(w) w, 
						error=function(e) e)
	temp_file4k <- tryCatch(read.table(temp_fname4k, header=T, sep='\t'), warning=function(w) w, 
						error=function(e) e)
	
 	# Sometimes GBJ will return NA due to lack of numerical precision in R, but we can
 	# tell if the p-value would be low enough from previous runs which did not run into 
 	# numerical precision issues. Corresponds to rerunning and getting an upper bound on the p-value.
	all_genes1k <- unique(temp_file1k$start_row)
	all_genes2k <- unique(temp_file2k$start_row)
	all_genes4k <- unique(temp_file4k$start_row)
	for (j in 1:length(all_genes1k)) {
	    temp_gene1k <- all_genes1k[j]
	    temp_gene2k <- all_genes2k[j]
	    temp_gene4k <- all_genes4k[j]
	    gene_data1k <- temp_file1k[which(temp_file1k$start_row==temp_gene1k), ]
	    gene_data2k <- temp_file2k[which(temp_file2k$start_row==temp_gene2k), ]
	    gene_data4k <- temp_file4k[which(temp_file4k$start_row==temp_gene4k), ]
	    good_GBJ_data1k <- gene_data1k[which(gene_data1k$GBJ_err == 0), ]
	    good_GBJ_data2k <- gene_data2k[which(gene_data2k$GBJ_err == 0), ]
	    good_GBJ_data4k <- gene_data4k[which(gene_data4k$GBJ_err == 0), ]
	    GBJ_threshold1k <- good_GBJ_data1k$GBJ[which(good_GBJ_data1k$GBJ_p < alpha)]
	   	GBJ_threshold2k <- good_GBJ_data2k$GBJ[which(good_GBJ_data2k$GBJ_p < alpha)]
	   	GBJ_threshold4k <- good_GBJ_data4k$GBJ[which(good_GBJ_data4k$GBJ_p < alpha)]
	    
	    # If we have good runs from before that did not run into numerical precision errors.
	    if (length(GBJ_threshold1k) > 0) {
	        GBJ_threshold1k <- min(GBJ_threshold1k)
	        fixable_rows1k <- which(gene_data1k$GBJ_err != 0 & gene_data1k$GBJ > GBJ_threshold1k)
	        if (length(fixable_rows1k) > 0) {
	            gene_data1k$GBJ_p[fixable_rows1k] <- alpha * 0.1
	            temp_file1k[which(temp_file1k$start_row==temp_gene1k), ] <- gene_data1k
	        }
	    }
	    if (length(GBJ_threshold2k) > 0) {
	        GBJ_threshold2k <- min(GBJ_threshold2k)
	        fixable_rows2k <- which(gene_data2k$GBJ_err != 0 & gene_data2k$GBJ > GBJ_threshold2k)
	        if (length(fixable_rows2k) > 0) {
	            gene_data2k$GBJ_p[fixable_rows2k] <- alpha * 0.1
	            temp_file2k[which(temp_file2k$start_row==temp_gene2k), ] <- gene_data2k
	        }
	    }
	    if (length(GBJ_threshold4k) > 0) {
	        GBJ_threshold4k <- min(GBJ_threshold4k)
	        fixable_rows4k <- which(gene_data4k$GBJ_err != 0 & gene_data4k$GBJ > GBJ_threshold4k)
	        if (length(fixable_rows4k) > 0) {
	            gene_data4k$GBJ_p[fixable_rows4k] <- alpha * 0.1
	            temp_file4k[which(temp_file4k$start_row==temp_gene4k), ] <- gene_data4k
	        }
	    }
	}
		
	# Sometimes GHC will return NA due to lack of numerical precision in R, but we can
 	# tell if the p-value would be low enough from previous runs which did not run into 
 	# numerical precision issues. Corresponds to rerunning and getting an upper bound on the p-value.
	for (j in 1:length(all_genes1k)) {
	    temp_gene1k <- all_genes1k[j]
	    temp_gene2k <- all_genes2k[j]
	    temp_gene4k <- all_genes4k[j]
	    gene_data1k <- temp_file1k[which(temp_file1k$start_row==temp_gene1k), ]
	    gene_data2k <- temp_file2k[which(temp_file2k$start_row==temp_gene2k), ]
	    gene_data4k <- temp_file4k[which(temp_file4k$start_row==temp_gene4k), ]
	    good_GHC_data1k <- gene_data1k[which(gene_data1k$GHC_err == 0), ]
	    good_GHC_data2k <- gene_data2k[which(gene_data2k$GHC_err == 0), ]
	    good_GHC_data4k <- gene_data4k[which(gene_data4k$GHC_err == 0), ]
	    GHC_threshold1k <- good_GHC_data1k$GBJ[which(good_GHC_data1k$GHC_p < alpha)]
	   	GHC_threshold2k <- good_GHC_data2k$GBJ[which(good_GHC_data2k$GHC_p < alpha)]
	   	GHC_threshold4k <- good_GHC_data4k$GBJ[which(good_GHC_data4k$GHC_p < alpha)]
	    
	  	# If we have good runs from before that did not run into numerical precision errors.
	    if (length(GHC_threshold1k) > 0) {
	        GHC_threshold1k <- min(GHC_threshold1k)
	        fixable_rows1k <- which(gene_data1k$GHC_err != 0 & gene_data1k$GHC > GHC_threshold1k)
	        if (length(fixable_rows1k) > 0) {
	            gene_data1k$GHC_p[fixable_rows1k] <- alpha * 0.1
	            temp_file1k[which(temp_file1k$start_row==temp_gene1k), ] <- gene_data1k
	        }
	    }
	    if (length(GHC_threshold2k) > 0) {
	        GHC_threshold2k <- min(GHC_threshold2k)
	        fixable_rows2k <- which(gene_data2k$GHC_err != 0 & gene_data2k$GHC > GHC_threshold2k)
	        if (length(fixable_rows2k) > 0) {
	            gene_data2k$GHC_p[fixable_rows2k] <- alpha * 0.1
	            temp_file2k[which(temp_file2k$start_row==temp_gene2k), ] <- gene_data2k
	        }
	    }
	    if (length(GHC_threshold4k) > 0) {
	        GHC_threshold4k <- min(GHC_threshold4k)
	        fixable_rows4k <- which(gene_data4k$GHC_err != 0 & gene_data4k$GHC > GHC_threshold4k)
	        if (length(fixable_rows4k) > 0) {
	            gene_data4k$GHC_p[fixable_rows4k] <- alpha * 0.1
	            temp_file4k[which(temp_file4k$start_row==temp_gene4k), ] <- gene_data4k
	        }
	    }
	}
	
	# Attach to results
	if (is.null(results1k)) {
	    results1k <- temp_file1k
	   	results2k <- temp_file2k
	    results4k <- temp_file4k
	} else {
	    mylist1k <- list(results1k, temp_file1k)
	    results1k <- rbindlist(mylist1k, use.names=F)
		mylist2k <- list(results2k, temp_file2k)
	    results2k <- rbindlist(mylist2k, use.names=F)
	    mylist4k <- list(results4k, temp_file4k)
	    results4k <- rbindlist(mylist4k, use.names=F)	    			
	}
}



################################################################################
# First plot n=1k results
################################################################################
alpha <- 3.34 * 10^(-6)
temp_results <- results1k
power_mat <- matrix(data=NA, nrow=10, ncol=7)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni', 'BJ')
for (i in 1:8)
{
	# Find the power for each ncausal
	temp_data <- temp_results[which(temp_results$ncausal==i),]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	# Record
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
main_tex = expression(paste('40 SNPs on Chr5, All Data, 1000 subjects', sep=''))

postscript("n1000_d40_chr5_all.eps", width=720, height=500)
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
legend(x=6, y=0.42, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
       lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')
dev.off()




################################################################################
# Next plot n=2k results
################################################################################
temp_results <- results2k
power_mat <- matrix(data=NA, nrow=10, ncol=7)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni', 'BJ')
for (i in 1:8)
{
	# Find the power for each ncausal
	temp_data <- temp_results[which(temp_results$ncausal==i),]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	# Record
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
main_tex = expression(paste('40 SNPs on Chr5, All Data, 2000 subjects', sep=''))

postscript("n2000_d40_chr5_all.eps", width=720, height=500)
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
legend(x=6, y=0.42, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
       lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')
dev.off()



################################################################################
# Last plot n=4k results
################################################################################
temp_results <- results4k
power_mat <- matrix(data=NA, nrow=10, ncol=7)
power_mat <- data.frame(power_mat)
colnames(power_mat) <- c('ncausal','GBJ', 'GHC',
							'minP', 'skat', 'omni', 'BJ')
for (i in 1:8)
{
	# Find the power for each ncausal
	temp_data <- temp_results[which(temp_results$ncausal==i),]
	power_mat$ncausal[i] <- temp_data$ncausal[1]
	
	# Record
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
main_tex = expression(paste('40 SNPs on Chr5, All Data, 4000 subjects', sep=''))

postscript("n4000_d40_chr5_all.eps", width=720, height=500)
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
legend(x=6, y=0.42, c('GBJ', 'GHC', 'MinP', 'SKAT', 'OMNI'), col=c(2,3,4,1,6), 
       lwd=c(3,3,3,3,3), lty=c(1, 2, 3, 4, 5), cex=2, bty='n')
dev.off()




