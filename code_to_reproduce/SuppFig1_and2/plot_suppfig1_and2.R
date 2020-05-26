# The purpose of this script is to plot Supplementary Figure 1.

# We read in the output produced by gen_data_suppfig4.R, write the correlation matrices
# for the genes, check the analytical p-value accuracy, and plot the accuracy.
# Output is the file pval_acc.png.

# Make sure all the jobs from gen_data_fig4.R have finished running,
# otherwise the necessary data will not exist and there will be errors.


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'
library(data.table, lib=libs)
library(ggplot2, lib=libs)
library(dplyr, lib=libs)
library(GBJ, lib=libs)

# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

###################################################################
# No manual changes need to be made after this point
###################################################################


########################################################################################
# Read in all the simulation results
setwd(data_dir)

# Should be 1000 files for each of 6 runs
chr5_dat <- fread('chr5_pval_acc_S2_aID1.txt', header=T)
for (file_it in 2:6000) {
    temp_fname <- paste0('chr5_pval_acc_S2_aID', file_it, '.txt')
    temp_file <- tryCatch(fread(temp_fname, header=T), warning=function(w) w,
                          error=function(e) e)
    if (class(temp_file)[1] %in% c('simpleWarning', 'simpleError')) {next}
    
    mylist <- list(chr5_dat, temp_file)
    chr5_dat <- rbindlist(mylist)
}
########################################################################################


########################################################################################
# Helper function, use this to get most extreme test statistic value that would
# generate a GBJ objective function of gbj, for a set of size d and with pairwise_cors.
stat_generating_gbj <- function(x, d, pairwise_cors, gbj) {
    GBJ_objective(t_vec=x, d=d, k_vec=c(1), pairwise_cors=pairwise_cors) - gbj
}
########################################################################################


########################################################################################
# These parameters help us locate our simulations for a specific Snum
leg_file <- fread('hapmap3.r2.b36.chr5.legend')
param_mat <- data.frame(start_row=c(1,1,75878,83510,22969,65400),
                        end_row=c(1,1,75897,83520,22987,65407)) %>%
    mutate(start_bp=leg_file$position[start_row]) %>%
    mutate(end_bp=leg_file$position[end_row]) %>%
    mutate(size = end_bp - start_bp) %>%
    mutate(Gene = c('FGFR2_low', 'FGFR2_high', 'RNF145', 'BNIP1', 'SEPP1','HARS2'))
keep_snps_list <- list(c(1), c(1), 2:20,
                       c(1,3,4,5,6,7,8,10,11),  c(3,4,8,10,12,14,19), c(2,3,4,8))
########################################################################################



########################################################################################
# Record the accuracy results in a data.frame for plotting
# For p-values of 0.001 to 0.01
numS <- 6
p_quantiles <- seq(from=0.99, to=0.999, length.out=10)
results_df <- data.frame(pvalue = rep(1-p_quantiles, numS), 
                         Setting=(rep(1:numS, each=length(p_quantiles))),
                         gbj=NA, gbj_p=NA, Gene=NA)
# loop through each parameter setting in results_df
for (row_it in 1:nrow(results_df)) {

    # temporary parameter settings
    temp_p <- results_df$pvalue[row_it]
    temp_s <- results_df$Setting[row_it]
    
    # select only the simulation runs for this Snum
    # also grab the correlation structure
    cor_struct <- read.table(paste0('cor_mat_S1_aID', results_df$Setting[row_it], '.txt', header=F))
    temp_idx <- which(chr5_dat$Gene == param_mat$Gene[temp_s])
    temp_data <- chr5_dat[temp_idx, ]
    temp_gbj <- quantile(temp_data$gbj, 1-temp_p)
  
    # which most extreme test statistic can get us this observed gbj
    # (the one that simulation tells us is at the 95th/99th/etc quantile)
    pairwise_cors <- cor_struct[lower.tri(cor_struct)]
    temp_Z <- uniroot(stat_generating_gbj, pairwise_cors=pairwise_cors, gbj=temp_gbj,
                      d=nrow(cor_struct), interval=c(0, 5))$root
    temp_stats <- c(temp_Z, rep(0, nrow(cor_struct)-1))
    
    # what is the analytical p-value for this 95th/99th/etc quantile gbj
    GBJ_output <- GBJ(test_stats=temp_stats, pairwise_cors=pairwise_cors)
    results_df$gbj[row_it] <- GBJ_output$GBJ
    results_df$gbj_p[row_it] <- GBJ_output$GBJ_p
    results_df$Gene[row_it] <- param_mat$Gene[temp_s]
    cat(row_it)
}


########################################################################################
# Plot p-values of 0.001 to 0.01

# The x-axis is Monte Carlo p-value and y-axis analytical p-value
plot_data <- results_df 
ggplot(data=plot_data, aes(x=pvalue, y=gbj_p,  
                           linetype=as.factor(Gene), color=as.factor(Gene))) + 
    geom_line() + geom_point(aes(shape=as.factor(Gene))) + 
    scale_x_continuous(breaks=seq(from=0.001, to=0.01, by=0.001)) +
    scale_y_continuous(breaks=seq(from=0.001, to=0.015, by=0.001)) +
    labs(color = "Gene", linetype = "Gene", shape="Gene") + 
    xlab('Simulation GBJ p-value') + ylab('Analytical GBJ p-value') + 
    geom_abline(slope=1, intercept=0)

# save
ggsave('pvalue_acc_01.png')


########################################################################################
# Record the accuracy results in a data.frame for plotting
# For p-values of 0.00001 to 0.0001
numS <- 6
p_quantiles <- seq(from=0.9999, to= 0.99999, length.out=10)
results_df <- data.frame(pvalue = rep(1-p_quantiles, numS), 
                         Setting=(rep(1:numS, each=length(p_quantiles))),
                         gbj=NA, gbj_p=NA, Gene=NA)
# loop through each parameter setting in results_df
for (row_it in 1:nrow(results_df)) {

    # temporary parameter settings
    temp_p <- results_df$pvalue[row_it]
    temp_s <- results_df$Setting[row_it]
    
    # select only the simulation runs for this Snum
    # also grab the correlation structure
    cor_struct <- read.table(paste0('cor_mat_S1_aID', results_df$Setting[row_it], '.txt', header=F))
    temp_idx <- which(chr5_dat$Gene == param_mat$Gene[temp_s])
    temp_data <- chr5_dat[temp_idx, ]
    temp_gbj <- quantile(temp_data$gbj, 1-temp_p)
  
    # which most extreme test statistic can get us this observed gbj
    # (the one that simulation tells us is at the 95th/99th/etc quantile)
    pairwise_cors <- cor_struct[lower.tri(cor_struct)]
    temp_Z <- uniroot(stat_generating_gbj, pairwise_cors=pairwise_cors, gbj=temp_gbj,
                      d=nrow(cor_struct), interval=c(0, 5))$root
    temp_stats <- c(temp_Z, rep(0, nrow(cor_struct)-1))
    
    # what is the analytical p-value for this 95th/99th/etc quantile gbj
    GBJ_output <- GBJ(test_stats=temp_stats, pairwise_cors=pairwise_cors)
    results_df$gbj[row_it] <- GBJ_output$GBJ
    results_df$gbj_p[row_it] <- GBJ_output$GBJ_p
    results_df$Gene[row_it] <- param_mat$Gene[temp_s]
    cat(row_it)
}
########################################################################################






########################################################################################
# Plot p-values of 0.00001 to 0.0001

# The x-axis is Monte Carlo p-value and y-axis analytical p-value
plot_data <- results_df 
ggplot(data=plot_data, aes(x=pvalue, y=gbj_p,  
                           linetype=as.factor(Gene), color=as.factor(Gene))) + 
    geom_line() + geom_point(aes(shape=as.factor(Gene))) + 
    scale_x_continuous(breaks=seq(from=0.00001, to=0.0001, by=0.00001)) +
    scale_y_continuous(breaks=seq(from=0.00001, to=0.00015, by=0.00001)) +
    labs(color = "Gene", linetype = "Gene", shape="Gene") + 
    xlab('Simulation GBJ p-value') + ylab('Analytical GBJ p-value') + 
    geom_abline(slope=1, intercept=0)

# save
ggsave('pvalue_acc_0001.png')





