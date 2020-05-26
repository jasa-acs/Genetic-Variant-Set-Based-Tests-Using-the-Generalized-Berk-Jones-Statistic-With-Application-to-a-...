# The purpose of this script is to plot Supplementary Figures 7 and 8 of
# the submitted manuscript.

# We read in the output produced by gen_data_tab2.R and set up the plotting windows with various options.
# For example, we change options like the title of the plot as well as the colors and textures
# of the lines.
# Outputs are 5 png files (all QQplots).

# Make sure all the jobs from gen_data_tab2.R have finished running,
# otherwise the necessary data will not exist and there will be errors.

# Recall that since we are not allowed to share the original CGEMS data, it will not be possible to reproduce
# our results exactly unless you have procured the CGEMS dataset from dbGaP and used that data to run gen_data_tab2.R.
# Using the fake data we provide will not give you the output needed to reproduce our QQplots.


########################################################################################
# Change this location to the directory that holds the necessary R libraries
# (must be installed first).
libs <- '/n/home/user/Rlibrary/3.3.1'
library(data.table, lib=libs)
library(dplyr, lib=libs)

# Change data_dir your working directory directory which holds all the necessary data files.
data_dir <- '/users/user/desktop'

###################################################################
# No manual changes need to be made after this point
###################################################################

# Read in the CGEMS results
cgems_results <- read.table('cgems_results_S1_1.txt', header=T)
for(i in 2:181)
{
    temp_fname <- paste('cgems_results_S1_', i, '.txt', sep='')
    temp_file <- read.table(temp_fname, header=T)
    cgems_results <- rbind(cgems_results, temp_file)
}


###################################################################
# Function to draw QQplots
make_qqplot <- function(mainName, fileName, GBJ_p, GHC_p, SKAT_p, minP_p, OMNI_p) {
    # Sort the -log10 pvalues
    exp1 <- -log10(1:length(GBJ_p) /(length(GBJ_p)+1)) 
    exp2 <- -log10(1:length(GHC_p) /(length(GHC_p)+1)) 
    exp3 <- -log10(1:length(SKAT_p) /(length(SKAT_p)+1)) 
    exp4 <- -log10(1:length(minP_p) /(length(minP_p)+1)) 
    exp5 <- -log10(1:length(which(!is.na(OMNI_p))) /(length(which(!is.na(OMNI_p)))+1)) 
    o1 <- -log10(sort(GBJ_p))
    o2 <- -log10(sort(GHC_p))
    o3 <- -log10(sort(SKAT_p))
    o4 <- -log10(sort(minP_p))
    o5 <- -log10(sort(OMNI_p))
    
    # Define axes
    xlabel = expression(paste(-log[10], "(expected p-value)", sep=""))
    ylabel = expression(paste(-log[10], "(observed p-value)", sep=""))
    xlimits = c(0,6)
    ylimits = c(0,6)
    
    # Open the file
    png(file = fileName, width = 820, height = 820)
    par(mar=c(4.5,5.5,3,1.5)+1.0)
    
    # Plot each of the p-values
    plot(exp1, o1, xlab ='', ylab = ylabel, main = mainName, xlim=xlimits, ylim=ylimits,
         cex.main = 4, cex.lab=3, cex.axis=2, family = "serif", col = "red2", pch = 1)
    title(xlab=xlabel, line=4, cex.lab=3)
    par(new=T)
    plot(exp2, o2, xlab = '', ylab = ylabel, main = '', xlim=xlimits, ylim=ylimits,
         cex.main = 4, cex.lab=3, cex.axis=2, family = "serif", col = "blue", pch = 2)
    par(new=T)
    plot(exp3, o3, xlab ='', ylab = ylabel, main = '', xlim=xlimits, ylim=ylimits,
         cex.main = 4, cex.lab=3, cex.axis=2, family = "serif", col = "green2", pch = 3)
    par(new=T)
    plot(exp4, o4, xlab ='', ylab = ylabel, main = '', xlim=xlimits, ylim=ylimits,
         cex.main = 4, cex.lab=3, cex.axis=2, family = "serif", col = "orange", pch = 4)
    par(new=T)
    plot(exp5, o5, xlab ='', ylab = ylabel, main = '', xlim=xlimits, ylim=ylimits,
         cex.main = 4, cex.lab=3, cex.axis=2, family = "serif", col = "darkgrey", pch = 5)
    abline(0,1)
    
    # Draw the legend
    leg1 = bquote("GBJ")
    leg2 = bquote("GHC")
    leg3 = bquote("SKAT")
    leg4 = bquote("MinP")
    leg5 = bquote("OMNI")
    legend(x=4, y=2, c("GBJ", "GHC", "SKAT", "MinP", "OMNI"), 
           col=c('red2', 'blue', 'green2', 'orange', 'darkgrey'), lty=1, pch=1:5, lwd=6, cex=2.5, bty='n')
    
    dev.off()
}


###################################################################
# Loop through four quantiles, draw the QQplots
# Empirical distribution of gene p-values in CGEMS analysis
alpha_vec <- c(0.05, 0.01, 0.001)
quantile_vec <- c(0, quantile(cgems_results$num_snps, prob=c(0.25, 0.5, 0.75, 1)))
alpha_tab <- c()
for (quantile_idx in 1:4) {
    # split results into quantiles
    temp_results<- cgems_results %>% 
            filter(num_snps > quantile_vec[quantile_idx] & num_snps <= quantile_vec[quantile_idx+1])
      
    # QQplot
    make_qqplot(mainName=paste0('Q', quantile_idx), fileName=paste0('Q', quantile_idx, 'plot.png'), 
                GBJ_p=temp_results$gBJ_p, GHC_p=temp_results$GHC_p, 
                SKAT_p=temp_results$skat_p, minP_p=temp_results$minP_pvalue, 
                OMNI_p=temp_results$omni_p)
}

# overall QQ-plot
make_qqplot(mainName='All Genes', fileName='allQQ.png', 
            GBJ_p=cgems_results$gBJ_p, GHC_p=cgems_results$GHC_p, 
            SKAT_p=cgems_results$skat_p, minP_p=cgems_results$minP_pvalue, 
            OMNI_p=cgems_results$omni_p)

