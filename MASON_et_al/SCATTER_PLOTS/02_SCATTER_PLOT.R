#!/usr/bin/env Rscript

suppressPackageStartupMessages(library('ggplot2'))
suppressPackageStartupMessages(library('reshape2'))
suppressPackageStartupMessages(library('reshape'))
suppressPackageStartupMessages(library('gridExtra'))
suppressPackageStartupMessages(library('plyr'))
suppressPackageStartupMessages(library('gtools'))

cols <- c("AA" = "orange",
		"TT" = "pink",
		"CC" = "yellow",
		"GG" = "green",
		'N' = 'black',
		'AC' = "#89C5DA",
		'CA' = "#89C5DA",
		'AG' = "#DA5724",
		'GA' = "#DA5724",
		'AT' = "#74D944",
		'TA' = "#74D944",
		'CG' = "#CE50CA",
		'GC' = "#CE50CA",
		'CT' = "#5F4991",
		'TC' = "#5F4991",
		'TG' = "#C0717C",
		'GT' = "#C0717C")


lm_eqn <- function(df){
   m <- lm(y ~ x, df);
   f <- summary(m)$fstatistic
   p <- pf(f[1],f[2],f[3],lower.tail=F)
   attributes(p) <- NULL
   label=paste('r2:', format(summary(m)$r.squared, digits = 3), 'f-statistic:', round(f[1], digits=1), sep = ' ')
   as.character(label);
}


scatter <- function(profile, num_entries, width = 12, height = 12){
    if(invalid(num_entries)){
        num_entries = 0
    }

    variability_profile <- as.data.frame(read.csv(profile, header=TRUE, sep="\t"))

	samples = levels(variability_profile$sample_id)
	if(!(length(samples) == 2)){
		print('There must be only two samples in the data frame...')
		exit()
	}

	# keep the number of eligable positions of nucleotide variation somewhere:
	total_num_unique_pos <- length(unique(variability_profile$unique_pos_identifier))
	
	# setup the transversion_vs_transition data frame:
	# get the df for competing nt identitis before subsampling:
	competing = variability_profile[variability_profile$n2n1ratio > 0, ]
	competing_nts_df <- count(competing, "competing_nts")
	competing_nts_df$x <- 'all'

	transitions = c('AG', 'CT', 'GA', 'TC')
	transversion_vs_transition <- data.frame('mtype' = character(), 'freq' = integer(), stringsAsFactors=FALSE)
	transversion_vs_transition[1, ] <- c('transversion', 0)
	transversion_vs_transition[2, ] <- c('transition', 0)
	transversion_vs_transition[transversion_vs_transition$mtype == "transition", ]$freq = sum(competing_nts_df[competing_nts_df$competing_nts %in% transitions, ]$freq)
	transversion_vs_transition[transversion_vs_transition$mtype == "transversion", ]$freq = sum(competing_nts_df[!(competing_nts_df$competing_nts %in% transitions), ]$freq)
	transversion_vs_transition$freq <- as.numeric(transversion_vs_transition$freq)
	transversion_vs_transition$x <- 'all'
	transversion_vs_transition_ratio <- transversion_vs_transition$freq[2] / transversion_vs_transition$freq[1]

    # subsample the df randomly if necessary:
    if(num_entries > 0){
        if(num_entries > length(unique(variability_profile$unique_pos_identifier))){
            variability_profile_subsampled <- variability_profile
        } else {
            variability_profile_subsampled <- variability_profile[variability_profile$unique_pos_identifier %in% sample(unique(variability_profile$unique_pos_identifier), num_entries), ]
        }
    } else {
        variability_profile_subsampled <- variability_profile
    }

	
    new_df = data.frame('pos' = character(), 'x' = numeric(), 'y' = numeric(), 'var' = numeric(), competing_nts = character(), stringsAsFactors = FALSE)
    N = 1
    for(pos in unique(variability_profile_subsampled$unique_pos_identifier)){
        m = variability_profile_subsampled[variability_profile_subsampled$unique_pos_identifier == pos, ]
        new_df[N, ] <- c(m$unique_pos_identifier[1],
					     x=m[m$sample_id == samples[1], ]$n2n1ratio,
					     y=m[m$sample_id == samples[2], ]$n2n1ratio,
						 var = 1 - sd(m$coverage / max(m$coverage)),
						 competing_nts = as.character(m$competing_nts[1]))
        N = N + 1
    }
	new_df$x <- as.numeric(new_df$x)
	new_df$y <- as.numeric(new_df$y)
	new_df$var <- as.numeric(new_df$var)

	points  <- data.frame(x = c(0, 0, 1, 1), y = c(0, 1, 0, 1))

    g <- ggplot(new_df, aes(x=x, y=y))
    g <- g + geom_point(aes(alpha = var, size=var, color=competing_nts)) 
	g <- g + scale_fill_manual(values = cols)
	g <- g + scale_color_manual(values = cols, guide = guide_legend(override.aes=aes(fill=NA)))
    g <- g + stat_smooth(method = "glm", color="black", size=3, formula = y ~ x) 
    g <- g + annotate('text', hjust=0, x = 0.05, y = 0.9, label = paste(profile, lm_eqn(new_df), paste('N:', total_num_unique_pos, 'n:', nrow(new_df), 'k:', round(transversion_vs_transition_ratio, digits=2), sep=' '), sep='\n'), size=7, color='red') 
	g <- g + xlab(samples[1]) + ylab(samples[2])
    g <- g + scale_y_log10() + scale_x_log10()
    g <- g + theme(legend.position = c(0.1, 0.8), 
			legend.justification = c(0.1, 0.8), 
			legend.background = element_rect(colour = NA, fill = "white"), axis.text.y = element_text(angle = 90))
	g <- g + geom_point(data = points, colour = "red", size=0)
    
	
    pdf(paste(profile, '.pdf', sep=""), width=width, height=height)
    print(g)
    dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
profile <- args[1]
num_entries <- as.integer(args[2])

scatter(profile, num_entries)
