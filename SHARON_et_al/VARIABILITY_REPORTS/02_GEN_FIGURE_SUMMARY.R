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

transversion_vs_transition_colors <- c('transversion' = "#451910",
								       'transition' = '#8C8983')


draw <- function (profile, density, num_entries, genome_size, width = 10, height = 20){
	# profile = 'PROFILE_E_faecalis.txt'
	# density = 'DENSITY_E_faecalis.txt'
	# num_entries = 158
	# genome_size = 2870000
    if(invalid(num_entries)){
        num_entries = 0
    }

    variability_profile <- as.data.frame(read.csv(profile, header=TRUE, sep="\t"))
	variability_density <- as.data.frame(read.csv(density, header=TRUE, sep="\t"))

	# keep the number of eligable positions of nucleotide variation somewhere:
	total_num_unique_pos <- length(unique(variability_profile$unique_pos_identifier))

	# get the df for competing nt identitis before subsampling:
	competing = variability_profile[variability_profile$n2n1ratio > 0, ]
	competing_nts_df <- count(competing, "competing_nts")
	competing_nts_df$x <- 'all'

	# get the df for number of SNPs per kb in each sample
	samples <- levels(variability_profile$sample_id)
	variability_density_per_sample_df <- count(variability_density, "sample_id")
	variability_density_per_sample_df$freq <- variability_density_per_sample_df$freq / (genome_size / 1000)
	for(sample in samples){
		if(!(sample %in% variability_density_per_sample_df$sample_id)){
			dd <- data.frame('sample_id' = sample, 'freq' = 0.0, stringsAsFactors = FALSE)
			variability_density_per_sample_df <- rbind(variability_density_per_sample_df, dd)
		}
	}
	

	# setup the transversion_vs_transition data frame:
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

	# competing nucleotide identities
	b <- ggplot(competing_nts_df, aes(x=factor(x),y=freq, fill=competing_nts, color=competing_nts)) + geom_bar(stat="identity")
	b <- b + scale_color_manual(values = cols, guide = guide_legend(override.aes=aes(fill=NA)))
	b <- b + scale_fill_manual(values = cols)
	b <- b + coord_flip()
	b <- b + xlab(paste('n:', total_num_unique_pos, sep=" "))
    b <- b + theme_bw() + theme(axis.text.y=element_blank(),
                                axis.text.y = element_blank(),
                                axis.ticks.y=element_blank(),
                                legend.position="none",
                                axis.line.y = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),
                                axis.title.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.text.x = element_blank(),
                                plot.margin=unit(c(-0.20,-0.10,-0.20,-0.10), "cm"))

	# transversion vs translation
	d <- ggplot(transversion_vs_transition, aes(x=factor(x), y=freq)) + geom_bar(aes(fill=mtype), stat="identity")
	d <- d + scale_color_manual(values = transversion_vs_transition_colors, guide = guide_legend(override.aes=aes(fill=NA)))
	d <- d + scale_fill_manual(values = transversion_vs_transition_colors)
	d <- d + coord_flip()
	d <- d + xlab(paste('k:', round(transversion_vs_transition_ratio, digits=2), sep=" "))
	d <- d + theme_bw() + theme(axis.text.y=element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y=element_blank(),
			legend.position="bottom",
			axis.line.y = element_line(colour = "black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_blank(),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x = element_blank(),
            plot.margin=unit(c(-0.20,-0.10,-0.20,-0.10), "cm"))

	g <- ggplot(data=variability_profile_subsampled, aes(x = sample_id, y = factor(unique_pos_identifier)))
    g <- g + geom_tile(aes(group=competing_nts, fill=competing_nts, alpha=sqrt(n2n1ratio)))
    g <- g + scale_x_discrete(expand = c(0, 0))
    g <- g + scale_y_discrete(expand = c(0, 0))
    g <- g + theme_bw() + theme(axis.text.y=element_blank(),
                                axis.text.y = element_text(size = 6),
                                axis.ticks.y=element_blank(),
                                legend.position="bottom",
                                axis.line.y = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),
                                axis.title.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12),
                                plot.margin=unit(c(-0.2,1,1,1), "cm"))
    g <- g + scale_color_manual(values = cols, guide = guide_legend(override.aes=aes(fill=NA)))
    g <- g + scale_fill_manual(values = cols)
    g <- g + scale_alpha(limits = c(0.0, mean(variability_profile_subsampled[variability_profile_subsampled$n2n1ratio > 0, ]$n2n1ratio)))


    p <- ggplot(data=variability_profile, aes(x = sample_id, y = coverage))
    p <- p + geom_violin(alpha=0.9)
    p <- p + geom_jitter(data=variability_profile, aes(x = sample_id, y = coverage), position = position_jitter(width=c(0.4)), alpha=0.3, size=0.8)
    p <- p + scale_y_log10()
    p <- p + theme_bw() + theme(axis.line.y = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_rect(fill = "#EAEAEA"),
                                axis.text.x=element_blank(),
                                axis.title.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                axis.ticks.y=element_blank(),
            					plot.margin=unit(c(-0.20,1,-0.20,1), "cm"))

	# variability denisty bars
	x <- ggplot(variability_density_per_sample_df, aes(x=factor(sample_id), y=freq)) + geom_bar(stat="identity")
	x <- x + ylab(paste('max:', round(max(variability_density_per_sample_df$freq), digits=2), sep=" "))
	x <- x + theme_bw() + theme(
			axis.text.y=element_blank(),
			axis.text.y = element_blank(),
			axis.ticks.y=element_blank(),
			legend.position="none",
			axis.line.y = element_line(colour = "black"),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_rect(fill = "#EAEAEA"),
			axis.title.x=element_blank(),
			axis.ticks.x=element_blank(),
			axis.text.x = element_blank(),
            plot.margin=unit(c(-0.20,1,-0.20,1), "cm"))


    xo <- ggplot_gtable(ggplot_build(x))
    po <- ggplot_gtable(ggplot_build(p))
    go <- ggplot_gtable(ggplot_build(g))
    bo <- ggplot_gtable(ggplot_build(b))
    do <- ggplot_gtable(ggplot_build(d))
    maxWidth = unit.pmax(po$widths[2:3], go$widths[2:3], bo$widths[2:3], do$widths[2:3], xo$widths[2:3])
    po$widths[2:3] <- maxWidth
    go$widths[2:3] <- maxWidth
    bo$widths[2:3] <- maxWidth
    do$widths[2:3] <- maxWidth
    xo$widths[2:3] <- maxWidth

    pdf(paste(profile, '.pdf', sep=""), width=width, height=height)
    grid.arrange(po, xo, go, bo, do, heights=c(2/15, 2/15, 9/15, 1/15, 1/15), ncol=1)
    dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
profile <- args[1]
density <- args[2]
num_entries <- as.integer(args[3])
genome_size <- as.integer(args[4])

draw(profile, density, num_entries, genome_size)

#setwd('/Users/meren/papi-stuff/anvio-methods-paper-analyses/SHARON_et_al/variability_analysis')
#draw('E_faecalis.txt', 158, 2870000)
#draw('S_epidermidis_pan.txt', 158, 2610000)
#draw('S_aureus.txt', 158, 2720000)

