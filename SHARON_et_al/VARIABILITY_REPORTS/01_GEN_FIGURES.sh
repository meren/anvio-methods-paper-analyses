#!/bin/bash

set -e

# clean the directory
./00_CLEAN_THIS_DIR.sh

# get the merged profile:
wget http://files.figshare.com/2196633/INFANT_CLC_MERGED.tar.gz
tar -zxvf INFANT_CLC_MERGED.tar.gz

path_to_the_merged_dir="INFANT-CLC-MERGED"
collection="SUPERVISED"
bins="E_faecalis S_epidermidis_pan S_aureus"

C() {
    echo -e "\033[0;30m\033[46m$1\033[0m"
}

INFO() { 
    echo
    C "#"
    C "#"
    C "# $1"
    C "#"
    C "#"
    echo
}

for bin in $bins
do
    # first we generate a varaibiliity profile from SNPs with a scattering
    # power of three (-m 3)
    INFO "generating the profile for the analysis of variation for $bin"
    anvi-gen-variability-profile -p $path_to_the_merged_dir/PROFILE.db \
                                 -a $path_to_the_merged_dir/ANNOTATION.db \
                                 -c $collection \
                                 -b $bin \
                                 -n 5 \
                                 -o PROFILE_"$bin".txt \
                                 -m 3 \
                                 --quince \
                                 -S SAMPLES.txt

    # the reason we create this one is because we want to learn about all
    # reported splits during the profiling to understand the variation density
    # per kbp for each bin.
    INFO "generating the profile for the analysis of density for $bin"
    anvi-gen-variability-profile -p $path_to_the_merged_dir/PROFILE.db \
                                 -a $path_to_the_merged_dir/ANNOTATION.db \
                                 -c $collection \
                                 -b $bin \
                                 -n 0 \
                                 -o DENSITY_"$bin".txt \
                                 -m 0 \
                                 -S SAMPLES.txt
done


# time to visualize:
INFO "visualizing E_faecalis"
./02_GEN_FIGURE_SUMMARY.R PROFILE_E_faecalis.txt DENSITY_E_faecalis.txt 158 2870000

INFO "visualizing S_epidermidis_pan"
./02_GEN_FIGURE_SUMMARY.R PROFILE_S_epidermidis_pan.txt DENSITY_S_epidermidis_pan.txt 158 2610000

INFO "visualizing S_aureus"
./02_GEN_FIGURE_SUMMARY.R PROFILE_S_aureus.txt DENSITY_S_aureus.txt 158 2720000
       # this is how many unique nucleotide positions |              ^     ^
       # will be randomly sampled from the profile to |              |     |
       # make things a bit more visually comparable   |--------------/     |
       #                                                                   |
       # this is the genome size that is used to comp |                    |
       # variation density per kbp in a given bin, so |                    |
       # also must be changed manually... sorry...    |-------------------/
