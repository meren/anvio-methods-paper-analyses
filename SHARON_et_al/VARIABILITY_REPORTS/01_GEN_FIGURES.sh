#!/bin/bash

set -e

# this file will first generate reports of variabliilty profiles for
# each bin in `bins` described under the `collection` id in the merged
# profile database that will be found in `path_to_the_merged_dir`.
path_to_the_merged_dir="/Users/meren/papi-stuff/INFANT-CLC-MERGED"
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
rm -rf *DENSITY*.txt
rm -rf *PROFILE*.txt
rm -rf *pdf

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
                                 -S 00_SAMPLES

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
                                 -S 00_SAMPLES
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
