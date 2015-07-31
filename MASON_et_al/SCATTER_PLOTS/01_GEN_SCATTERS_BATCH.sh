#!/bin/bash
set -e

# clean all
./00_CLEAN_THIS_DIR.sh

# get the anvi'o merged profile for mason metagenomes
wget http://files.figshare.com/2196643/MASON_YERGEAU_MG_MERGED.tar.gz
tar -zxvf MASON_YERGEAU_MG_MERGED.tar.gz

# edit these three varaibles to specify which `bins` to focus on,
# in a particular `collection` that can be found in the profile
# database in a merged directory (`path_to_the_merge`).
path_to_the_merge="MASON-YERGEAU-MG-MERGED"
collection="SUPERVISED"
bins="DWH_O_desum DWH_Cryptic DWH_Unknown"

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

rm -rf *PROFILE_*.txt *PROFILE_*.pdf

for comparison in `ls SAMPLES_* | awk 'BEGIN{FS="SAMPLES_"}{print $2}'`
do
    for bin in $bins
    do
        INFO "generating variability profile for $bin / $comparison"
        anvi-gen-variability-profile -p $path_to_the_merge/PROFILE.db \
                                     -a $path_to_the_merge/ANNOTATION.db \
                                     -c $collection \
                                     -b $bin \
                                     -n 5 \
                                     -o VARIABILITY_PROFILE_"$comparison"_"$bin".txt \
                                     -m 0 \
                                     -x 2\
                                     -S SAMPLES_"$comparison" \
                                     --quince


        INFO "generating scatter plot for $bin / $comparison"
        ./02_SCATTER_PLOT.R VARIABILITY_PROFILE_"$comparison"_"$bin".txt 1000
    done
done
