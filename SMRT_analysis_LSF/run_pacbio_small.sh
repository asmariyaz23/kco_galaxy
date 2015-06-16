#!/bin/bash -l

umask 0002

. /broad/software/scripts/useuse

use Samtools 
use .bedtools-2.17.0 
use GCC-4.9 
use LSF 

unset PYTHONPATH 
source /seq/regev_genome_portal/SOFTWARE/SMRT_analysis_LSF/current/etc/setup.sh 

##Store input parameters 
LIBRARY=$1
LIBRARYNAME=$2
LIBSUMMARY=$3
OUTPUTDIR=$4
##echo $OUTPUTDIR
##echo $LIBRARYNAME
##echo $LIBRARY

##Logic for appending a / to output directory if it doesn't exist and later concatenating directory path with pipe.out and pipe.err files

OUTPUTDIRLASTCHAR="${OUTPUTDIR: -1}"
if [ $OUTPUTDIRLASTCHAR != "/" ]
then
OUTPUTDIR="$OUTPUTDIR/"
fi

##initializing names of out and err files
PIPEOFILENAME=$LIBRARYNAME".o"
PIPEEFILENAME=$LIBRARYNAME".e"

##concatenating output directory and file name
PIPERRFILE="$OUTPUTDIR$PIPEEFILENAME"
PIPEOUTFILE="$OUTPUTDIR$PIPEOFILENAME"

echo $PIPERRFILE
echo $PIPEOUTFILE

##PacBio pipeline command
bsub -q regevlab -P pacbio -e $PIPERRFILE -o $PIPEOUTFILE -R "rusage[mem=10]" -K python /seq/regev_genome_portal/SOFTWARE/galaxy/tools/SMRT_analysis_LSF/launch_library_small.py --lib_summary $LIBSUMMARY --lib_table $LIBRARY --out_path $OUTPUTDIR

##for i in $(find $OUTPUTDIR -maxdepth 1); do
  ##if [ $i = "junctions_rarefaction.png" || $i = "loci_rarefaction.png" || $i = "lib_summary.txt" ]
     ##then
     ##mv $i ..
  ##fi
##done


