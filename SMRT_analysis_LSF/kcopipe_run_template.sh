#!/bin/bash -l

umask 0002

. /broad/software/scripts/useuse

reuse GCC-4.4
reuse LSF
reuse Perl-5.8
reuse Python-2.7
reuse  Java-1.7
reuse .coreutils-8.22
reuse .samtools-0.1.19

/home/unix/bhaas/SVN/trinityrnaseq/trunk/Trinity --seqType fq --left 10M.left.fq --right 10M.right.fq --max_memory 20G  --CPU 5 --SS_lib_type RF --grid_conf /home/unix/bhaas/SVN/trinityrnaseq/trunk/htc_conf/BroadInst_LSF.regev.10.conf



