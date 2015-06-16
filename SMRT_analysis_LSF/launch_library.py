#!/usr/bin/env python


# -*- coding: utf-8 -*-
"""
Created on Thu Dec 04 01:27:31 2014

@author: Eran
"""

import os, sys
import argparse
import subprocess
import pandas
import re
from time import sleep
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg') #Must come before import pyplot to save figures on server
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()


parser.add_argument('--lib_table', type=str, help='Path to library to cells \
				table')



parser.add_argument('--out_path', type=str, help='Path where library output \
				folder will be created')


parser.add_argument('--sum_lib', action='store_true', help='Specify this to \
only repeat library summary and avoid running full analysis. All appropriate \
cell output must be in place. Make sure to first remove all existing library \
summary data in the library folder')
parser.add_argument('--cell_script', type=str, default=\
'/broad/moothalab/sandbox/eranmick/PacBio/IsoSeq/launch_isoseq.py')
parser.add_argument('--queue', type=str, default=\
'regevlab')
parser.add_argument('--mem', type=str, default=\
'50')
parser.add_argument('--cufflinks', type=str, default=\
'/broad/moothalab/third_party/cufflinks-2.2.1/cufflinks-2.2.1.Linux_x86_64')
parser.add_argument('--lib_summary', type=str, default=\
'lib_summary.txt')


parser.add_argument('--junctions_rarefaction', type=str, default=\
'junctions_rarefaction.png')
parser.add_argument('--loci_rarefaction', type=str, default=\
'loci_rarefaction.png')


args = parser.parse_args()


def run_cmd(cmd):
	sys.stderr.write(cmd+'\n\n')
	out=subprocess.check_output(cmd, shell=True)
	return(out)

mouse = ['/broad/moothalab/sandbox/eranmick/PacBio/IsoSeq/\
mouse_brain/gmap_genome_dir', 'mm10', '/broad/moothalab/sandbox/eranmick/\
PacBio/IsoSeq/mouse_brain/compare_to_annotation/mm10.refseq.genes.from.ucsc.gtf']

human = ['/broad/moothalab/sandbox/eranmick/PacBio/IDP/\
gmap_genome_dir', 'hg19_ucsc', '/broad/moothalab/sandbox/eranmick/PacBio/\
IsoSeq/HEK293_1_2_UMASS/hg19.refseq.genes.from.ucsc.gtf']

cell_template = 'bsub -n 4 -q %s -P pacbio -o %s.out -e %s.err -R\"rusage[mem=\
%s] span[hosts=1]\" \"python %s %s %s %s %s %s %s %s\"'

JOB_ID_RE = re.compile('Job <(\d+)>')
JOB_STAT_RE = re.compile('Status <(\w+)>')

def make_cell_cmd(row):
	params = [args.queue, str(row.Library)+'_'+str(row.Cell), \
	str(row.Library)+'_'+str(row.Cell), args.mem, args.cell_script, \
	row.Data, args.out_path, row.Library, str(row.Cell)]
	if(row.Organism.lower() == 'mouse'):
		params = params + mouse
	if(row.Organism.lower() == 'human'):
		params = params + human
	cmd = cell_template % tuple(params)
	job_submit = run_cmd(cmd)
	m = JOB_ID_RE.match(job_submit)
	job_id = m.group(1)
	return(job_id)

def check_lsf(job_dict, interval):
	while(True):
		sleep(interval)
		remaining_jobs = [x for x in job_dict if (job_dict[x] != 'DONE')]
		check_cmds = ['bjobs -l '+x for x in remaining_jobs]
		status_str = [subprocess.check_output(x, shell=True) for x in \
											check_cmds]
		matches = [JOB_STAT_RE.search(x) for x in status_str]
		status = [x.group(1) for x in matches]
		for i,x in enumerate(remaining_jobs):
			job_dict[x] = status[i]
		if(all([job_dict[x]=='DONE' for x in job_dict])):
			return(status)

lib_table = pandas.read_table(args.lib_table)
lib_name = lib_table.ix[0,'Library']
out_path = args.out_path 
out_dir = args.out_path+'/'+lib_name
os.makedirs(out_dir)
os.chdir(out_dir)

if(not(args.sum_lib)):
	cell_job_ids = list(lib_table.apply(make_cell_cmd, axis=1).values)
	job_dict = dict([(x,'RUN') for x in cell_job_ids])
	check_lsf(job_dict, 1800)

cells = lib_table.Cell.values.astype('str')
lib_summary = pandas.concat([pandas.read_table(x+'/SMRTanalysis/summary.txt') \
										   for x in cells])
lib_summary_path = os.path.join(out_path,args.lib_summary)
lib_summary.to_csv(lib_summary_path,sep='\t',header=True,index=False, \
				quoting=False)

# No need to make any rarefaction plots if there's only one cell
if(len(cells) == 1):
	sys.stderr.write('Finished!\n')
	sys.exit()

def plot_and_save(xvals, yvals, title, xlab, ylab, save_path):

	f = plt.figure()
	ax = plt.subplot(1,1,1)

	xlabels = xvals.copy().astype('int').astype('str')
	xvals = xvals.astype('float')
	xvals[0] += 0.1; xvals[-1] -= 0.1

	ax.plot(xvals, yvals, linestyle='-', linewidth=3, \
			marker='o', markersize=6, markeredgecolor='none')

	ax.set_xticks(xvals)
	ax.set_xticklabels(xlabels)

	ax.set_title(title, fontsize=16)
	ax.set_xlabel(xlab, fontsize=14)
	ax.set_ylabel(ylab, fontsize=14)

	ax.set_axis_bgcolor('#F8ECE0')
	f.savefig(save_path, bbox_inches='tight',transparent=False, \
			frameon=True,edgecolor='none',dpi=250)
	f.clear()


## Plot cumulative loci
lib_cuffcmp = out_path+'/lib_cuffcmp'
os.makedirs(lib_cuffcmp)
cell_gtfs = ' '.join([out_dir+'/'+x+'/SMRTanalysis/align/final.gtf' for x in cells])
organism = lib_table.ix[0,'Organism'].lower()
ref_gtf = human[2] if(organism == 'human') else mouse[2]
cuffcmp_cmd = '%s/cuffcompare -r %s -T -R -o %s/cuffcmp %s' % (args.cufflinks,\
			 ref_gtf, lib_cuffcmp, cell_gtfs)
run_cmd(cuffcmp_cmd)

cols = ['t_id','l_id','ref_id','class_code']
for i,cell in enumerate(cells):
	cols += ['Cell'+str(i)]
cuffcmp_all = pandas.read_table(lib_cuffcmp+'/cuffcmp.tracking', names=cols)

in_cell = cuffcmp_all.ix[:,4:].apply(lambda x: x != '-', axis=0)

cumulative_loci = []
for i,cell in enumerate(cells):
	#cumulative_loci.append(len(cuffcmp_all.loc[in_cell.ix[:,range(0, i+1)].\
		#apply(lambda x: any(x), axis=1),'l_id'].unique()))
	cumulative_loci.append(len(cuffcmp_all.loc[(cuffcmp_all.class_code.isin(\
		['=','j','c'])) & (cuffcmp_all.ref_id != '-') & \
		(in_cell.ix[:,range(0, i+1)].apply(lambda x: any(x), axis=1)),\
		'ref_id'].apply(lambda x: x.partition('|')[0]).unique()))
xvals = np.arange(0,len(cells)+1)
yvals = np.array([0]+cumulative_loci)
loci_path = os.path.join(out_path,args.loci_rarefaction)

loci_path = args.loci_rarefaction+".png"
plot_and_save(xvals, yvals, 'Cumulative unique loci', '#Cells', '#Loci',\
			        loci_path)
os.rename(loci_path,args.loci_rarefaction)

## Plot cumulative junctions

juncs = []
cols = ['Chr','Start','End','Strand']
for i,cell in enumerate(cells):
	df = pandas.read_table(out_dir+'/'+cell+'/SMRTanalysis/align/junctions.bed', names=\
		cols)
	df['Cell'] = i
	juncs.append(df)
juncs = pandas.concat(juncs)

cumulative_uniq_juncs = []
for i in range(0, len(cells)):
	cumulative_uniq_juncs.append(juncs.loc[juncs.Cell.isin(range(0,i+1)),:].\
						drop_duplicates(subset = cols).shape[0])

xvals = np.arange(0,len(cells)+1)
yvals = np.array([0]+cumulative_uniq_juncs)
junction_path = os.path.join(out_path,args.junctions_rarefaction)


plot_and_save(xvals, yvals, 'Cumulative unique junctions', '#Cells', '#Juncs',\
			        args.junctions_rarefaction + ".png")
os.rename(args.junctions_rarefaction + ".png", args.junctions_rarefaction)

sys.stderr.write('Finished!\n')
