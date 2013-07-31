"""
run_impute2.py
author: Alicia Martin (armartin@stanford.edu)

Helper script used to schedule Impute2 jobs. This script chunks the reference set into multiple
chromosome sections and schedule Impute2 jobs for each chunk in parallel.
"""

#using reference haplotypes, run imputation on study haplotypes. Requires phasing to be done already.
#also need iteration number so files don't overlap
#can feed in test_hap, ref_hap, ref_legend, and map with wrong chromosome and it will be replaced by chr option

from optparse import  OptionParser
import sys
import os
import subprocess
import commands
import numpy as np

#from lib import common

USAGE = """
run_imputation.py   --test_hap <study haplotypes from shapeit2>
                    --ref_hap <reference 1 haplotypes from shapeit2>
                    --ref_legend <reference 1 legend from shapeit2>
                    --map <recombination map>
                    --out <output file prefix>
                    --job <job id to hold on>
                    --chr
                    --first_last <first and last line of ref vcf>
                    --wd <working directory>
            OPTIONAL
                    --log <log prefix>
                    --int_size <interval size>
                    --impute2_bin <impute2 binary path>
"""
parser = OptionParser(USAGE)
parser.add_option('--first_last')
parser.add_option('--test_hap')
parser.add_option('--ref_hap')
parser.add_option('--ref_legend')
parser.add_option('--out')
parser.add_option('--wd')
parser.add_option('--chr')
parser.add_option('--job')
parser.add_option('--log')
parser.add_option('--map')
parser.add_option('--int_size',default=5000000)
parser.add_option('--impute2_bin')
parser.add_option('--qsub_job_id')

(options,args)=parser.parse_args()

#set global variables
first_last = options.first_last
test_hap = options.test_hap
ref_hap = options.ref_hap
ref_legend = options.ref_legend
wd = options.wd
out = options.out
map = options.map
int_size = options.int_size
impute2_bin = options.impute2_bin

chrom_dict = {}
first_last = open(first_last)
first_last.readline()
for line in first_last:
    myLine = line.strip().split(',')
    chrom_dict[myLine[0].lstrip('chr')] = myLine[1]

first_pos = float(1)
last_pos = float(chrom_dict[options.chr])

even_space = (last_pos - first_pos)/float(int_size)
intervals = np.linspace(first_pos-1, last_pos, even_space)
intervals2 = intervals + 1

# Helper function to run commands, handle return values and print to log file
log = open(os.path.join(options.log, os.path.basename(out) + '.log'), 'a')
def write_log(cmd, logFile, info, command):
    date=exit_status, stdout = commands.getstatusoutput('date')
    logFile.write(info + '\n')
    logFile.write(date[1] + '\n')
    logFile.write(cmd + '\n\n')
    print command #use this to get job ids
    job = command.split()[2]
    return job #if starting at a random step, need to ignore job

#########################Run impute2 for 5 Mb regions over chr 22######################
def impute2(intervals, intervals2):
    all_jobs = []
    all_iters = []
    all_info = []
    for i in range(len(intervals)-1):
        cmd = ('qsub -wd %s -hold_jid %s -w e -b y -N impute2 -e %s -o %s %s -known_haps_g %s -m %s -int %s %s -h %s -l %s -o %s'
               % (options.wd, options.job, options.log, options.log, impute2_bin, test_hap, map, str(int(intervals2[i])), str(int(intervals[i+1])), ref_hap, ref_legend, out + '_iter' + str(i)))
        exit_status, stdout = commands.getstatusoutput(cmd)
        all_jobs.append(write_log(cmd, log, 'Step 3e, Iteration ' + str(i + 1), stdout))
        all_iters.append(out + '_iter' + str(i))
        all_info.append(out + '_iter' + str(i) + '_info')
    
    cmd = ('''qsub -wd %s -hold_jid %s -w e -e %s -o %s -b y -N cat_output "cat %s | gzip > %s.gz"'''
           % (options.wd, ','.join(all_jobs), options.log, options.log, ' '.join(all_iters), out))
    exit_status, stdout = commands.getstatusoutput(cmd)
    job1 = write_log(cmd, log, 'Step 3e, Iteration ' + str(i + 1), stdout)
    
    cmd = ('''qsub -wd %s -hold_jid %s -w e -e %s -o %s -b y -N cat_info "tail -n +2 %s | sed '/^$/d' | grep -v == > %s.gz.info"'''
           % (options.wd, ','.join(all_jobs), options.log, options.log, ' '.join(all_info), out))
    exit_status, stdout = commands.getstatusoutput(cmd)
    job2 = write_log(cmd, log, 'Step 3e, Iteration ' + str(i + 1), stdout)

impute2(intervals, intervals2)
