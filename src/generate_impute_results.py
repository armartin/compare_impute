"""
generate_impute_results.py
author: Gerard Tse (gerardtse@gmail.com)

Provides functionalities for job management of phase and impute pipelines.
Jobs can be started, stopped, and inspected by using different parameters.
"""
import commands
import errno
import inspect
from optparse import OptionParser
import os
import re
import shutil
import time

# Our own libraries
from lib import analysis
from lib import common

USAGE = """
generate_impute_results.py
No flags required by default.
Please refer to source code for command-line options.
"""

DEFAULT_SCRIPT_LOG_LOCATION = analysis.SAMPLED_DATA_BASE + "/generate_impute_results_script.log"
CURRENT_SCRIPT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
parser = OptionParser(USAGE)
parser.add_option('--script_log', default=DEFAULT_SCRIPT_LOG_LOCATION, help='script log file')
parser.add_option('--script_dir', default=CURRENT_SCRIPT_DIR, help='script log file')
parser.add_option('--steps', default='1,2', help='steps_to_run')
parser.add_option('--types', default='a,b,c,d,e', help="analysis types to run")
parser.add_option('--arrays', default="affy", help='test arrays to process, see lib/analysis.py for valid types')
parser.add_option('--runonly', default='small', help="Can either be small, large, an execution ID, or if empty, runs them all")
parser.add_option('--dry_run', default='1', help="if set, no jobs are scheduled")
parser.add_option('--reset', default='0', help="if set, reset job if currently failed.")
parser.add_option('--reset_running', default='0', help="if set, reset job if currently incomplete, too.")
parser.add_option('--reset_completed', default='0', help="if set, reset job if currently completed, too.")

# Pass through
parser.add_option('--qsub_use_testq', default='0', help='use the test queue for scheduling')
parser.add_option('--qsub_run_locally', default='0', help='run jobs locally instead of using qsub')

(options,args) = parser.parse_args()
if options.script_log is None:
    parser.error('no script log file given')
common.initialize_config(options)

def process_sample(item, chrom, iteration):
    ref_gzvcf = item.vcf_filename(analysis.T_REF, chrom)
    if not os.path.exists(ref_gzvcf):
        print "%s: Reference VCF not found. Skipping" % item.display_name()
        return
    
    for step in options.steps.split(","):
        for test_set in options.arrays.split(","):
            test_set = "test_" + test_set
            test_gzvcf = item.vcf_filename(test_set, chrom)
            if not os.path.exists(test_gzvcf):
                print "%s: Test VCF not found for %s. Skipping this test set." % (item.display_name(), test_set)
                continue
            start_pipeline(item, chrom, iteration, step, test_set)

def run_instance_id(step, chrom, test_set):
    return "step%s-ch%02d-%s" % (step, chrom, test_set.replace("test_", ""))

def parse_instance_id(id):
    return {
        'step': id[4],
        'chrom': int(id[8:10]),
        'test_set': "test_" + id[11:]
    }
    
def start_pipeline(item, chrom, iteration, step, test_set):
    base_dir = item.dirname(iteration)
    job_base_dir = os.path.join(base_dir, "imputation", run_instance_id(step, chrom, test_set))

    job_out_dir = os.path.join(job_base_dir, "results")
    job_log_dir = os.path.join(job_base_dir, "logs")
    job_tmp_dir = os.path.join(job_base_dir, "tmp")

    job_id_file = os.path.join(job_log_dir, "job_id.txt")
    was_run = os.path.exists(job_id_file)
    if was_run:
        _, create_time = commands.getstatusoutput("find %s -type f -printf '%%T@\n' | sort | head -n 1" % job_base_dir)
        _, last_mod_time = commands.getstatusoutput("find %s -type f -printf '%%T@\n' | sort -r | head -n 1" % job_base_dir)
        create_time = time.ctime(float(create_time.strip()))
        last_mod_time = time.ctime(float(last_mod_time.strip()))

    completed_file = os.path.join(job_out_dir, "completed.flag")

    log_prefix = "%s/%s" % (item.display_name(), os.path.basename(job_base_dir))
    log_prefix = "%-52s" % log_prefix
    if os.path.exists(completed_file) and not (int(options.reset) and int(options.reset_completed)):
        print "%s: OK done at %s" % (log_prefix, last_mod_time)
        return

    running_jobs = []
    
    if was_run:
        f = open(job_id_file, 'r')
        jobs = [l.strip() for l in f]
        f.close()
        
        code, _ = commands.getstatusoutput("qstat -j %s" % ",".join(jobs))
        # qstat returns success, meaning some job in that list is still running.
        if not code:
            running_jobs = jobs

    if int(options.reset):
        if was_run:
            if running_jobs and int(options.reset_running):
                if int(options.dry_run):
                    print "%s: WILL RESET kill running jobs: %s" % (log_prefix, ",".join(running_jobs))
                else:
                    print "%s: RESET kill running jobs: %s" % (log_prefix, ",".join(running_jobs))
                    commands.getstatusoutput("qdel " + ",".join(running_jobs))
                running_jobs = []
            if not running_jobs:
                if int(options.dry_run):
                    print "%s: WILL RESET deleting data" % log_prefix
                else:
                    print "%s: RESET deleting data" % log_prefix
                    common.mkdirp(job_out_dir, empty=True)
                    common.mkdirp(job_log_dir, empty=True)
                    common.mkdirp(job_tmp_dir, empty=True)
            else:
                print "%s: RUNNING (scheduled %s)" % (log_prefix, create_time)
        # When resetting. no jobs are submitted
        else:
            print "%s: NOT RUN" % log_prefix
        return

    if int(options.dry_run):
        if was_run:
            if running_jobs:
                print "%s: RUNNING (scheduled %s)" % (log_prefix, create_time)
            else:
                print "%s: WILL RUN. Last run failed (%s)" % (log_prefix, last_mod_time)
        else:
            print "%s: WILL RUN" % log_prefix
        return
            
    if running_jobs:
        print "%s: RUNNING (scheduled %s)" % (log_prefix, create_time)
        return

    print "%s: OK scheduling" % log_prefix
    common.mkdirp(job_out_dir, empty=True)
    common.mkdirp(job_log_dir, empty=True)
    common.mkdirp(job_tmp_dir, empty=True)

    # Read job_name
    f = open(os.path.join(base_dir, 'job_name.cfg'), 'r')
    job_name = f.read().strip()
    f.close()

    # Count individuals so pipeline can get a better estimate on memory usage.
    ind_count = 0
    f = open(item.inds_filename(analysis.T_REF), 'r')
    ind_count += len(f.readlines())
    f.close()
    f = open(item.inds_filename(analysis.T_TEST_TRUTH), 'r')
    ind_count += len(f.readlines())
    f.close()
    
    command = """
      python %(script_dir)s/phase_impute_pipeline.py
        --ind_ref=%(ind_ref)s
        --ind_count=%(ind_count)d
        --ref_gzvcf=%(ref_gzvcf)s
        --test_gzvcf=%(test_gzvcf)s
        --step=%(step)s
        --logs=%(logdir)s
        --script_log=%(script_log)s
        --script_dir=%(script_dir)s
        --job_name_prefix=%(job_name)s
        --intermediate_dir=%(intermediate_dir)s
        --outdir=%(out_dir)s
        --completed_file=%(completed_file)s
        --chr=%(chrom)d
        --qsub_id_file=%(qsub_id_file)s
        --qsub_use_testq=%(qsub_use_testq)s
        --qsub_run_locally=%(qsub_run_locally)s
    """ % {
        'ind_ref': item.inds_filename(analysis.T_REF),
        'ind_count': ind_count,
        'ref_gzvcf': item.vcf_filename(analysis.T_REF, chrom),
        'test_gzvcf': item.vcf_filename(test_set, chrom),
        'step': step,
        'logdir': job_log_dir,
        'script_log': os.path.join(job_log_dir, "script.log"),
        'script_dir': options.script_dir,
        'job_name': job_name,
        'intermediate_dir': job_tmp_dir,
        'out_dir': job_out_dir,
        'completed_file': completed_file,
        'chrom': chrom,
        'qsub_id_file': job_id_file,
        'qsub_use_testq': options.qsub_use_testq,
        'qsub_run_locally': options.qsub_run_locally }
    command = re.sub(r'\s+', ' ', command).strip()
    try:
      common.run(item.display_name(), command)
    except Exception, e:
      print "%s: Pipeline command failed. Exiting" % item.display_name()
      raise e

def process_sample_size_analysis():
    for item in analysis.SampleSizeAnalysis.ALL:
        for chrom, iteration in item.iterate():
            process_sample(item, chrom, iteration)

def process_hold_one_pop_analysis():
    for item in analysis.HoldOnePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            process_sample(item, chrom, iteration)

def process_hold_none_pop_analysis():
    for item in analysis.HoldNonePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            process_sample(item, chrom, iteration)

def process_region_biased_analysis():
    for item in analysis.RegionBiasedAnalysis.ALL:
        for chrom, iteration in item.iterate():
            process_sample(item, chrom, iteration)

def process_region_biased_control_analysis():
    for item in analysis.RegionBiasedControlAnalysis.ALL:
        for chrom, iteration in item.iterate():
            process_sample(item, chrom, iteration)

def main():
    if options.runonly == "small":
        # In fast test mode, process the first sample size analysis item (small)
        process_sample(analysis.SampleSizeAnalysis.ALL[0], 22, 1)
    elif options.runonly == "large":
        # In large test mode, process the last sample size analysis item (large)
        process_sample(analysis.SampleSizeAnalysis.ALL[-1], 22, 1)
    elif options.runonly:
        run_id = options.runonly
        run_analysis_type = run_id[:run_id.rfind('/')]
        run_instance_id = run_id[run_id.rfind('/') + 1:]
        item = [ item for item in analysis.Analysis.ALL if item.path() == run_analysis_type ][0]
        parsed = parse_instance_id(run_instance_id)
        start_pipeline(item, parsed['chrom'], 1, parsed['step'], parsed['test_set'])
    else:
        # Pipeline mode. Process everything
        types = options.types.split(",")
        if 'a' in types:
            process_hold_one_pop_analysis()
        if 'b' in types:
            process_hold_none_pop_analysis()
        if 'c' in types:
            process_sample_size_analysis()
        if 'd' in types:
            process_region_biased_analysis()
        if 'e' in types:
            process_region_biased_control_analysis()
    if int(options.dry_run):
        print "Dry run mode. set --dry_run=0 for skipped items to excute"

if __name__ == "__main__":
    main()
