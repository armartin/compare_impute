"""
generate_vcfs.py
author: Gerard Tse (gerardtse@gmail.com)

Generates test and reference sets in the format of VCF from the 1000 Genome data. This script
generates sets sample size and population based analysis. This data this script produces is
consumed by generate_impute_results.py.
"""
import errno
from optparse import OptionParser
import os
import random
import re
import stat
import sys
import tempfile

# Our own libraries
from lib import analysis
from lib import common

USAGE = """
generate_vcf.py
No flags required by default.
Please refer to source code for command-line options.
"""

DEFAULT_SCRIPT_LOG_LOCATION = analysis.SAMPLED_DATA_BASE + "/generate_vcfs_script.log"
parser = OptionParser(USAGE)
parser.add_option('--script_log', default=DEFAULT_SCRIPT_LOG_LOCATION, help='script log file')
parser.add_option('--vcftools', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/tools/bin/vcftools')
parser.add_option('--inds_data', default='/srv/gs1/projects/bustamante/reference_panels/1kG_DATA/integrated/20120317_new_phase1_integrated_genotypes_version_3/PopInds')
parser.add_option('--java', default='/srv/gs1/projects/bustamante/scg3_inhousebin/java')
parser.add_option('--split_vcf_jar', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/tools/bin/splitVCFref.jar')
parser.add_option('--markers_affy', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/affy6_markers.txt')
parser.add_option('--markers_illumina', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/illumina3_markers.txt')
parser.add_option('--markers_omni', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/omni2_5_markers.txt')
parser.add_option('--markers_exome', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/exome_chip/annotatedList.txt')
parser.add_option('--markers_affy_exome', default='/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/data/exome_chip/affy/Axiom_bbv09_content.core_50k_YRI.rfmt.txt')

(options,args) = parser.parse_args()
if options.script_log is None:
    parser.error('no script log file given')
common.initialize_config(options)

# Helper class for loading individuals and sampling them
class Reservior(object):
    def __init__(self, items):
        self.items = random.sample(list(items), len(items))

    def consume(self, count):
        count = min(count, len(self.items))
        consumed = self.items[:count]
        self.items = self.items[count:]
        return consumed

    def half(self):
        return self.consume(len(self.items) / 2)

    def all(self):
        return self.consume(len(self.items))

class Individuals(object):
    def __init__(self):
        self.load()

    def load(self):
        filenames = os.listdir(options.inds_data)
        samples = {}
        all_inds = []
        for fn in filenames:
            pop = re.sub(r'\.inds$', '', fn)
            f = open(options.inds_data + "/" + fn, 'r')
            inds = [line.strip() for line in f]
            samples[pop] = inds
            all_inds += inds
            
        self.samples = samples
        self.pops = samples.keys()
        self.all_inds = all_inds

    def create_sampling(self, include_pops=None, exclude_pops=None):
        assert include_pops is None or exclude_pops is None
        if include_pops is not None:
            for pop in include_pops:
                assert pop in self.pops
        if exclude_pops is not None:
            for pop in exclude_pops:
                assert pop in self.pops

        if include_pops is None and exclude_pops is None:
            return Reservior(self.all_inds)
        pops = None
        if include_pops:
            pops = include_pops
        elif exclude_pops:
            pops = [pop for pop in self.pops if pop not in exclude_pops]
        inds = []
        for pop in pops:
            inds += self.samples[pop]
        return Reservior(inds)

def gzinput(chrom):
    return '/srv/gs1/projects/bustamante/reference_panels/1kG_DATA/integrated/20120317_new_phase1_integrated_genotypes_version_3/orig_files/ALL.chr%d.phase1_release_v3.20101123.snps_indels_svs.genotypes.vcf.gz' % chrom

def process_sample_size_analysis(individuals):
    for item in analysis.SampleSizeAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "samplesize-%04d-%04d-ch%02d-run-%02d" % (item.testsize(), item.refsize(), chrom, iteration)
            sampling = individuals.create_sampling()
            ref_inds = sampling.consume(item.refsize())
            test_inds = sampling.consume(item.testsize())
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

def process_population_analysis(individuals):
    for item in analysis.HoldOnePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-excluded-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            test_inds = individuals.create_sampling(include_pops = [item.pop()]).all()
            ref_inds = individuals.create_sampling(exclude_pops = [item.pop()]).consume(900)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    for item in analysis.HoldNonePopAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-included-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(exclude_pops = [item.pop()])
            ref_inds = (pop_inds.half() + other_inds.all())[:900]
            test_inds = pop_inds.all()  # All the remaining, which is the other half
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    # Americas is the smallest group with 3 pops, 60, 66, 55 individuals each.
    # We leave up to 33 in the test set and 60 + 66/2 + 55 = 148 in the ref set and apply it all over.
    # For Iberian people the set has size 14. In this case all of it becomes test case, and we take
    # 148 from the other pops.
    # so we are comparing apples to apples
    for item in analysis.RegionBiasedControlAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-unbiased-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(exclude_pops = [item.pop()])
            test_inds = pop_inds.consume(30)
            ref_inds = random.sample(pop_inds.all() + other_inds.all(), 148)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

    for item in analysis.RegionBiasedAnalysis.ALL:
        for chrom, iteration in item.iterate():
            job_name_suffix = "pop-biased-%s-ch%02d-run-%02d" % (item.pop(), chrom, iteration)
            pop_inds = individuals.create_sampling(include_pops = [item.pop()])
            other_inds = individuals.create_sampling(include_pops = item.regional_pops())
            test_inds = pop_inds.consume(33)
            ref_inds = random.sample(pop_inds.all() + other_inds.all(), 148)
            create_vcfs(job_name_suffix, item, chrom, iteration, gzinput(chrom), ref_inds, test_inds)

def create_inds_file(filename, samples):
    f = open(filename, 'w')
    f.writelines(["%s\n" % sample for sample in samples])
    f.close()

def create_vcf(jobname, outfile, vcftools_args, wait_jobs = []):
    if os.path.exists(outfile):
        # No need to recreate
        print "Output already exists, skipping: %s" % jobname
        return None
    # Make a log dir if it doesn't exists
    logdir = os.path.join(os.path.dirname(outfile), 'logs')
    common.mkdirp(logdir)
    return common.qsub_shell_commands(
        jobname, "Create VCF: " + jobname,
        """
        mkdir %(jobname)s &&
        cd %(jobname)s &&
        %(vcftools)s %(vcftools_args)s --recode-to-stream | gzip -c - > %(outfile)s.inprogress &&
        cat out.log &&
        cd .. &&
        mv %(outfile)s.inprogress %(outfile)s &&
        rm -Rf %(jobname)s
        """ % {
            'jobname' : jobname,
            'vcftools' : options.vcftools,
            'vcftools_args' : vcftools_args,
            'outfile' : outfile
        },
        wait_jobs = wait_jobs,
        logdir = logdir
        )

def chunk_vcf(jobname, vcffile, wait_jobs=[], logdir=None, core_size=5000*1000, flanking_size=250*1000):
    if os.path.exists(vcffile + ".splitlog.csv"):
        # No need to recreate
        print "Output already exists, skipping: %s" % jobname
        return None
    
    logdir = os.path.join(os.path.dirname(vcffile), 'logs')
    cmd = """
      #
      # To determine the chunk size
      # First find first and last locations in the VCF
      # Then calculate the range size and divide by desired chunk size for number of chunks
      # Use the ceil of that to recalculate actual chunk size with that many chunks
      #
      MIN_AND_MAX=$(zcat %(vcf)s | grep -v -E "^#" | awk '{print $2}' | sort |  awk 'NR==1; END{print}')  &&
      MAX=$(echo "$MIN_AND_MAX" | tail -n 1)  &&
      MIN=$(echo "$MIN_AND_MAX" | head -n 1)  &&
      CHUNKSIZE=$(python -c "from math import ceil; range=$MAX-$MIN; chunks=ceil(range/float(%(core_size)s)); print int(ceil(range / chunks))") &&
      echo "Using min:$MIN max:$MAX desired_chunk_size:%(core_size)s actual_chunk_size:$CHUNKSIZE" &&
      # To work around a bug in base-pair count parsing in the jar, add period (It normally expects Kb Mb)
      %(java)s -Xmx4g -jar %(split_vcf_jar)s --vcfref %(vcf)s --coresize ${CHUNKSIZE}.b --flankingsize %(flanking_size)s.b
    """
    return common.qsub_shell_commands(
        jobname, "Chunk VCF: " + jobname,
        cmd % {
            'jobname' : jobname,
            'java': options.java,
            'split_vcf_jar' : options.split_vcf_jar,
            'vcf' : vcffile,
            'core_size' : core_size,
            'flanking_size' : flanking_size
        },
        wait_jobs = wait_jobs,
        logdir = logdir,
        vmem = 7
        )

def create_vcfs(job_name_suffix, item, chrom, iteration, src_vcf, ref_inds, test_inds):
    job_dir = item.dirname(iteration)
    common.mkdirp(job_dir)
    # Drop in a file called job_name.cfg
    f = open(os.path.join(job_dir, "job_name.cfg"), 'w')
    f.write(job_name_suffix)
    f.close()
    
    ref_vcf = item.vcf_filename(analysis.T_REF, chrom, iteration)
    test_vcf = item.vcf_filename(analysis.T_TEST_TRUTH, chrom, iteration)
    
    state = os.path.exists(ref_vcf) + os.path.exists(test_vcf)
    if state == 1:
        if os.path.exists(ref_vcf):
            os.remove(ref_vcf)
        if os.path.exists(test_vcf):
            os.remove(test_vcf)
    
    ref_inds_file = None
    test_inds_file = None
    if state != 2:
        ref_inds_file = item.inds_filename(analysis.T_REF)
        test_inds_file = item.inds_filename(analysis.T_TEST_TRUTH)
        create_inds_file(ref_inds_file, ref_inds)
        create_inds_file(test_inds_file, test_inds)
    
    make_ref_jid = create_vcf(
        "make-ref-" + job_name_suffix,
        item.vcf_filename(analysis.T_REF, chrom, iteration),
        "--gzvcf %s --keep %s --remove-indels" % (src_vcf, ref_inds_file))
    
    chunk_vcf(
        "chunk-ref-" + job_name_suffix,
        item.vcf_filename(analysis.T_REF, chrom, iteration),
        wait_jobs = [make_ref_jid])
    make_test_jid = create_vcf(
        "make-test-" + job_name_suffix,
        test_vcf,
        "--gzvcf %s --keep %s --remove-indels" % (src_vcf, test_inds_file))
    
    create_vcf(
        "make-test-affy-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_AFFY, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_affy),
        wait_jobs = [make_test_jid])
    
    create_vcf(
        "make-test-illumina-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_ILLUMINA, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_illumina),
        wait_jobs = [make_test_jid])
    
    create_vcf(
        "make-test-omni-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_OMNI, chrom, iteration),
        "--gzvcf %s --snps %s" % (test_vcf, options.markers_omni),
        wait_jobs = [make_test_jid])
    
    create_vcf(
        "make-test-exome-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_EXOME, chrom, iteration),
        "--gzvcf %s --positions %s" % (test_vcf, options.markers_exome),
        wait_jobs = [make_test_jid])
    
    create_vcf(
        "make-test-affy-exome-" + job_name_suffix,
        item.vcf_filename(analysis.T_TEST_AFFY_EXOME, chrom, iteration),
        "--gzvcf %s --positions %s" % (test_vcf, options.markers_affy_exome),
        wait_jobs = [make_test_jid])

def main():
    individuals = Individuals()
    process_sample_size_analysis(individuals)
    process_population_analysis(individuals)

if __name__ == "__main__":
    main()
