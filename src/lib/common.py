"""
common.py
author: Gerard Tse (gerardtse@gmail.com)

Common utility methods for our Python scripts. In particular
this provides helper methods for working with qsub.
"""
import commands
import errno
import os
import shutil
import stat
import sys
import time
import tempfile

class Config(object):
    SCRIPT_LOG = None
    LOGS = None
    QSUB_ID_FILE = None
    QSUB_USE_TESTQ = None
    QSUB_RUN_LOCALLY = None

CONFIG = Config()

def initialize_config(options):
    # Read globals from options
    CONFIG.SCRIPT_LOG = options.script_log
    logFile = open(CONFIG.SCRIPT_LOG, 'a')
    logFile.write(("="*100 + "\n")*3 + "New execution at " + time.ctime() + "\n" + " ".join(sys.argv) + "\n\n")
    logFile.close()

    if hasattr(options, 'logs'):
        CONFIG.LOGS = options.logs

    if hasattr(options, 'qsub_id_file'):
        CONFIG.QSUB_ID_FILE = options.qsub_id_file
    if hasattr(options, 'qsub_use_testq'):
        CONFIG.QSUB_USE_TESTQ = int(options.qsub_use_testq)
    if hasattr(options, 'qsub_run_locally'):
        CONFIG.QSUB_RUN_LOCALLY = int(options.qsub_run_locally)

# Helper function to record commands, handle return values and print to log file
def run(info, cmd):
    start_time = time.time()
    header = "====\nInfo: %s\nStart Time: %s\nCommand: %s" % (info, time.ctime(start_time), cmd)

    logFile = open(CONFIG.SCRIPT_LOG, 'a')
    logFile.write(header + "\n")
    logFile.close()

    exit_status, stdout = commands.getstatusoutput(cmd)
    end_time = time.time()

    summary = "Status: %s\nElapsed: %s ms\n%s" % (exit_status, int((end_time - start_time) * 1000), stdout)

    logFile = open(CONFIG.SCRIPT_LOG, 'a')
    logFile.write(summary + "\n\n")
    logFile.close()

    if exit_status:
        print header
        print summary
        raise RuntimeError("Command failed. See above for info")

    print stdout
    return stdout

def qsub(job_name, comment, command, wait_jobs = [], qsub_args = [], is_script = False, vmem = None, logdir = None, queue = None):
    if CONFIG.QSUB_RUN_LOCALLY:
        print "TEST ONLY! Running locally: %s\n%s" % (comment, command)
        run(comment, command)
        return None
    
    logdir = logdir or CONFIG.LOGS
    assert logdir is not None

    wait_jobs = [jid for jid in wait_jobs if jid is not None]

    cmd = 'qsub -wd %(logdir)s -w e -V -N %(job_name)s -e %(logdir)s -o %(logdir)s' % { 'logdir' : logdir, 'job_name' : job_name }
    if wait_jobs:
        cmd += " -hold_jid " + ",".join(wait_jobs)
    if not is_script:
        cmd += " -b y"
    if vmem:
        cmd += " -l h_vmem=%dg -R y" % vmem
    if CONFIG.QSUB_USE_TESTQ:
        print "Warning: Command running on test queue"
        cmd += " -l testq=1"
    elif queue:
        cmd += " -q " + queue
    if qsub_args:
        cmd += " " + qsub_args
    cmd += " " + command
    output = run(comment, cmd)
    jobid = output.split()[2]

    if CONFIG.QSUB_ID_FILE:
        f = open(CONFIG.QSUB_ID_FILE, 'a')
        f.write("%s\n" % jobid)
        f.close()
    
    return jobid

def qsub_shell_commands(job_name, comment, command, wait_jobs = [], vmem = None, logdir = None, queue = None):
    fd, fn = tempfile.mkstemp(suffix=".sh", prefix="qsub_script.")
    f = os.fdopen(fd, 'w')
    f.write("#!/bin/bash\n")
    f.write(command)
    f.close()
    os.chmod(fn, stat.S_IXUSR | stat.S_IRUSR)
    result = qsub(job_name, comment, fn, wait_jobs, is_script = True, vmem = vmem, logdir = logdir, queue = queue)
    os.remove(fn)
    return result

def mkdirp(dirname, empty = False):
    if empty:
        try:
             shutil.rmtree(dirname)
        except OSError, e:
             if e.errno == errno.ENOENT:
                pass
             else: raise

    # Make sure directory is writable
    try:
        os.makedirs(dirname)
    except OSError, e:
        if e.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else: raise

# Common setup work
def prepend_env(var, path):
    if os.environ.has_key(var):
        os.environ[var] = path + ":" + os.environ[var]
    else:
        os.environ[var] = path

def prepare_env():
    # For new libz
    prepend_env(
        'LD_LIBRARY_PATH', "/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/tools/lib")

    # For patched SLRP
    prepend_env(
        'PYTHONPATH',
        "/srv/gs1/projects/bustamante/armartin_projects/rare_imputation/tools/lib/python2.7/site-packages")

# Run this on module load time.
prepare_env()
