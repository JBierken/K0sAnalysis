#####################################
# resubmit failed skimming jobs     #
#####################################

# Finds failed skimming jobs from their condor log files and resubmits them
# with the job template of skim_submit.py (local scratch + checked stage-out).
#
# - Only logs modified after --since are considered (older submission rounds are ignored).
# - Jobs that are still in the condor queue (idle, running or held) are never touched.
# - Current-format jobs (set -e, '# skimjob-*' headers) count as failed if condor's
#   _log_ file shows no normal termination with return value 0. Their stderr is not
#   scanned: a stage-out attempt can crash and be retried successfully.
# - Old-format jobs always exited 0, so for them a job counts as failed if its _err_ log
#   has unmatched ###starting###/###done### tags or contains a known error string
#   (see tools/jobcheck.py).
# - The inputs and output of a failed job are read from its executable. Both the
#   current format ('# skimjob-*' header lines) and the old format (per-file skims
#   written next to the output, then hadd + rm) are understood.
# - Before resubmitting, the broken output file (and, for old-format jobs, leftover
#   per-file skims) are deleted from the output directory.
# - A failed job with more than --nfilesperjob inputs is split; its output
#   skimmed_<N>.root is replaced by skimmed_<N>-<k>.root.
# - Resubmitted executables get a '# skimjob-superseded' line, so running this
#   script again does not resubmit the same failure twice.
#
# Run from the skimming directory, e.g.:
#   python3 skim_resubmit.py --since 2026-11-23 --dryrun
#   python3 skim_resubmit.py --since 2026-11-23 -n 10

import os
import sys
import argparse
import datetime
import re
import subprocess
from six.moves import input
sys.path.append('../tools')
import condortools as ct
import jobcheck as jc
import listtools as lt
from skim_submit import CMSSW, make_job_commands, run_local

SUPERSEDED_TAG = '# skimjob-superseded'


def get_active_jobs():
  ### return the set of 'cluster.proc' ids of this user's jobs in the condor queue
  # (raises if condor_q fails, so that running jobs are never resubmitted by accident)
  out = subprocess.check_output(['condor_q', '-af', 'ClusterId', 'ProcId'])
  out = out.decode() if isinstance(out, bytes) else out
  return set('.'.join(line.split()) for line in out.splitlines() if line.strip())


def parse_executable(shfile, skimmer):
  ### read inputs, output, veto tagger and old temporary files from a job executable.
  # returns None if the executable cannot be understood.
  with open(shfile) as f: lines = [l.rstrip('\n') for l in f]
  info = {'shfile': shfile, 'inputs': [], 'output': None, 'vetotagger': None,
          'tempfiles': [], 'superseded': SUPERSEDED_TAG in lines}
  for line in lines:
    if line.startswith('# skimjob-input: '): info['inputs'].append(line.split(': ',1)[1])
    elif line.startswith('# skimjob-output: '): info['output'] = line.split(': ',1)[1]
    elif line.startswith('# skimjob-vetotagger: '): info['vetotagger'] = line.split(': ',1)[1]
  if info['output'] is not None: return info
  # old format: './skimmer <input> <temp output> -1 [tagger]', then 'hadd -f <output> <temps>'
  for line in lines:
    parts = line.split()
    if len(parts)>=4 and parts[0]=='./'+skimmer:
      info['inputs'].append(parts[1])
      info['tempfiles'].append(parts[2])
      if len(parts)>=5: info['vetotagger'] = parts[4]
    elif len(parts)>=3 and parts[0]=='hadd':
      info['output'] = parts[2] if parts[1]=='-f' else parts[1]
  if info['output'] is None or len(info['inputs'])==0: return None
  return info


def is_new_format(shfile):
  ### True if the executable was made by skim_submit.make_job_commands
  with open(shfile) as f:
    for line in f:
      if line.startswith('# skimjob-output: '): return True
  return False


def get_return_value(logfile):
  ### return value of a job from its condor _log_ file, or None if it did not terminate normally
  if not os.path.exists(logfile): return None
  with open(logfile) as f: text = f.read()
  match = re.findall(r'Normal termination \(return value (\d+)\)', text)
  return int(match[-1]) if match else None


def remove(path, dryrun):
  if not os.path.exists(path): return
  print('  removing {}'.format(path))
  if not dryrun: os.remove(path)


if __name__=='__main__':

  parser = argparse.ArgumentParser( description = 'Resubmit failed skimming jobs' )
  parser.add_argument('--since', required=True,
                      help='only consider logs modified after this date (YYYY-MM-DD or "YYYY-MM-DD HH:MM")')
  parser.add_argument('--logdir', default=CMSSW+'/src/K0sAnalysis/log_automatic_jobs')
  parser.add_argument('-s', '--skimmer', default='skim_ztomumu')
  parser.add_argument('-n', '--nfilesperjob', default=10, type=int)
  parser.add_argument('--runmode', default='condor', choices=['local', 'condor'])
  parser.add_argument('--dryrun', action='store_true', help='only print what would be done')
  args = parser.parse_args()

  # checks
  if not os.path.exists('./'+args.skimmer):
    raise Exception('ERROR: skimmer ./{} not found; run this script from the skimming directory'.format(args.skimmer))
  fmt = '%Y-%m-%d %H:%M' if ':' in args.since else '%Y-%m-%d'
  since = datetime.datetime.strptime(args.since, fmt)

  # jobs still in the queue
  active = get_active_jobs()
  print('{} of your jobs are still in the condor queue; they will be skipped.'.format(len(active)))

  # collect error logs of this submission round
  efiles = []
  for fname in os.listdir(args.logdir):
    if not (fname.startswith('cjob_skim_err_')): continue
    path = os.path.join(args.logdir, fname)
    if datetime.datetime.fromtimestamp(os.path.getmtime(path)) < since: continue
    efiles.append(path)
  print('found {} error logs modified after {}.'.format(len(efiles), since))

  # latest attempt per executable: (cluster, proc, failed)
  attempts = {}
  nactive = 0
  nunreadable = 0
  for efile in efiles:
    cluster, proc = os.path.basename(efile).split('_')[3:5]
    if '{}.{}'.format(cluster, proc) in active:
      nactive += 1
      continue
    ofile = efile.replace('_err_', '_out_')
    shfile = None
    if os.path.exists(ofile):
      with open(ofile) as f: firstline = f.readline()
      if firstline.startswith('###exename###: '):
        shfile = firstline.replace('###exename###: ', '').strip()
    if shfile is None or not os.path.exists(shfile):
      print('WARNING: cannot find the executable for {}; skipping it.'.format(efile))
      nunreadable += 1
      continue
    if is_new_format(shfile):
      failed = (get_return_value(efile.replace('_err_', '_log_'))!=0)
    else:
      failed = (jc.check_start_done(efile, verbose=False) or jc.check_error_content(efile, verbose=False))
    key = (int(cluster), int(proc))
    if shfile not in attempts or key > attempts[shfile][0]:
      attempts[shfile] = (key, failed)

  # failed jobs to resubmit
  failedjobs = []
  for shfile in sorted(attempts):
    if not attempts[shfile][1]: continue
    info = parse_executable(shfile, args.skimmer)
    if info is None:
      print('WARNING: cannot parse executable {}; skipping it.'.format(shfile))
      continue
    if info['superseded']: continue
    failedjobs.append(info)

  # make the new jobs
  allcmds = []
  for info in failedjobs:
    parts = lt.split(info['inputs'], args.nfilesperjob)
    if len(parts)==1: outputs = [info['output']]
    else:
      base = info['output'][:-len('.root')]
      outputs = ['{}-{}.root'.format(base, k) for k in range(len(parts))]
    info['newoutputs'] = outputs
    for part, output in zip(parts, outputs):
      allcmds.append(make_job_commands(args.skimmer, part, output, vetotagger=info['vetotagger']))

  # summary
  print('skipped {} jobs still in the queue and {} logs without executable.'.format(nactive, nunreadable))
  print('found {} failed jobs ({} input files) -> {} new jobs with up to {} files each:'.format(
        len(failedjobs), sum(len(i['inputs']) for i in failedjobs), len(allcmds), args.nfilesperjob))
  for info in failedjobs:
    print(' - {} ({} files) -> {} new job(s) for {}'.format(
          os.path.basename(info['shfile']), len(info['inputs']), len(info['newoutputs']), info['output']))
  if len(allcmds)==0: sys.exit()
  if args.dryrun:
    print('dry run: nothing was removed or submitted.')
    sys.exit()
  print('remove the broken outputs and resubmit (y/n)?')
  if input()!='y': sys.exit()

  # remove broken outputs and leftover per-file skims of old-format jobs
  for info in failedjobs:
    remove(info['output'], args.dryrun)
    for f in info['tempfiles']: remove(f, args.dryrun)

  # submit
  if args.runmode=='local':
    for cmds in allcmds:
      if run_local(cmds)!=0: print('### WARNING ###: job failed: {}'.format(cmds[-1]))
  else:
    ct.submitCommandsAsCondorJobsBatched(args.logdir+'/cjob_skim', allcmds, cmssw_version=CMSSW)

  # mark the old executables so they are not resubmitted again
  for info in failedjobs:
    with open(info['shfile'], 'a') as f: f.write(SUPERSEDED_TAG+'\n')
  print('done.')
