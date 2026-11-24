#####################################################
# python script to prepare and submit skimming jobs #
#####################################################

# Runs on a sample list containing an arbitrary number of samples,
# see the samplelist directory for examples.
# Within the provided output directory, a subdirectory per sample will be created.
# Output files are named 'skimmed_<index>.root' and placed within the correct subdirectory.
# The number of output files is always one per job;
# if the number of input files per job is larger than one, the outputs are hadded into one.
# The per-file skims and the hadd run on the worker's local scratch disk;
# only the merged file is copied to the output directory (see make_job_commands).


# import external modules
import os
import sys
import argparse
import subprocess
from six.moves import input
# import framework modules
sys.path.append('../tools')
import condortools as ct
CMSSW = '/user/jbierken/CMSSW_15_0_17'
import listtools as lt
import sampletools as st


# stage-out helper used inside the job scripts: copy, read back and compare, retry up to 3 times
STAGE_OUT_FUNCTION = '''stage_out() {
  local src="$1" dst="$2" i
  for i in 1 2 3; do
    rm -f "$dst"
    if cp "$src" "$dst" && cmp -s "$src" "$dst"; then
      echo "copied $src to $dst (attempt $i)"
      return 0
    fi
    echo "WARNING: copy of $src to $dst failed or does not match (attempt $i)" >&2
    sleep 60
  done
  rm -f "$dst"
  return 1
}'''


def make_job_commands(skimmer, inputfiles, outputfile, vetotagger=None):
  ### make the bash commands for one skimming job.
  # - the '# skimjob-*' header lines record the inputs and output (read back by skim_resubmit.py).
  # - per-file skims and the hadd are done in a temporary directory on the worker's
  #   local disk; only the merged file is copied to outputfile and checked after copying.
  # - any failing command stops the job with a non-zero exit code and writes
  #   an '###error###' tag to stderr (recognized by tools/jobcheck.py).
  cmds = []
  for f in inputfiles: cmds.append('# skimjob-input: '+f)
  cmds.append('# skimjob-output: '+outputfile)
  if vetotagger is not None: cmds.append('# skimjob-vetotagger: '+vetotagger)
  cmds.append('set -e')
  cmds.append('''trap 'echo "###error###: exit code $? at line $LINENO: $BASH_COMMAND" >&2' ERR''')
  cmds.append('WORKDIR=$(mktemp -d "${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}/skim_XXXXXX")')
  cmds.append('''trap 'rm -rf "$WORKDIR"' EXIT''')
  cmds.append(STAGE_OUT_FUNCTION)
  tempfiles = []
  for i,f in enumerate(inputfiles):
    tempfile = '"$WORKDIR"/skim_{}.root'.format(i)
    cmd = './{} {} {} -1'.format(skimmer, f, tempfile)
    if vetotagger is not None: cmd += ' '+vetotagger
    cmds.append(cmd)
    tempfiles.append(tempfile)
  cmds.append('hadd -f "$WORKDIR"/merged.root '+' '.join(tempfiles))
  cmds.append('stage_out "$WORKDIR"/merged.root '+outputfile)
  return cmds


def run_local(cmds):
  ### run the commands of one job in a local bash shell; returns the exit code
  return subprocess.call(['bash', '-c', '\n'.join(cmds)])


if __name__=='__main__':
  
  # read command line arguments
  parser = argparse.ArgumentParser( description = 'Submitter for skimming jobs' )
  parser.add_argument(  '-i',       '--samplelist',                             required=True,  type=os.path.abspath)
  parser.add_argument(  '-o',       '--outputdir',                              required=True,  type=os.path.abspath)
  parser.add_argument(  '-s',       '--skimmer',        default='skim_ztomumu')
  parser.add_argument(  '-n',       '--nfilesperjob',   default=5,                              type=int)
  parser.add_argument(              '--recoverytag',    default=None)
  parser.add_argument(              '--runmode',        default='condor',       choices=['local', 'condor'])
  args = parser.parse_args()

  # check if samplelist exists
  if not os.path.exists(args.samplelist):
    raise Exception('ERROR: sample list '+args.samplelist+' does not seem to exist')

  # check if executable exists
  exe = './'+args.skimmer
  if not os.path.exists(exe):
    raise Exception('ERROR: skimmer '+exe+' does not seem to exist')

  # define help function to find all files in a sample directory
  def get_files(sampledir):
    # check if sample directory exists
    if not os.path.exists(sampledir):
        print('### WARNING ###: directory '+sampledir+' does not seem to exist;')
        print('                 continue ignoring this directory? (y/n)')
        go = input()
        if go=='y': return
        else: sys.exit()
    # check if sample directory is empty
    if len(os.listdir(sampledir))==0:
        print('### WARNING ###: directory '+sampledir+' is empty;')
        print('                 continue ignoring this directory? (y/n)')
        go = input()
        if go=='y': return
        else: sys.exit()
    # check if this directory has multiple subdirectories;
    # this would correspond to multiple processings of the same dataset,
    # in which case choose most recent one.
    if len(os.listdir(sampledir))>1:
        print('### WARNING ###: directory '+sampledir+' has multiple subdirectories;')
        print('                 choosing most recent processing of this dataset.')
        subdir = sorted(os.listdir(sampledir))[-1]
    else: subdir = os.listdir(sampledir)[0]
    thisinputdir = os.path.join(sampledir,subdir)
    thisfiles = []
    # list all root files in this directory
    print('Finding all files in sample {}.'.format(thisinputdir))
    for root,dirs,files in os.walk(thisinputdir):
        for fname in files:
            if fname[-5:]=='.root':
                thisfiles.append(os.path.join(root,fname))
    return (thisinputdir, thisfiles)

  # find input files for all samples
  inputfiles = []
  ninputfiles = 0
  ndirs = 0
  samples_info = st.readsamplelist(args.samplelist)
  print('Found {} samples in sample list.'.format(len(samples_info)))
  for i,sampleinfo in enumerate(samples_info):
    line = sampleinfo['line']
    temp = get_files(line)
    if temp is None: continue
    thisinputdir, thisfiles = temp
    print('Found {} files.'.format(len(thisfiles)))
    inputfiles.append({})
    inputfiles[-1]['inputdir'] = thisinputdir
    inputfiles[-1]['files'] = thisfiles
    inputfiles[-1]['outputdir'] = os.path.join(args.outputdir,
                                   sampleinfo['sample_name']+'_'+sampleinfo['version_name'])
    ninputfiles += len(thisfiles)
    ndirs += 1
    # find recovery task
    if args.recoverytag is not None:
      recoverydir = line+'_'+args.recoverytag
      if os.path.exists(recoverydir):
          recoverydir, recoveryfiles = get_files(recoverydir)
          print('Found {} additional files in recovery task'.format(len(recoveryfiles)))
          inputfiles.append({})
          inputfiles[-1]['inputdir'] = recoverydir
          inputfiles[-1]['files'] = recoveryfiles
          inputfiles[-1]['outputdir'] = os.path.join(args.outputdir,
            sampleinfo['sample_name']+'_'+sampleinfo['version_name']+'_'+args.recoverytag)
          ninputfiles += len(recoveryfiles)
          ndirs += 1

  # optional checks
  dochecks = False
  if dochecks:
    for inputstruct in inputfiles:
      for f in inputstruct['files']:
        if not os.path.exists(f):
          print('ERROR: file {} does not exist'.format(f))
        basename = os.path.basename(f)
        skimname, fn = basename.replace('.root','').split('_',1)
        try: fn = int(fn)
        except: print('WARNING: unexpected file name {}'.format(f))

  print('Found a total of {} input files in {} directories'.format(ninputfiles,ndirs))
  print('Summary of skimming jobs:')
  for inputstruct in inputfiles:
    print('  {} -> {}'.format(inputstruct['inputdir'], inputstruct['outputdir']))
  print('- Skimmer: '+args.skimmer)
  print('- Estimated total number of jobs: '+str(int(ninputfiles/args.nfilesperjob)))
  print('Continue (y/n)?')
  go = input()
  if not go=='y': sys.exit()

  # loop over samples
  allcmds = []
  for inputstruct in inputfiles:
    thisoutputdir = inputstruct['outputdir']
    # check if output directory exists; if so, clean it
    if os.path.exists(thisoutputdir):
        os.system('rm -r '+thisoutputdir)
    # make output directory
    os.makedirs(thisoutputdir)
    # get input files and split in partitions
    thisinputfiles = inputstruct['files']
    thisinputfiles = lt.split(thisinputfiles,args.nfilesperjob)
    # loop over partitions
    for counter, partlist in enumerate(thisinputfiles):
        outputfile = os.path.join(thisoutputdir, 'skimmed_'+str(counter)+'.root')
        allcmds.append(make_job_commands(args.skimmer, partlist, outputfile))

  # run or submit the jobs
  if args.runmode=='local':
    for cmds in allcmds:
      if run_local(cmds)!=0: print('### WARNING ###: job failed: {}'.format(cmds[-1]))
  elif args.runmode=='condor':
    store_dir = CMSSW + '/src/K0sAnalysis/log_automatic_jobs/'
    ct.submitCommandsAsCondorJobsBatched(store_dir + 'cjob_skim', allcmds, cmssw_version=CMSSW)
