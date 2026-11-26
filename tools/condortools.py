###################################################################################
# functionality for making job description files for the condor submission system #
###################################################################################

# general use:
# two ingredients are needed for a condor job:
# - a job description (.txt) file
# - an executable
# the functions in this tool allow creating an executable bash script
# and its submission via a job description file

import os
import sys

# next index to try per base name, so that submitting many jobs in a loop
# does not re-scan all existing names for every job
_unique_next = {}

def makeUnique(fname, maxapp=100000):
    ### make a file name unique by appending a number to it,
    # e.g. test.txt -> test1.txt (in case test.txt already exists)
    if not os.path.exists(fname): return fname
    [name,ext] = os.path.splitext(fname)
    app = _unique_next.get(fname, 1)
    while app < maxapp:
        tryname = name+str(app)+ext
        if not os.path.exists(tryname):
            _unique_next[fname] = app+1
            return tryname
        app += 1
    msg = 'ERROR ###: already {} files named {} exist.'.format(maxapp, fname)
    msg += ' consider choosing more specific names, splitting in folders, etc.'
    raise Exception(msg)

def initJobScript(name, 
                  home=None,
                  cmssw_version=None,
                  proxy=None):
    ### initialize an executable bash script with:
    # - setting HOME environment variable
    #   note: use 'auto' to extract the HOME from os.environ
    # - sourcing the shared cms software
    # - setting correct cmssw release
    # - exporting proxy
    
    # parse filename
    name = os.path.splitext(name)[0]
    fname = name+'.sh'
    if os.path.exists(fname): os.system('rm {}'.format(fname))
    cwd = os.path.abspath(os.getcwd())
    # parse home
    if home=='auto': home = os.environ['HOME']
    # write script
    with open(fname,'w') as script:
        # write bash shebang
        script.write('#!/bin/bash\n')
        # write echo script name
        script.write("echo '###exename###: {}'\n".format(fname))
        # write export home
        if home is not None:
            script.write('export HOME={}\n'.format(home))
        # write sourcing of common software
        script.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
        # write setting correct cmssw release
        if cmssw_version is not None:
            script.write('cd {}\n'.format( os.path.join( cmssw_version,'src' ) ) )
            script.write('eval `scram runtime -sh`\n')
        # write export proxy
        if proxy is not None:
            script.write('export X509_USER_PROXY={}\n'.format( proxy ))
        script.write('cd {}\n'.format( cwd ) )
    # make executable
    os.system('chmod +x '+fname)
    print('initJobScript created {}'.format(fname))

def makeJobDescription(name, exe, argstring=None, 
                       stdout=None, stderr=None, log=None,
                       cpus=1, mem=1024, disk=10240, 
                       proxy=None, jobflavour=None, queue='queue'):
    ### create a single job description txt file
    # note: exe can for example be a runnable bash script
    # note: argstring is a single string containing the arguments to exe (space-separated)
    # note: for job flavour: see here: https://batchdocs.web.cern.ch/local/submit.html
    # note: queue is the (single) queue statement that ends the file; exe and argstring
    #       may refer to its item variable, e.g. exe='$(jobexe)', queue='queue jobexe from (...)'
    
    # parse arguments
    name = os.path.splitext(name)[0]
    fname = name+'.txt'
    if os.path.exists(fname): os.system('rm {}'.format(fname))
    if stdout is None: stdout = name+'_out_$(ClusterId)_$(ProcId)'
    if stderr is None: stderr = name+'_err_$(ClusterId)_$(ProcId)'
    if log is None: log = name+'_log_$(ClusterId)_$(ProcId)'
    # write file
    with open(fname,'w') as f:
        f.write('executable = {}\n'.format(exe))
        if argstring is not None: f.write('arguments = "{}"\n\n'.format(argstring))
        f.write('output = {}\n'.format(stdout))
        f.write('error = {}\n'.format(stderr))
        f.write('log = {}\n\n'.format(log))
        f.write('request_cpus = {}\n'.format(cpus))
        f.write('request_memory = {}\n'.format(mem))
        f.write('request_disk = {}\n'.format(disk))
        if proxy is not None: 
            f.write('x509userproxy = {}\n'.format(proxy))
            f.write('use_x509userproxy = true\n\n')
        #f.write('should_transfer_files = yes\n\n') 
        # (not fully sure whether to put 'yes', 'no' or omit it completely)
        if jobflavour is not None:
            f.write('+JobFlavour = "{}"\n\n'.format(jobflavour))

        f.write('{}\n'.format(queue))
    print('makeJobDescription created {}'.format(fname))

def _queueFrom(var, items):
    ### a single 'queue <var> from (...)' statement with one job per item.
    # with one variable, condor assigns each full line to it, so items may contain spaces.
    for item in items:
        if '\n' in item or item.strip()!=item or item=='' or item.startswith('#') or item.startswith(')'):
            raise Exception('ERROR: cannot put {!r} in a condor queue statement'.format(item))
    return 'queue {} from (\n{}\n)'.format(var, '\n'.join(items))

def submitCondorJob(jobDescription):
    ### submit a job description file as a condor job
    fname = os.path.splitext(jobDescription)[0]+'.txt'
    if not os.path.exists(fname):
        print('ERROR: job description file {} not found'.format(fname))
        sys.exit()
    # maybe later extend this part to account for failed submissions etc!
    os.system('condor_submit {}'.format(fname))

def submitCommandAsCondorJob(name, command, stdout=None, stderr=None, log=None,
                        cpus=1, mem=1024, disk=10240,
                        home=None,
                        proxy=None, 
                        cmssw_version=None,
                        jobflavour=None):
    ### submit a single command as a single job
    # command is a string representing a single command (executable + args)
    submitCommandsAsCondorJobs(name, [[command]], stdout=stdout, stderr=stderr, log=log,
            cpus=cpus, mem=mem, disk=disk,
            home=home,
            proxy=proxy,
            cmssw_version=cmssw_version,
            jobflavour=jobflavour)

def submitCommandsAsCondorCluster(name, commands, stdout=None, stderr=None, log=None,
                        cpus=1, mem=1024, disk=10240,
                        home=None,
                        proxy=None,
                        cmssw_version=None,
                        jobflavour=None):
    ### run several similar commands within a single cluster of jobs
    # note: each command must have the same executable and number of args, only args can differ!
    # note: commands can be a list of commands (-> a job will be submitted for each command)
    
    # parse arguments
    name = os.path.splitext(name)[0]
    shname = makeUnique(name+'.sh')
    jdname = name+'.txt'
    [exe,argstring] = commands[0].split(' ',1) # exe must be the same for all commands
    nargs = len(argstring.split(' ')) # nargs must be the same for all commands
    # first make the executable
    initJobScript(shname, home=home, cmssw_version=cmssw_version, proxy=proxy)
    with open(shname,'a') as script:
        script.write(exe)
        script.write(' $@')
        script.write('\n')
    # collect the arguments of all jobs
    argstrings = []
    for command in commands:
        [thisexe,thisargstring] = command.split(' ',1)
        thisnargs = len(thisargstring.split(' '))
        if( thisexe!=exe or thisnargs!=nargs):
            print('### ERROR ###: commands are not compatible to put in same cluster')
            return
        argstrings.append(thisargstring)
    # then make the job description, with a single queue statement over all jobs
    # (multiple queue statements per submit file are deprecated in condor 24 and removed in 25)
    makeJobDescription(name,shname,argstring='$(jobargs)',stdout=stdout,stderr=stderr,log=log,
                       cpus=cpus,mem=mem,disk=disk,proxy=proxy,
                       jobflavour=jobflavour, queue=_queueFrom('jobargs', argstrings))
    # finally submit the job
    submitCondorJob(jdname)

def submitCommandsAsCondorJob(name, commands, stdout=None, stderr=None, log=None,
                        cpus=1, mem=1024, disk=10240, 
                        home=None,
                        proxy=None,
                        cmssw_version=None,
                        jobflavour=None):
    ### submit a set of commands as a single job
    # commands is a list of strings, each string represents a single command (executable + args)
    # the commands can be anything and are not necessarily same executable or same number of args.
    submitCommandsAsCondorJobs(name, [commands], stdout=stdout, stderr=stderr, log=log,
                        cpus=cpus, mem=mem, disk=disk, 
                        home=home,
                        proxy=proxy,
                        cmssw_version=cmssw_version,
                        jobflavour=jobflavour)

def submitCommandsAsCondorJobs(name, commands, stdout=None, stderr=None, log=None,
            cpus=1, mem=1024, disk=10240,
            home=None,
            proxy=None,
            cmssw_version=None,
            jobflavour=None):
    ### submit multiple sets of commands as jobs (one job per set)
    # commands is a list of lists of strings, each string represents a single command
    # the commands can be anything and are not necessarily same executable or number of args.
    for commandset in commands:
        # parse arguments
        name = os.path.splitext(name)[0]
        shname = makeUnique(name+'.sh')
        jdname = name+'.txt'
        # first make the executable
        initJobScript(shname, home=home, cmssw_version=cmssw_version, proxy=proxy)
        with open(shname,'a') as script:
             for cmd in commandset: script.write(cmd+'\n')
        # then make the job description
        makeJobDescription(name,shname,stdout=stdout,stderr=stderr,log=log,
                            cpus=cpus,mem=mem,disk=disk,proxy=proxy,
                            jobflavour=jobflavour)
        # finally submit the job
        submitCondorJob(jdname)

def submitCommandsAsCondorJobsBatched(name, commands, stdout=None, stderr=None, log=None,
            cpus=1, mem=1024, disk=10240,
            home=None,
            proxy=None,
            cmssw_version=None,
            jobflavour=None):
    ### same as submitCommandsAsCondorJobs (one job per set of commands),
    # but all jobs go into one job description file and one condor_submit call.
    # note: all jobs end up in a single cluster (one ProcId per job).
    # returns the list of executables that were created.
    if len(commands)==0: return []
    name = os.path.splitext(name)[0]
    jdname = name+'.txt'
    shnames = []
    for commandset in commands:
        shname = makeUnique(name+'.sh')
        initJobScript(shname, home=home, cmssw_version=cmssw_version, proxy=proxy)
        with open(shname,'a') as script:
            for cmd in commandset: script.write(cmd+'\n')
        shnames.append(shname)
    # one job description with a single queue statement over all executables
    # (multiple queue statements per submit file are deprecated in condor 24 and removed in 25)
    makeJobDescription(name,'$(jobexe)',stdout=stdout,stderr=stderr,log=log,
                       cpus=cpus,mem=mem,disk=disk,proxy=proxy,
                       jobflavour=jobflavour, queue=_queueFrom('jobexe', shnames))
    submitCondorJob(jdname)
    return shnames
