#########################################
# python script to merge built v0 files #
#########################################

import sys 
import os
sys.path.append('../tools')
import mergetools as mt

# this script will hadd all root files with appendix '_selected.root' in a folder
# command line arguments (in sequence):
# - (top-level) input directory

if __name__=='__main__':

    if len(sys.argv)==2:
        # find all subdirectories in provided directory
        # that contain at least one root file named '*_selected.root'
        input_directory = os.path.abspath( sys.argv[1] )
        indirs = {}
        for root,dirs,files in os.walk(input_directory):
            for thisdir in dirs:
                thisdir = os.path.join(root,thisdir)
                filelist = ([os.path.join(thisdir,f) for f in os.listdir(thisdir)
                            if f.endswith('_selected.root')])
                if len(filelist)>1: indirs[thisdir] = filelist
        for indir in indirs.keys():
            print('now running on {} ({} files)'.format(indir, len(indirs[indir])))
            outfile = os.path.join(indir, 'selected.root')
            if os.path.exists(outfile): os.system('rm '+outfile)
            mt.mergefiles(indirs[indir], outfile, removeinput=False, runjob=False)

    else:
        print('### ERROR ###: need different number of command line arguments')
        print('               usage: python merge.py <input_directory>')
