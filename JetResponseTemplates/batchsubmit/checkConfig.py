#!/bin/python

import sys, os
from glob import glob

### checks a condor submission config to see if all output files are present
### if not all present, creates a resubmission config with the missing ones
### assumes that all outputs from one config live in the same directory
### usage: python checkConfig.py <original_config>

if len(sys.argv) < 2:
    print 'no input file specified, exiting'
    exit

input_file = sys.argv[1]
input_lines = []
print 'RUNNING',input_file
try:
    f = open(input_file,'r')
    input_lines = f.readlines()
    f.close()
except IOError:
    print 'Could not open input file, exiting:',input_file
    exit

found_jobs = False
preamble = []
executable_command = ''
input_dir = ''
output_dir = ''
expected_output_files = []
for line in input_lines:
    # save preamble at the beginning of the file, for resubmission script
    #  all lines until we find "executable="
    line = line.replace("= ","=")
    if not found_jobs:
        if 'executable=' in line:
            found_jobs = True
        else:
            preamble.append(line)
            continue
    if found_jobs:
        if len(executable_command) == 0 and 'executable=' in line:
            executable_command = line
            continue
        if not 'arguments' in line:
            continue
        tokens = line.split()
        if len(input_dir) == 0:
            k = tokens[1].rfind('/')
            input_dir = tokens[1][:k]
        if len(output_dir) == 0:
            output_dir = tokens[2]
        output_file = tokens[0].split('=')[1]
        expected_output_files.append(output_file)

existing_output_files = glob(output_dir+'/*.root')

missing_files = []
for f in expected_output_files:
    expf = f+'.root'
    if not any(expf in s for s in existing_output_files):
        print 'missing file:',expf
        missing_files.append(f)

if len(missing_files) == 0:
    print 'no missing files out of %d, hurray!' % len(expected_output_files)
    exit
else:
    print 'missing %d/%d files, creating resubmission config' % (len(missing_files),len(expected_output_files))

    new_config_dir = input_file.split('/')[0]+'_resubmit'
    if not os.path.exists(new_config_dir):
        os.mkdir(new_config_dir)
    resubmit_filename = new_config_dir + '/' + input_file.split('/')[-1].replace('.cmd','_resubmit.cmd')

    # if missing files, make a resubmission config
    outf = open(resubmit_filename,'w')
    for line in preamble:
        outf.write(line)
    for mf in missing_files:
        outf.write('\n')
        outf.write(executable_command)
        outf.write('transfer_executable=True\n')
        outf.write('arguments=%s %s/%s.root %s\n' % (mf,input_dir,mf,output_dir))
        outf.write('queue\n')
        outf.write('\n')

