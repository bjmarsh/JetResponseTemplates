import glob
import os
import subprocess
import socket
import time
import sys

tag = "94x_HT_v6"

submitJobs = True  # set to False if you've already submitted jobs and just want to monitor, or if you just want to check for missing files
sweeprootExisting = False # set to False to skip sweeprooting on existing files. For if you just want to check for new cms4/resubmit
removeLogs = False # remove all condor log files after successful completion

def getCompletedFiles():
    global tag, samples
    files = set()
    for s in samples:
        files.update(set(glob.glob("/hadoop/cms/store/user/bemarsh/JRTbabies/{0}/{1}/*.root".format(tag,s))))
    return files

def sweepRoot(files):
    global samples, tag

    # cmd = "sweeproot/sweepRoot -b -o mt2 -t mt2"
    cmd = "sweeproot/sweepRoot -b -o Events"
    for fname in files:
        cmd += " "+fname
        sweeprooted.add(fname)

    # print '\n'+cmd
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    out, err = p.communicate()
    # print out
    for line in out.split('\n'):
        if line.startswith("BAD FILE"):
            print '\n'+line
            sweeprooted.remove(line.split(":")[1].strip())
            subprocess.call("hadoop fs -rm -r "+line.split(":")[1].replace("/hadoop",""), shell=True)
    

def deleteLogFiles():
    global tag
    subprocess.call("rm -r /data/tmp/bemarsh/condor_job_logs/{0}*/*.out".format(tag), shell=True)
    subprocess.call("rm -r /data/tmp/bemarsh/condor_job_logs/{0}*/*.err".format(tag), shell=True)
    subprocess.call("rm -r /data/tmp/bemarsh/condor_submit_logs/{0}*/*.log".format(tag), shell=True)
    
    
def watchCondorJobs(sweeprootWhileWaiting=False):
    global sweeprooted, haveDeletedLogs
    print "* Waiting and watching..."
    NJOBS = 1
    while NJOBS > 0:
        if "uaf-10" in socket.gethostname():
            out = subprocess.check_output(["condor_q", "bemarsh"]).split('\n')[-4].split(":")[1].strip()
            NJOBS = int(out.split()[0])
            nidle = int(out.split(";")[1].split(",")[2].strip().split()[0])
            nrunning = int(out.split(";")[1].split(",")[3].strip().split()[0])
            nheld = int(out.split(";")[1].split(",")[4].strip().split()[0])
        else:
            out = subprocess.check_output(["condor_q", "bemarsh"]).split('\n')[-2]
            NJOBS = int(out.split()[0])
            nidle = int(out.split(";")[1].split(",")[2].strip().split()[0])
            nrunning = int(out.split(";")[1].split(",")[3].strip().split()[0])
            nheld = int(out.split(";")[1].split(",")[4].strip().split()[0])
        print "  {0} jobs left. {1} idle, {2} running, {3} held".format(NJOBS, nidle, nrunning, nheld)
        if nheld > 0:
            print "* Releasing held jobs..."
            os.system("condor_release bemarsh")
        if NJOBS > 0:
            print "  Checking again in 5 min..."
            curTime = time.time()
            if sweeprootWhileWaiting:
                newFiles = getCompletedFiles().difference(sweeprooted)
                if len(newFiles) > 0:
                    print "  Sweeprooting {0} files while waiting...".format(len(newFiles)),
                    sweepRoot(newFiles)
                    if removeLogs:
                        deleteLogFiles()
                    print "Done!"
            remainingTime = 5*60 - (time.time()-curTime)
            # remainingTime = 10 - (time.time()-curTime)
            if remainingTime > 0:
                time.sleep(remainingTime)
    print "* condor jobs finished!"


        
print "* Running babymaker for tag \"{0}\"".format(tag)
samples = []
for fname in glob.glob("config_files/{0}/*.cmd".format(tag)):
    fname = fname.split("/")[-1]
    samp = fname[7:].split(".")[0]
    samples.append(samp)

print "* Found {0} samples:".format(len(samples))
for s in samples:
    print "   ", s

if not sweeprootExisting:
    sweeprooted = getCompletedFiles().copy()
else:
    sweeprooted = set()

if submitJobs:
    print "* Submitting condor jobs"
    for s in samples:
        print s
        subprocess.call("condor_submit config_files/{0}/condor_{1}.cmd".format(tag, s), shell=True)


NresubmitFiles = 1
while NresubmitFiles > 0:
    watchCondorJobs(sweeprootWhileWaiting=True)

    newFiles = getCompletedFiles().difference(sweeprooted)
    if len(newFiles)>0:
        print "* sweeprooting remaining files..."
        sweepRoot(newFiles)
                
    subprocess.call("mkdir -p config_files_resubmit/", shell=True)
    if os.listdir("config_files_resubmit/") != []:
        subprocess.call("rm config_files_resubmit/*", shell=True)
    subprocess.call("./checkAllConfig.sh config_files/{0}".format(tag), shell=True)

    resubmitFiles = glob.glob("config_files_resubmit/*.cmd".format(tag))
    NresubmitFiles = len(resubmitFiles)
    if NresubmitFiles > 0:
        for f in resubmitFiles:
            subprocess.call("condor_submit "+f+" > /dev/null", shell=True)

print "* WOO! All files are good!"
