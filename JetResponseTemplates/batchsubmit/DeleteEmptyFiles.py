import os
import glob

base_dir = "/hadoop/cms/store/user/bemarsh/JRTbabies/v2"

dirs = glob.glob(os.path.join(base_dir,"qcd_pt*"))

for dir in dirs:
    for x in glob.glob(os.path.join(dir,"*.root")):
        if 1000 < os.path.getsize(x) < 10000:
            hadoop_name = x[x.find("/cms"):]
            cmd = "hadoop fs -rm " + hadoop_name
            print cmd
            os.system(cmd)
