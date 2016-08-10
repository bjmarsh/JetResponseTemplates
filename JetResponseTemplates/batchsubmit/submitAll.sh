
for FILE in `ls config_files/condor*.cmd`; do
    condor_submit $FILE
done