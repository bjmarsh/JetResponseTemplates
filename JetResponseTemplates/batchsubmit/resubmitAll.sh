
for FILE in `ls config_files/condor*resubmit.cmd`; do
    condor_submit $FILE
done