## Jet Response Templates

CMSSW package for deriving jet response templates (distributions of pT(reco)/pT(gen)). Runs over AOD MC.

Check this repository out in the `src` directory of a CMSSW release (tested using 8_0_11, but will probably work with others). In the directory where this README is, you should be able to do a `cmsenv` and a `scram b` to compile everything. This will put the babymaker and the looper in your PATH.

#### Directory contents

bin/

* `JRTbabymaker.cc` - The babymaker. Uses the python config files in the `test` directory.
* `JRTlooper.cc` - Looper that runs over the babies and produces the jet response templates. Usage is `JRTlooper <input files>`.
* `JRTprintjets.cc` - Used for debugging, prints contents of the jets of a specified event.

src/

* `JRTTree.cc` - class for the handling of the TTree in the babies. Can be used in "write mode" (in the babymaker) or in "read mode" (in the looper).

interface/

* `JRTTree.h` - header for the above.

test/

* `test_cfg.py` - sample cfg file for testing the babymaker. Run with `JRTbabymaker test_cfg.py`.
* `condor_template_cfg.py` - don't touch this. Template cfg file used by the condor jobs.

jecs/

* Contains all of the JEC text files used by the babymaker.

batchsubmit/

* All of the condor submission and merging scripts. See the directory for instructions.
