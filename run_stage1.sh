#!/bin/bash
set -e
cd /afs/cern.ch/user/b/brunetd/CMSSW_14_1_6/src
cmsenv
cd ZplusXpython
exec python3 main_FR_CR.py
