#!/bin/bash
set -e
VAR_NAME="${1:-ptl3}"
BATCH_INDEX="$2"
if [ -z "$BATCH_INDEX" ]; then
  echo "Usage: run_stage2_batch.sh <varName> <batch_index>"
  exit 1
fi
cd /afs/cern.ch/user/b/brunetd/CMSSW_14_1_6/src
cmsenv
cd ZplusXpython
exec python3 main_estimateZX_ntuples.py --var "$VAR_NAME" --condor-batch "$BATCH_INDEX"
