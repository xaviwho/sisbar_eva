#!/bin/bash
# Orchestration helper used in the cloud workspace: waits until the three replicates'
# RQ1 balanced eps=0.05 potentials exist, then runs 07_rq2_forecast.py for one cell set and
# transition, serialised with flock so that only one RQ2 job runs at a time.
# usage: wait_rq2.sh CELLSET TRANSITION_TAG TRANSITION  -> waits for all 3 reps' bal_0.05 (done) then runs RQ2
cd /home/claude/sisbar_eva
while true; do
  ok=1
  for r in 1 2 3; do f=data/interim/rq1/${1}_${2}_r${r}_bal_0.05.npz
    python3 -c "import numpy as np,sys; sys.exit(0 if bool(np.load('$f')['done']) else 1)" 2>/dev/null || ok=0; done
  [ $ok = 1 ] && break
  sleep 60
done
flock data/interim/rq1/logs/rq2.lock python3 -W ignore scripts/07_rq2_forecast.py run --cellset $1 --transition "$3" > data/interim/rq1/logs/rq2_${1}_${2}.log 2>&1
echo "exit $? $(date +%T)" >> data/interim/rq1/logs/rq2_${1}_${2}.log
