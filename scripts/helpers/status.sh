# Progress monitor used in the cloud workspace while RQ1 ran (prints finished jobs,
# running PIDs, newest cached potentials, memory). Produces no result.
cd /home/claude/sisbar_eva
date +%T; echo "done: $(ls data/interim/rq1/*.done 2>/dev/null | xargs -n1 basename | tr '\n' ' ')"
pgrep -fa 06_rq1 | awk '{print $1}' | tr '\n' ' '; echo
ls --time-style=+%T -lt data/interim/rq1 | grep -E "npz|emb" | head -2 | awk '{print $6,$7}'
tail -1 data/interim/rq1/logs/laneA.log 2>/dev/null; free -m | sed -n 2p | awk '{print "mem used",$3}'
