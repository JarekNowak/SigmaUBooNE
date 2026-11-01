#!/bin/bash
# Raise the running limit of the geniefix arrays when the cluster is empty, drop it back when
# others are waiting.  Every INTERVAL seconds:
#   free = sum over usable nodes of min(idle CPUs, unallocated memory / 2 GB)   (our jobs: 1 CPU, 2G)
#   if no other user has PENDING jobs:  limit = running + free - RESERVE, clamped to [BASE, MAX]
#   else:                                limit = BASE
# and sets ArrayTaskThrottle on every geniefix array still in the queue. Lowering the limit never
# kills a running job; it only stops new starts. Exits when no geniefix job is left.
# Usage: nohup ./throttle_daemon.sh > throttle.log 2>&1 &
BASE=${BASE:-100}; MAX=${MAX:-240}; RESERVE=${RESERVE:-16}; INTERVAL=${INTERVAL:-600}
EXCL='py-batch1[3-6]'          # submit_genie_hypfix.sh excludes py-batch[13-16]
MEMJOB=2000
last=""
while :; do
  arrays=$(squeue -u "$USER" -h -o '%F %j' | awk '$2 ~ /^geniefix/ {print $1}' | sort -u)
  [ -n "$arrays" ] || { echo "$(date '+%F %T') no geniefix jobs left, exiting"; exit 0; }
  running=$(squeue -u "$USER" -h -t RUNNING -o '%j' | grep -c '^geniefix')
  others=$(squeue -h -t PENDING -o '%u' | grep -vcx "$USER")
  free=0
  while read -r n st cpus mem amem; do
    [[ "$n" =~ $EXCL ]] && continue
    case "$st" in idle*|mixed*) ;; *) continue ;; esac
    idle=$(echo "$cpus" | cut -d/ -f2)
    m=$(( (mem - amem) / MEMJOB )); [ $m -lt 0 ] && m=0
    free=$(( free + (idle < m ? idle : m) ))
  done < <(sinfo -h -N -p main -o '%n %t %C %m %e' | while read -r n st c m e; do
             am=$(scontrol show node "$n" | grep -o 'AllocMem=[0-9]*' | cut -d= -f2); echo "$n $st $c $m ${am:-0}"; done)
  if [ "$others" -eq 0 ]; then lim=$(( running + free - RESERVE )); else lim=$BASE; fi
  [ $lim -lt $BASE ] && lim=$BASE; [ $lim -gt $MAX ] && lim=$MAX
  if [ "$lim" != "$last" ]; then
    for a in $arrays; do scontrol update JobId="$a" ArrayTaskThrottle="$lim"; done
    echo "$(date '+%F %T') running $running, free slots $free, others pending $others -> limit $lim ($(echo $arrays | tr '\n' ' '))"
    last=$lim
  fi
  sleep "$INTERVAL"
done
