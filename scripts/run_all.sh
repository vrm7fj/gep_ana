#!/bin/bash
# usage: ./run_all.sh ../config/gep5/conf_gep5.json [ncores]
# Runs every run in "runnums", ncores at a time, printing all output to the terminal.
# Ctrl-C stops everything (running and queued). Logs still go to ../out/runs/logs.
CFG=$1
NCORES=${2:-5}
OUT=../out/runs
mkdir -p $OUT/logs

# run list = "runnums" from the config
RUNS=$(python3 -c "import json,sys; print(' '.join(map(str, json.load(open('$CFG'))['runnums'])))")

# remove logs of a previous attempt
rm -f $OUT/logs/run_*.log

# Ctrl-C: kill this whole process group (xargs, every root job, tee/sed)
trap 'echo; echo "Ctrl-C: killing all jobs"; trap - INT TERM; kill -KILL 0' INT TERM

echo "$RUNS" | tr ' ' '\n' | xargs -P $NCORES -I{} "$(dirname "$0")/run_one.sh" "$CFG" "$OUT" {}
echo "all jobs finished"
