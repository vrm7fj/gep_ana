#!/bin/bash
# usage: run_one.sh CFG OUTDIR RUN   (called by run_all.sh, from the scripts folder)
# Prints the job's output to the terminal with a [run] prefix, and also writes it to the log.
CFG=$1; OUT=$2; RUN=$3
LOG=$OUT/logs/run_$RUN.log

root -l -b -q "gep_elastic_ana.cpp(\"$CFG\",\"$OUT/gep5_$RUN\",$RUN)" 2>&1 | tee "$LOG" | sed -u "s/^/[$RUN] /"
rc=${PIPESTATUS[0]}
echo "EXIT $rc" >> "$LOG"
echo "[$RUN] ---- finished, exit $rc"
