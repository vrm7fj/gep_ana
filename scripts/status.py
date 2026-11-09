#!/usr/bin/env python3
# usage: python3 status.py [logdir] [refresh_seconds] [config.json]
#   default logdir ../out/runs/logs, refresh 5 s.
#   With the config, runs that have not started yet are counted as "queued".
import glob, json, os, re, sys, time

logdir = sys.argv[1] if len(sys.argv) > 1 else "../out/runs/logs"
every  = float(sys.argv[2]) if len(sys.argv) > 2 else 5
runs   = [str(r) for r in json.load(open(sys.argv[3]))["runnums"]] if len(sys.argv) > 3 else None
prog   = re.compile(r"file (\d+)/(\d+).*FPP ok (\d+)")
errors = ("Break", "Error", "Nothing to analyze")

while True:
    n = dict(done=0, running=0, FAILED=0, queued=0)
    lines = []
    logs = {os.path.basename(p)[4:-4]: p for p in glob.glob(os.path.join(logdir, "run_*.log"))}
    for run in (runs if runs else sorted(logs)):
        if run not in logs:
            n["queued"] += 1
            continue
        text  = open(logs[run], errors="replace").read()
        m     = prog.findall(text)
        f, nf, fpp = m[-1] if m else ("-", "-", "-")
        final = re.search(r"(\d+) elastic with a good FPP", text)
        ex    = re.search(r"^EXIT (\d+)", text, re.M)       # written by run_all.sh when root exits
        if "Output written to" in text:
            state = "done"; f = nf
            if final: fpp = final.group(1)
        elif any(s in text for s in errors) or ex:          # error message, or exited without output
            state = "FAILED"
        else:
            state = "running"
        n[state] += 1
        info = f"file {f}/{nf}  FPP ok {fpp}"
        if not m and state != "done":                        # no progress line yet: show where it is
            last = [l for l in text.splitlines() if l.strip() and not l.startswith("EXIT ")]
            flagged = [l for l in last if "**!**" in l]       # the macro marks its own errors with **!**
            info = ((flagged or last)[-1].strip()[:90] if last else "(log empty)")
            if ex: info += f"  [exit {ex.group(1)}]"
        lines.append(f"{run:>6}  {state:<8} {info}")
    os.system("clear")
    total = sum(n.values())
    print(f"{time.strftime('%H:%M:%S')}   done {n['done']}   running {n['running']}   "
          f"failed {n['FAILED']}   queued {n['queued']}   (of {total})\n")
    print("\n".join(lines) or "no logs yet")
    if total and n["running"] == 0 and n["queued"] == 0:
        break
    time.sleep(every)
