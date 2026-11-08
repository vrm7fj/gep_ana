#!/usr/bin/env python3
# usage: python3 status.py [logdir] [refresh_seconds]      (default: ../out/runs/logs, 5 s)
import glob, os, re, sys, time

logdir = sys.argv[1] if len(sys.argv) > 1 else "../out/runs/logs"
every  = float(sys.argv[2]) if len(sys.argv) > 2 else 5
prog   = re.compile(r"file (\d+)/(\d+).*FPP ok (\d+)")

while True:
    lines, ndone, nfail, nrun = [], 0, 0, 0
    for log in sorted(glob.glob(os.path.join(logdir, "run_*.log"))):
        run  = os.path.basename(log)[4:-4]
        text = open(log, errors="replace").read()
        m    = prog.findall(text)
        f, n, fpp = m[-1] if m else ("-", "-", "-")
        final = re.search(r"(\d+) elastic with a good FPP", text)
        if "Output written to" in text:
            state = "done";    ndone += 1; f = n
            if final: fpp = final.group(1)
        elif any(s in text for s in ("Break", "Error", "Nothing to analyze")):
            state = "FAILED";  nfail += 1
        else:
            state = "running"; nrun += 1
        lines.append(f"{run:>6}  {state:<8} file {f}/{n}  FPP ok {fpp}")
    os.system("clear")
    print(f"{time.strftime('%H:%M:%S')}   done {ndone}   running {nrun}   failed {nfail}\n")
    print("\n".join(lines) or "no logs yet")
    if lines and nrun == 0:
        break
    time.sleep(every)
