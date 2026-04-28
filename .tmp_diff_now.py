import re
from pathlib import Path
from collections import Counter

rmsh = Path("crates/py/examples/out_all_dimensions_step/all_dimensions_demo.step")
gmsh = Path("crates/py/examples/out_all_dimensions_step_gmsh/all_dimensions_demo_gmsh.step")
pat = re.compile(r"^\s*#\d+\s*=\s*([A-Z0-9_]+)\s*\(", re.I)

def counts(p):
    c = Counter()
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        m = pat.match(line)
        if m:
            c[m.group(1).upper()] += 1
    return c

rc, gc = counts(rmsh), counts(gmsh)
mm = [k for k in sorted(set(rc) | set(gc)) if rc.get(k, 0) != gc.get(k, 0)]
print(f"DIFF_TYPE_COUNT={len(mm)}")
for k in mm:
    print(f"{k}: rmsh={rc.get(k,0)} gmsh={gc.get(k,0)} delta={rc.get(k,0)-gc.get(k,0)}")
