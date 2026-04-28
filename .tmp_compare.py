import re
from pathlib import Path
rmsh = Path('crates/py/examples/out_all_dimensions_step_viewer_safe/all_dimensions_demo_viewer_safe.step')
gmsh = Path('crates/py/examples/out_all_dimensions_step_gmsh/all_dimensions_demo_gmsh.step')
pat = re.compile(r'^\s*#\d+\s*=\s*([A-Z0-9_]+)\s*\(', re.I)
def counts(p):
    c={}
    for line in p.read_text(encoding='utf-8', errors='ignore').splitlines():
        m=pat.match(line)
        if m:
            t=m.group(1).upper(); c[t]=c.get(t,0)+1
    return c
rc=counts(rmsh); gc=counts(gmsh)
keys=sorted(set(rc)|set(gc))
mm=[k for k in keys if rc.get(k,0)!=gc.get(k,0)]
print(f'RMSH_FILE={rmsh}')
print(f'GMSH_FILE={gmsh}')
print(f'DIFF_TYPE_COUNT={len(mm)}')
for k in mm[:40]:
    print(f'{k}: rmsh={rc.get(k,0)} gmsh={gc.get(k,0)} delta={rc.get(k,0)-gc.get(k,0)}')
