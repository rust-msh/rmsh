import re
from pathlib import Path
rmsh = Path('crates/py/examples/out_all_dimensions_step/all_dimensions_demo_pipeline_rmsh.step')
gmsh = Path('crates/py/examples/out_all_dimensions_step/all_dimensions_demo_pipeline_gmsh_rewrite.step')
pat = re.compile(r'^\s*#\d+\s*=\s*([A-Z0-9_]+)\s*\(', re.I)
def counts(p):
    c={}
    for line in p.read_text(encoding='utf-8', errors='ignore').splitlines():
        m=pat.match(line)
        if m:
            t=m.group(1).upper(); c[t]=c.get(t,0)+1
    return c
rc=counts(rmsh); gc=counts(gmsh)
mm=[k for k in sorted(set(rc)|set(gc)) if rc.get(k,0)!=gc.get(k,0)]
print(f'DIFF_TYPE_COUNT_PIPELINE={len(mm)}')
