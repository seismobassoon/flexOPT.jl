"""Compare the pure-Python reader (HDF5 export) with the Julia reader (JLD2 catalogue) for a few cases.
Run on a compute node:  srun -p ncpushort -n1 -c1 --mem=4G -t 0:30:00 python3 check_python_vs_julia.py 00001
"""
import sys
import numpy as np
import moon_rsgt as jl          # Julia wrapper (synth.jl)
import moon_rsgt_py as py       # pure Python

mid = sys.argv[1] if len(sys.argv) > 1 else py.model_ids()[0]
cases = [dict(stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=300.0),
         dict(stla=-3.0, stlo=-23.4, evla=-25.0, evlo=70.0, evdp=913.0),
         dict(stla=-3.0, stlo=-23.4, evla=10.0, evlo=-60.0, evdp=9.0)]
M = dict(mrr=3e14, mtt=-1e15, mpp=7e14, mrt=2e14, mrp=-4e14, mtp=5e14)
for c in cases:
    t1, d1, _ = jl.moon_synthetics(mid, **c, **M, comps="ZRT", quantity="velocity", launcher=[], verbose=False)
    t2, tr = py.synthetics(mid, **c, **M, comps="ZRT", quantity="velocity")
    for i, comp in enumerate("ZRT"):
        a, b = d1[0, 0, i], tr[comp]
        print(f"evdp {c['evdp']:6.1f} km  {comp}:  relative difference {np.linalg.norm(a - b) / np.linalg.norm(a):.2e}")
