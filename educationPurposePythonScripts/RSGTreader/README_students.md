# Moon synthetics for any 1-D model of the catalogue (Python)

Quickest start: open `moonRSGT_tutorial_python.ipynb` (or `moonRSGT_tutorial_julia.ipynb` with
`moon_rsgt_h5.jl`), set the catalogue path in the first cell, and run all cells.

You need Python 3 with `numpy` and `h5py` (`pip install numpy h5py`), the file
`moon_rsgt_py.py`, and the catalogue folder (`rsgt_00001.h5` … one file per model,
≈ 0.6 GB each).

```python
import moon_rsgt_py as moon
moon.CATALOGUE = "/path/to/catalogue_h5"           # the folder on the disk

# one model, one event, one station -- velocity (m/s), components Z, R, T
t, tr = moon.synthetics("00817",
                        stla=-3.0, stlo=-23.4,                 # station (deg)
                        evla=20.0, evlo=15.0, evdp=900.0,      # event (deg, km)
                        mrr=1e15, mtt=-5e14, mpp=-5e14,        # moment tensor, GCMT convention, N.m
                        mrt=0.0, mrp=0.0, mtp=0.0)
import matplotlib.pyplot as plt
for c in "ZRT":
    plt.plot(t, tr[c], label=c)
plt.legend(); plt.xlabel("time (s)"); plt.ylabel("velocity (m/s)"); plt.show()

# helpers: moon.double_couple(strike, dip, rake, m0) -> dict mrr..mtp ; moon.point_at(lat, lon, dist, az)

# every model
t, ids, data = moon.synthetics_all_models(stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                                          mrr=1e15, mtt=-5e14, mpp=-5e14)
# data[imodel, icomp, it]
```

Options: `comps="ZNE"` (or any of Z N E R T), `quantity="displacement" | "velocity" | "acceleration"`,
`dt=0.5`, `half_duration=15.0` (Gaussian source, s; keep ≥ 10 s), `force=(fr, ft, fp)` in N
for an impact (e.g. `force=(-1e10, 0, 0)` with `evdp=0`). `moon.max_depth(mid)` gives the
deepest allowed event depth (top of the fluid core) of a model.

Limits: periods ≥ 20 s (0 – 0.05 Hz), elastic models (no attenuation), window 3276.8 s.
Very shallow events (< 5 km) are reliable only beyond ≈ 2° from the station.
