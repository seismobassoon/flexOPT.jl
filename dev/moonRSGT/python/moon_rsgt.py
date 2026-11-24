"""moon_rsgt.py -- Moon synthetics from the flexOPT RSGT catalogue, from Python.

Thin wrapper around dev/moonRSGT/synth.jl (Julia does the work; Python only
writes the requests and reads the result).  Needs numpy and h5py.

    import sys; sys.path.append("/gpfs/users/nobuaki/github/flexOPT.jl/dev/moonRSGT/python")
    from moon_rsgt import moon_synthetics

    t, data, models = moon_synthetics(
        models="00817",                      # "all", "00001,00817", "1:50"
        stla=-3.0, stlo=-23.4,               # station (deg)
        evla=20.0, evlo=15.0, evdp=900.0,    # event (deg, km)
        mrr=1e15, mtt=-5e14, mpp=-5e14, mrt=0, mrp=0, mtp=0)   # GCMT order, N.m
    # data[ireq, imodel, icomp, it]  (m), comps "ZRT" by default

stla ... mtp may be arrays of equal length (one synthetic per entry); every
synthetic is computed for every model.  Julia must run on a compute node
(S-CAPAD/DANTE rule): inside a Slurm job or the Jupyterhub nothing is added,
otherwise the call is wrapped in `srun` (see `launcher`).
"""
import os
import subprocess
import tempfile

import h5py
import numpy as np

FLEXOPT = os.environ.get("FLEXOPT_ROOT", "/gpfs/users/nobuaki/github/flexOPT.jl")
CATALOGUE = os.environ.get("MOON_RSGT_DIR", "/gpfs/scratch/nobuaki/moonRSGT/catalogue")
SRUN = ["srun", "-p", "ncpushort", "-n1", "-c1", "--mem=4G", "-t", "1:00:00"]
COLS = ["stla", "stlo", "evla", "evlo", "evdp", "mrr", "mtt", "mpp", "mrt", "mrp", "mtp"]


def _str(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def moon_synthetics(models, stla, stlo, evla, evlo, evdp, mrr, mtt, mpp, mrt, mrp, mtp,
                    comps="ZRT", dt=1.0, half_duration=10.0, quantity="displacement",
                    launcher=None, workdir=None, verbose=True):
    """Return (t, data, model_ids); data[ireq, imodel, icomp, it] in m (or m/s, m/s^2).

    comps: any of "ZNERT"; quantity: displacement | velocity | acceleration;
    half_duration: Gaussian source time function (s), keep >= 10 s (catalogue fmax 0.05 Hz).
    launcher: list prefixed to the julia command (default: [] inside a Slurm job, else SRUN).
    """
    cols = np.broadcast_arrays(*[np.atleast_1d(np.asarray(v, dtype=float))
                                 for v in (stla, stlo, evla, evlo, evdp, mrr, mtt, mpp, mrt, mrp, mtp)])
    req = np.column_stack(cols)
    if launcher is None:
        launcher = [] if "SLURM_JOB_ID" in os.environ else SRUN
    workdir = workdir or tempfile.mkdtemp(prefix="moon_rsgt_")
    rfile = os.path.join(workdir, "requests.csv")
    ofile = os.path.join(workdir, "synth.jld2")
    np.savetxt(rfile, req, delimiter=",", header=",".join(COLS), comments="")
    cmd = launcher + ["julia", f"--project={FLEXOPT}", os.path.join(FLEXOPT, "dev", "moonRSGT", "synth.jl"),
                      f"--models={models}", f"--requests={rfile}", f"--out={ofile}", f"--dir={CATALOGUE}",
                      f"--comps={comps}", f"--dt={dt}", f"--hd={half_duration}", f"--quantity={quantity}"]
    if verbose:
        print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    with h5py.File(ofile, "r") as f:
        t = np.asarray(f["t"][()])
        data = np.asarray(f["data"][()])            # (nreq, nmodel, ncomp, nt)
        model_ids = _str(f["models"][()]).split(",")
    return t, data, model_ids
