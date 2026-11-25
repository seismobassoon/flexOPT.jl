"""moon_rsgt_py.py -- Moon synthetics from the RSGT catalogue, pure Python (numpy + h5py).

No Julia needed: reads the HDF5 export (rsgt_<id>.h5, made by dev/moonRSGT/export_h5.jl)
and does the same reconstruction as dev/moonRSGT/lib/moonSynthetics.jl (Julia):
reciprocity, explicit 3x3 rotations, Lagrange interpolation (order 3 in depth inside a
layer, order 5 in distance), one-sided inverse Fourier transform.  Only the few
(depth, distance) samples needed are read from disk, so a loop over 1000 models is fast.

    import moon_rsgt_py as moon
    moon.CATALOGUE = "/Volumes/MoonDisk/catalogue_h5"          # folder with rsgt_00001.h5 ...

    t, tr = moon.synthetics("00817", stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                            mrr=1e15, mtt=-5e14, mpp=-5e14, mrt=0, mrp=0, mtp=0)
    tr["Z"], tr["R"], tr["T"]          # m/s (default quantity="velocity")

    for mid in moon.model_ids():      # every model
        t, tr = moon.synthetics(mid, ...)

Conventions
-----------
* station (stla, stlo), event (evla, evlo) in degrees; evdp in km (0 ... top of the fluid core);
* moment tensor in GCMT order and convention (r = up, t = south, p = east), N.m;
  optional force (fr, ft, fp) in N (impacts: e.g. fr = -1e10 at evdp = 0);
* components: Z (up), N, E, R (horizontal, away from the event), T (normal to the great circle);
* quantity: "displacement" (m), "velocity" (m/s, default), "acceleration" (m/s^2);
* source time function: Heaviside (step) moment/force smoothed by a Gaussian of half
  duration `half_duration` (s); keep it >= 10 s (catalogue band 0 - 0.05 Hz).
"""
import os

import h5py
import numpy as np

CATALOGUE = os.environ.get("MOON_RSGT_H5_DIR", "/gpfs/scratch/nobuaki/moonRSGT/catalogue_h5")
FORMAT = "moonRSGT-h5 v1"


# ---------------------------------------------------------------------------- files

def model_ids(catalogue=None):
    """Sorted list of model ids ("00001", ...) present in the catalogue folder."""
    d = catalogue or CATALOGUE
    return sorted(f[5:-3] for f in os.listdir(d) if f.startswith("rsgt_") and f.endswith(".h5"))


def model_file(mid, catalogue=None):
    if isinstance(mid, str) and os.path.isfile(mid):
        return mid
    mid = f"{int(mid):05d}" if str(mid).isdigit() else str(mid)
    return os.path.join(catalogue or CATALOGUE, f"rsgt_{mid}.h5")


def _str(x):
    x = x[()] if hasattr(x, "shape") else x
    return x.decode() if isinstance(x, bytes) else str(x)


def max_depth(mid, catalogue=None):
    """Deepest source depth (km) of a model (top of its fluid core)."""
    with h5py.File(model_file(mid, catalogue), "r") as f:
        return float(f["max_depth_km"][()])


# ---------------------------------------------------------------------------- geometry
# (same formulas as dev/spherical1DelasticKernel/lib/sgtGeometry1D.jl)

def _sph_unit(th, ph):
    return np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])


def _local_basis(th, ph):
    """Columns: r-hat (up), theta-hat (south), phi-hat (east), Cartesian."""
    return np.column_stack([_sph_unit(th, ph),
                            [np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)],
                            [-np.sin(ph), np.cos(ph), 0.0]])


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _ry(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def _pole_rotation(th, ph):
    return _ry(-th) @ _rz(-ph)


def _receiver_direction(comp, src, rcv):
    """Unit vector in the station's local basis (r, theta, phi) for component comp."""
    if comp == "Z":
        return np.array([1.0, 0.0, 0.0])
    if comp == "N":
        return np.array([0.0, -1.0, 0.0])
    if comp == "E":
        return np.array([0.0, 0.0, 1.0])
    E = _local_basis(*rcv)
    s, rr = _sph_unit(*src), _sph_unit(*rcv)
    t = np.cross(s, rr)
    nt = np.linalg.norm(t)
    if nt < 1e-12:
        raise ValueError("event below the station: R/T undefined, use comps='ZNE'")
    t /= nt
    L = np.cross(t, rr)
    if comp == "R":
        return E.T @ L
    if comp == "T":
        return E.T @ t
    raise ValueError("component must be one of Z N E R T")


def _lagrange(x0, x):
    x = np.asarray(x, dtype=float)
    w = np.ones(len(x))
    for i in range(len(x)):
        for j in range(len(x)):
            if j != i:
                w[i] *= (x0 - x[j]) / (x[i] - x[j])
    return w


def epicentral_distance(lat1, lon1, lat2, lon2):
    a = _sph_unit(np.radians(90 - lat1), np.radians(lon1))
    b = _sph_unit(np.radians(90 - lat2), np.radians(lon2))
    return float(np.degrees(np.arccos(np.clip(a @ b, -1.0, 1.0))))


# ---------------------------------------------------------------------------- spectra

def spectra(mid, stla, stlo, evla, evlo, evdp, mrr=0.0, mtt=0.0, mpp=0.0, mrt=0.0, mrp=0.0, mtp=0.0,
            force=None, comps="ZRT", catalogue=None):
    """Complex spectra (m.s, Dirac source) at the catalogue frequencies.

    Returns (omega, omega_i, tlen, {comp: spectrum}, delta_deg).
    """
    with h5py.File(model_file(mid, catalogue), "r") as f:
        if _str(f["format"]) != FORMAT:
            raise ValueError(f"unexpected file format {_str(f['format'])}")
        R = float(f["R_km"][()])
        dmax = float(f["max_depth_km"][()])
        if not (0.0 <= evdp <= dmax + 1e-9):
            raise ValueError(f"evdp {evdp} km outside the catalogue (0 - {dmax:.2f} km)")
        r_km = f["r_km"][()]
        layer = f["layer"][()]
        th_deg = f["theta_deg"][()]
        omega = f["omega"][()]
        omega_i = float(f["omega_i"][()])
        tlen = float(f["tlen"][()])
        ps_src, ps_cmp = f["pairs_strain_src"][()] - 1, f["pairs_strain_comp"][()] - 1
        pd_src, pd_cmp = f["pairs_disp_src"][()] - 1, f["pairs_disp_comp"][()] - 1

        rcv = (np.radians(90.0 - stla), np.radians(stlo))
        src = (np.radians(90.0 - evla), np.radians(evlo))
        rq = R - evdp
        # receiver = pole source of the disk; event = field point
        Rs = _pole_rotation(*rcv)
        Bm = _local_basis(0.0, 0.0).T @ Rs @ _local_basis(*rcv)
        q = Rs @ _sph_unit(*src)
        thp = np.arccos(np.clip(q[2] / np.linalg.norm(q), -1.0, 1.0))
        php = np.arctan2(q[1], q[0])
        T = _local_basis(*src).T @ Rs.T @ _local_basis(thp, php)
        c, s = np.cos(php), np.sin(php)
        Qphi = np.array([[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]])

        # depth: block (physical layer) containing rq, the upper one on a discontinuity
        starts = [0] + [i for i in range(1, len(layer)) if layer[i] != layer[i - 1]]
        blocks = [(a, b) for a, b in zip(starts, starts[1:] + [len(layer)])]
        cand = [(a, b) for (a, b) in blocks if r_km[a] - 1e-9 <= rq <= r_km[b - 1] + 1e-9]
        a, b = cand[-1]
        rv = r_km[a:b]
        n = len(rv)
        qo = min(3, n - 1)
        j = int(np.clip(np.searchsorted(rv, rq, side="right"), 1, n))
        lo = int(np.clip(j - (qo - 1) // 2, 1, n - qo)) - 1
        r0, r1 = a + lo, a + lo + qo + 1
        wr = _lagrange(rq, r_km[r0:r1])
        # distance
        thd = np.degrees(thp)
        nt = len(th_deg)
        j = int(np.searchsorted(th_deg, thd, side="right"))
        if 1 <= j <= nt and abs(th_deg[j - 1] - thd) <= 1e-12 * max(1.0, abs(thd)):
            t0, t1, wt = j - 1, j, np.array([1.0])
        else:
            lo = int(np.clip(j - 2, 1, nt - 5)) - 1
            t0, t1 = lo, lo + 6
            wt = _lagrange(thd, th_deg[t0:t1])
        W = np.outer(wr, wt)
        S = f["strain_re"][:, :, r0:r1, t0:t1] + 1j * f["strain_im"][:, :, r0:r1, t0:t1]
        D = f["disp_re"][:, :, r0:r1, t0:t1] + 1j * f["disp_im"][:, :, r0:r1, t0:t1]

    Sw = np.einsum("kprt,rt->kp", S.astype(np.complex128), W)     # (nw, 10)
    Dw = np.einsum("kprt,rt->kp", D.astype(np.complex128), W)     # (nw, 5)
    # moment tensor (r, theta, phi) = GCMT (r, t, p), code units 1e18 N.m; force 1e15 N
    Mm = np.array([[mrr, mrt, mrp], [mrt, mtt, mtp], [mrp, mtp, mpp]], dtype=float) / 1e18
    F = np.zeros(3) if force is None else np.asarray(force, dtype=float) / 1e15
    sidx = np.array([[0, 3, 4], [3, 1, 5], [4, 5, 2]])            # strain (rr, tt, pp, rt, rp, tp) -> 3x3
    out = {}
    nw = len(omega)
    for comp in comps:
        nvec = _receiver_direction(comp, src, rcv)
        bvec = Qphi @ (Bm @ nvec)                                   # force at the pole, rotated by -phi'
        e6 = np.zeros((nw, 6), dtype=complex)
        for ip in range(len(ps_src)):
            e6[:, ps_cmp[ip]] += bvec[ps_src[ip]] * Sw[:, ip]
        u3 = np.zeros((nw, 3), dtype=complex)
        for ip in range(len(pd_src)):
            u3[:, pd_cmp[ip]] += bvec[pd_src[ip]] * Dw[:, ip]
        Ep = e6[:, sidx]                                            # (nw, 3, 3), pole frame
        Eq = np.einsum("ij,kjl,ml->kim", T, Ep, T)                  # event's local basis
        Uq = u3 @ T.T
        out[comp] = 1e3 * (np.einsum("kij,ij->k", Eq, Mm) + Uq @ F)
    delta = epicentral_distance(evla, evlo, stla, stlo)
    return omega, omega_i, tlen, out, delta


def time_series(omega, omega_i, tlen, spec, dt=1.0, stf="step", half_duration=10.0,
                quantity="velocity", taper=0.2):
    """Real time series of one spectrum (same transform as the Julia moonTimeSeries)."""
    p = {"displacement": 0, "velocity": 1, "acceleration": 2}[quantity]
    wc = omega - 1j * omega_i
    A = np.asarray(spec, dtype=complex).copy()
    if stf == "step":
        A /= 1j * wc
    elif stf != "impulse":
        raise ValueError("stf must be 'step' or 'impulse'")
    if p:
        A *= (1j * wc) ** p
    if half_duration > 0:
        sig = half_duration / 1.628
        A *= np.exp(-(omega * sig) ** 2 / 2)
    nw = len(omega)
    w = np.ones(nw)
    if taper > 0:
        k0 = int(np.floor((1 - taper) * (nw - 1)))
        k = np.arange(k0, nw)
        w[k] = 0.5 * (1 + np.cos(np.pi * (k - k0) / max(1, nw - 1 - k0)))
    t = np.arange(0.0, tlen - dt / 2, dt)
    ph = np.exp(1j * np.outer(t, omega[1:]))                      # (nt, nw-1)
    s = 0.5 * w[0] * A[0].real + (ph * (w[1:] * A[1:])).real.sum(axis=1)
    return t, s * (2 * np.pi / tlen / np.pi) * np.exp(omega_i * t)


def synthetics(mid, stla, stlo, evla, evlo, evdp, mrr=0.0, mtt=0.0, mpp=0.0, mrt=0.0, mrp=0.0, mtp=0.0,
               force=None, comps="ZRT", quantity="velocity", dt=1.0, half_duration=10.0, stf="step",
               catalogue=None):
    """Time-domain synthetics of one model: (t, {comp: trace}).  See the module docstring."""
    om, oi, tlen, sp, _ = spectra(mid, stla, stlo, evla, evlo, evdp, mrr, mtt, mpp, mrt, mrp, mtp,
                                  force=force, comps=comps, catalogue=catalogue)
    tr = {}
    t = None
    for c in comps:
        t, tr[c] = time_series(om, oi, tlen, sp[c], dt=dt, stf=stf, half_duration=half_duration,
                               quantity=quantity)
    return t, tr


def synthetics_all_models(stla, stlo, evla, evlo, evdp, ids=None, catalogue=None, **kw):
    """Same source/station for every model: (t, ids, data[imodel, icomp, it])."""
    ids = ids or model_ids(catalogue)
    comps = kw.get("comps", "ZRT")
    data, t = [], None
    for mid in ids:
        t, tr = synthetics(mid, stla, stlo, evla, evlo, evdp, catalogue=catalogue, **kw)
        data.append([tr[c] for c in comps])
    return t, ids, np.array(data)


def double_couple(strike, dip, rake, m0=1.0):
    """GCMT moment tensor (dict mrr ... mtp, N.m if m0 is) of a double couple; angles in degrees."""
    f, d, l = np.radians([strike, dip, rake])
    return dict(mrr=m0 * np.sin(2 * d) * np.sin(l),
                mtt=-m0 * (np.sin(d) * np.cos(l) * np.sin(2 * f) + np.sin(2 * d) * np.sin(l) * np.sin(f) ** 2),
                mpp=m0 * (np.sin(d) * np.cos(l) * np.sin(2 * f) - np.sin(2 * d) * np.sin(l) * np.cos(f) ** 2),
                mrt=-m0 * (np.cos(d) * np.cos(l) * np.cos(f) + np.cos(2 * d) * np.sin(l) * np.sin(f)),
                mrp=m0 * (np.cos(d) * np.cos(l) * np.sin(f) - np.cos(2 * d) * np.sin(l) * np.cos(f)),
                mtp=-m0 * (np.sin(d) * np.cos(l) * np.cos(2 * f) + 0.5 * np.sin(2 * d) * np.sin(l) * np.sin(2 * f)))


def point_at(lat, lon, delta_deg, az_deg):
    """(lat, lon) at distance delta_deg and azimuth az_deg (clockwise from north) from (lat, lon)."""
    la, lo, d, a = np.radians([lat, lon, delta_deg, az_deg])
    la2 = np.arcsin(np.sin(la) * np.cos(d) + np.cos(la) * np.sin(d) * np.cos(a))
    lo2 = lo + np.arctan2(np.sin(a) * np.sin(d) * np.cos(la), np.cos(d) - np.sin(la) * np.sin(la2))
    return float(np.degrees(la2)), float((np.degrees(lo2) + 180) % 360 - 180)


if __name__ == "__main__":
    import sys
    mid = sys.argv[1] if len(sys.argv) > 1 else model_ids()[0]
    t, tr = synthetics(mid, stla=-3.0, stlo=-23.4, evla=20.0, evlo=15.0, evdp=900.0,
                       mrr=1e15, mtt=-5e14, mpp=-5e14)
    for c, v in tr.items():
        print(c, f"max |v| = {np.abs(v).max():.3e} m/s")
