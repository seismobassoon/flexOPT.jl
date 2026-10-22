import numpy as np, os, sys
D = os.path.expanduser("~/mnt/MantleConvectionTakashi/op_old_full_mars_2025")
def _valid(m): return (1<=m<=12) or (101<=m<=112) or (301<=m<=312) or (401<=m<=412)
def read_scalar(path):
    buf = open(path, "rb").read(); off = 0
    def rd(dt, n=1):
        nonlocal off
        a = np.frombuffer(buf, dtype=dt, count=n, offset=off); off += a.nbytes
        return a if n > 1 else a[0]
    magic = int(rd("<i4")); ilen = 4
    if not _valid(magic) and magic > 8000:
        magic -= 8000; rd("<i4"); ilen = 8
    magic %= 100
    I = "<i4" if ilen == 4 else "<i8"; F = "<f4" if ilen == 4 else "<f8"
    nx, ny, nz = (int(v) for v in rd(I, 3)); nb = int(rd(I)) if magic >= 7 else 1
    rd(F, 2); nnx, nny, nnz = (int(v) for v in rd(I, 3)); nnb = int(rd(I)) if magic >= 8 else 1
    zg = rd(F, 2*nz+1); rcmb = float(rd(F)); rd(I); time = float(rd(F))
    if magic >= 5: rd(F)
    if magic >= 12: rd(F); rd(F, 2); rd(F, 3)
    if magic >= 6: rd(F)
    if magic >= 10: rd(F)
    if magic >= 11: rd(F)
    rd(F, nx); rd(F, ny); rd(F, nz)
    out = np.zeros((nx, ny, nz, nb))
    nxp, nyp, nzp, nbp = nx//nnx, ny//nny, nz//nnz, nb//nnb
    for ib in range(nnb):
        for iz in range(nnz):
            for iy in range(nny):
                for ix in range(nnx):
                    blk = rd(F, nxp*nyp*nzp*nbp).reshape((nxp, nyp, nzp, nbp), order="F")
                    out[ix*nxp:(ix+1)*nxp, iy*nyp:(iy+1)*nyp, iz*nzp:(iz+1)*nzp, ib*nbp:(ib+1)*nbp] = blk
    return out[0, :, :, 0], rcmb + zg[1::2], rcmb + zg[0::2], rcmb, time
snaps = range(201)
T = np.zeros((201, 1024, 128), np.float16); C = np.zeros_like(T); W = np.zeros_like(T); times = np.zeros(201)
for i in snaps:
    t, r, re, rcmb, time = read_scalar(f"{D}/test_t{i:05d}")
    c = read_scalar(f"{D}/test_bs{i:05d}")[0]; w = read_scalar(f"{D}/test_wtr{i:05d}")[0]
    T[i], C[i], W[i], times[i] = t, c, w, time
np.savez_compressed("stagyy_tbw_201.npz", T=T, C=C, W=W, times=times, r=r, r_edge=re, rcmb=rcmb)
print("times (Myr):", times[[0, 1, 100, 199, 200]] / 3.15576e13)
print("T", float(T.min()), float(T.max()), "C", float(C.min()), float(C.max()), "W", float(W.min()), float(W.max()))
