"""Sparse-bin (ANR style) neutrino maps + features, from the cached nu_*.jld2 runs.
Usage: python3 nu_sparse_features.py <dir with nu_*.jld2> <output dir> [exposure_Mt_yr]"""
import sys, glob, os, re, h5py, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
src, out = sys.argv[1], sys.argv[2]; EXPO = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
os.makedirs(out, exist_ok=True)
RUN_EXPO = 20.0                        # Mt·yr of the cached runs (NU.massKg 20e9, 1 yr)
R, RCMB, RICB, D670 = 6371.0, 3480.0, 1221.5, 670.0
models = ['PREM', 'G1m20_sb11_03140', 'G1m20_sb21_03240']; short = {'PREM': 'PREM', models[1]: 'sb11', models[2]: 'sb21'}
cols = {'PREM': 'black', models[1]: 'firebrick', models[2]: 'dodgerblue'}
dets = ['north', 'east', 'south', 'west']; sides = ['m', 'p']
detcol = {'north': 'gold', 'east': 'limegreen', 'south': 'magenta', 'west': 'cyan'}
D = {}
for f in glob.glob(os.path.join(src, 'nu_*.jld2')):
    m = re.match(r'nu_(.+)_(north|east|south|west)_([mp])_[0-9a-f]+\.jld2', os.path.basename(f))
    h = h5py.File(f, 'r')
    D[(m.group(1), m.group(2), m.group(3))] = {k: h[k][()].T for k in ('Ne', 'Nmu')} | {'E': h['energies'][()], 'c': h['cosθ'][()]}
E, c = D[('PREM', 'north', 'm')]['E'], D[('PREM', 'north', 'm')]['c']
le = np.log10(E); dle = le[1] - le[0]; hc = c[1] - c[0]
fe, fc = round(0.1 / dle), round(0.05 / hc); nE, nC = len(E) // fe, len(c) // fc
Ee = 10 ** np.concatenate([le[:fe * nE:fe] - dle / 2, [le[fe * nE - 1] + dle / 2]]); ce = np.concatenate([c[:fc * nC:fc] - hc / 2, [c[fc * nC - 1] + hc / 2]])
Ec, cc = np.sqrt(Ee[1:] * Ee[:-1]), (ce[1:] + ce[:-1]) / 2
coarse = lambda N: EXPO / RUN_EXPO * N[:fe * nE, :fc * nC].reshape(nE, fe, nC, fc).sum(axis=(1, 3))
C = {k: {ch: coarse(v[ch]) for ch in ('Ne', 'Nmu')} for k, v in D.items()}
def zclass(cz):
    rmin = R * np.sqrt(max(0, 1 - cz ** 2))
    return 'inner core' if rmin < RICB else 'outer core' if rmin < RCMB else 'lower mantle' if rmin < R - D670 else 'upper mantle'
cls = np.array([zclass(x) for x in cc]); classes = [k for k in ['upper mantle', 'lower mantle', 'outer core', 'inner core'] if (cls == k).any()]
bands = [(1, 3), (3, 8), (8, 20), (20, 100)]
chname = {'Ne': 'νe+ν̄e', 'Nmu': 'νμ+ν̄μ'}

def features(N, N0):
    """Features of an event map N (vs PREM N0) per zenith class."""
    F = {}
    for k in classes:
        s = cls == k
        if not s.any(): continue
        n = N[:, s]; n0 = N0[:, s]; tot = n.sum()
        F[('rate', k)] = tot
        F[('dN/N', k)] = tot / n0.sum() - 1
        F[('E centroid', k)] = (np.log10(Ec)[:, None] * n).sum() / tot          # log10 GeV
        F[('cosθ centroid', k)] = (cc[s][None, :] * n).sum() / tot
        r = n / np.maximum(n0, 1e-30) - 1
        F[('rms bin δN/N', k)] = np.sqrt(np.average(r ** 2, weights=n0))
        for b in bands:
            sel = (Ec >= b[0]) & (Ec < b[1])
            F[(f'dN/N {b[0]}–{b[1]} GeV', k)] = n[sel].sum() / max(n0[sel].sum(), 1e-30) - 1
    return F

def chi2(A, B, mask=None):
    m = 0.5 * (A + B); ok = m > 0
    if mask is not None: ok &= mask
    return ((A - B) ** 2 / np.where(ok, m, 1))[ok].sum()

# ---- feature table
rows = []
for ch in ('Ne', 'Nmu'):
    for d in dets:
        for s in sides:
            N0 = C[('PREM', d, s)][ch]
            for m in models:
                for (fn, k), v in features(C[(m, d, s)][ch], N0).items():
                    rows.append((ch, d, s, m, fn, k, v))
import csv
with open(os.path.join(out, 'nu_sparse_features.csv'), 'w') as fo:
    w = csv.writer(fo); w.writerow(['channel', 'detector', 'side', 'model', 'feature', 'class', f'value_{EXPO}Mtyr'])
    w.writerows(rows)
val = {(ch, d, s, m, fn, k): v for ch, d, s, m, fn, k, v in rows}
fnames = ['rate', 'dN/N', 'E centroid', 'cosθ centroid', 'rms bin δN/N'] + [f'dN/N {b[0]}–{b[1]} GeV' for b in bands]

# ---- bar plots (seismic style): per feature, x = zenith class; top: mean value (3 models);
#      bottom: scatter over the 8 detector sides of (model − PREM), × = sb11/sb21; and √Δχ²
for ch in ('Ne', 'Nmu'):
    fig, axs = plt.subplots(len(fnames), 2, figsize=(15, 3.1 * len(fnames)))
    for i, fn in enumerate(fnames):
        x = np.arange(len(classes)); w = 0.27
        for j, m in enumerate(models):
            mu = [np.mean([val[(ch, d, s, m, fn, k)] for d in dets for s in sides]) for k in classes]
            axs[i, 0].bar(x + (j - 1) * w, mu, w, color=cols[m], label=short[m])
        axs[i, 0].set_ylabel(fn); axs[i, 0].set_xticks(x); axs[i, 0].set_xticklabels(classes)
        if fn in ('E centroid', 'cosθ centroid'):
            lo = min(p.get_height() for p in axs[i, 0].patches); hi = max(p.get_height() for p in axs[i, 0].patches)
            axs[i, 0].set_ylim(lo - 0.02 * abs(lo) - 1e-3, hi + 0.02 * abs(hi) + 1e-3)
        w = 0.38; sc = []
        for j, m in enumerate(models[1:]):
            v = [np.std([val[(ch, d, s, m, fn, k)] - val[(ch, d, s, 'PREM', fn, k)] for d in dets for s in sides]) if fn not in ('dN/N', 'rms bin δN/N') and not fn.startswith('dN/N')
                 else np.sqrt(np.mean([val[(ch, d, s, m, fn, k)] ** 2 for d in dets for s in sides])) for k in classes]
            sc.append(v); axs[i, 1].bar(x + (j - 0.5) * w, v, w, color=cols[m], label=short[m])
        for kk in range(len(classes)):
            if sc[1][kk] > 0: axs[i, 1].text(kk, max(sc[0][kk], sc[1][kk]) * 1.03, f'×{sc[0][kk] / sc[1][kk]:.1f}', ha='center', fontsize=8)
        axs[i, 1].set_xticks(x); axs[i, 1].set_xticklabels(classes)
        axs[i, 1].set_ylabel('rms' if fn.startswith('dN/N') or fn == 'rms bin δN/N' else 'std (model−PREM)')
        if i == 0: axs[i, 0].legend(fontsize=8); axs[i, 1].legend(fontsize=8)
    axs[0, 0].set_title(f'{chname[ch]} CC, sparse bins, {EXPO:g} Mt·yr: mean over 4 detectors × 2 sides')
    axs[0, 1].set_title('spread of the deviation from PREM over the 8 detector sides  (×: sb11/sb21)')
    plt.tight_layout(); plt.savefig(os.path.join(out, f'nu_05_feature_bars_{"e" if ch == "Ne" else "mu"}.png'), dpi=110); plt.close()

# ---- √Δχ² sb11 vs sb21 per detector and zenith class (both sides), and vs PREM
fig, axs = plt.subplots(1, 2, figsize=(15, 5))
for a, ch in zip(axs, ('Ne', 'Nmu')):
    x = np.arange(len(classes)); w = 0.2
    for j, d in enumerate(dets):
        v = [np.sqrt(sum(chi2(C[(models[1], d, s)][ch], C[(models[2], d, s)][ch], (cls == k)[None, :].repeat(nE, 0)) for s in sides)) for k in classes]
        a.bar(x + (j - 1.5) * w, v, w, color=detcol[d], edgecolor='k', lw=0.5, label=d)
    a.axhline(1, color='gray', ls=':'); a.axhline(3, color='gray', ls='--')
    a.set_xticks(x); a.set_xticklabels(classes); a.set_ylabel('√Δχ²  sb11 vs sb21 (σ)')
    a.set_title(f'{chname[ch]}: separation per detector, sparse bins, {EXPO:g} Mt·yr'); a.legend()
plt.tight_layout(); plt.savefig(os.path.join(out, 'nu_06_chi2_per_detector.png'), dpi=110); plt.close()

# ---- sparse maps with features overlaid
for ch in ('Ne', 'Nmu'):
    for d in dets:
        fig, axs = plt.subplots(2, 4, figsize=(22, 9.5))
        for r, s in enumerate(sides):
            N0, N1, N2 = (C[(m, d, s)][ch] for m in models)
            sig = (N1 - N2) / np.sqrt(np.maximum(0.5 * (N1 + N2), 1e-30))
            panels = [(N1 - N0, 'sb11 − PREM', 'RdBu_r'), (N2 - N0, 'sb21 − PREM', 'RdBu_r'), (N1 - N2, 'sb11 − sb21', 'RdBu_r'),
                      (sig, '(sb11 − sb21)/√N  [σ per bin]', 'PuOr_r')]
            for k_, (Z, t, cm) in enumerate(panels):
                a = axs[r, k_]; lim = max(np.abs(Z).max(), 1e-12) if k_ < 3 else max(np.abs(Z).max(), 1.0)
                im = a.pcolormesh(Ee, ce, Z.T, cmap=cm, vmin=-lim, vmax=lim, rasterized=True)
                cs = a.contour(Ec, cc, N0.T, levels=np.quantile(N0[N0 > 0], [0.5, 0.8, 0.95]), colors='0.35', linewidths=0.6)
                a.clabel(cs, fmt='%.0f', fontsize=6)
                for rr, nm in ((RCMB, 'core'), (RICB, 'inner core'), (R - D670, '670')):
                    cz = -np.sqrt(1 - (rr / R) ** 2); a.axhline(cz, color='0.3', lw=0.6, ls='--')
                    if k_ == 0: a.text(1.05, cz + 0.008, nm, fontsize=7, color='0.3')
                a.set_xscale('log'); a.set_title(f'{d} side {"−" if s == "m" else "+"}: {t}', fontsize=10)
                if k_ == 0: a.set_ylabel('cos θ')
                if r == 1: a.set_xlabel('E (GeV)')
                plt.colorbar(im, ax=a, shrink=0.85, pad=0.01)
                # features: peak |σ| bin, per-class δN/N (sb11, sb21) and √Δχ²
                if k_ == 3:
                    i, j = np.unravel_index(np.argmax(np.abs(sig)), sig.shape)
                    a.plot(Ec[i], cc[j], marker='*', ms=13, mfc='none', mec='k', mew=1.2)
                    txt = [f'peak {sig[i, j]:+.2f}σ at {Ec[i]:.1f} GeV, cosθ {cc[j]:.2f}']
                    for kk in classes:
                        msk = (cls == kk)[None, :].repeat(nE, 0)
                        txt.append(f'{kk}: √Δχ² {np.sqrt(chi2(N1, N2, msk)):.2f}')
                    a.text(0.02, 0.98, '\n'.join(txt), transform=a.transAxes, va='top', fontsize=7.5,
                           bbox=dict(fc='white', ec='0.6', alpha=0.85))
                if k_ in (0, 1):
                    Nm = N1 if k_ == 0 else N2
                    txt = [f'{kk}: δN/N {100 * (Nm[:, cls == kk].sum() / N0[:, cls == kk].sum() - 1):+.2f} %' for kk in classes]
                    txt.append(f'√Δχ² vs PREM {np.sqrt(chi2(Nm, N0)):.2f}')
                    a.text(0.02, 0.98, '\n'.join(txt), transform=a.transAxes, va='top', fontsize=7.5,
                           bbox=dict(fc='white', ec='0.6', alpha=0.85))
        fig.suptitle(f'{chname[ch]} CC, {d} detector: sparse bins Δlog₁₀E = 0.1, Δcosθ = 0.05, {EXPO:g} Mt·yr perfect detector '
                     f'(grey contours: N_PREM per bin; ★ peak significance)', fontsize=13)
        plt.tight_layout(); plt.savefig(os.path.join(out, f'nu_04_sparse_{d}_{"e" if ch == "Ne" else "mu"}.png'), dpi=100); plt.close()
print('done', nE, nC)
