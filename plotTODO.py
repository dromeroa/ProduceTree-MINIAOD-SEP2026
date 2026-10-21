# Codigo del 20 de Octubre 17:47

#!/usr/bin/env python3
"""
PLOTmulti_v1.py  —  Fondo (QCD) + 5 señales en UN solo gráfico por variable.

Uso:
    python PLOTmulti_v1.py                                  # usa la lista SAMPLES de abajo
    python PLOTmulti_v1.py "QCD=qcd.csv" "Z'=zp.csv" "H=h.csv" ...   # etiqueta=archivo (el 1.º es el fondo)

Qué hace respecto a PLOTconPESOSv5.py:
  - Variables previas (JetPt, JetEta, JetPhi, JetMass, SDMass, Tau21, PF_Pt, PF_Eta, PF_Phi)
    + nuevas: PF_E, PF_PdgId, PF_Dxy, PF_Dz y cualquier columna PF_Dxy_* (p. ej. PF_Dxy_Err, PF_Dxy_Sig...).
  - Cada muestra se normaliza a área 1 (densidad) y se pondera con su GenWeight "limpio"
    (|w|, exclusión de los K jets más pesados SOLO en el fondo, clip a w_max según N_eff objetivo).
  - Banda = error estadístico (PF agregado por jet, como en v4/v5).

Optimizaciones de velocidad:
  [R1] Las 6 muestras se procesan EN PARALELO (un proceso por archivo).
  [R2] Solo se leen las columnas necesarias (usecols); no se exporta CSV ni se hacen tablas de diagnóstico.
  [R3] Parseo de listas PF vectorizado (np.fromstring sobre un único string) en vez de regex por jet,
       y solo para los jets con peso > 0.
  [R4] Binning uniforme por aritmética (floor) + bincount, en vez de np.histogram/searchsorted.
  [R5] Rangos de bins de las variables nuevas: se estiman con una muestra pequeña de cada archivo.
  [R6] Un solo histograma ponderado por variable (sin curva sin pesos / sin clip / resampleo).
"""
import os
import re
import sys
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

# ================= CONFIGURACIÓN =================
# (etiqueta, archivo). El PRIMERO es el fondo; los otros 5 son las señales. EDITA las rutas.
SAMPLES = [
    ("QCD (fondo)", "data_QCD1910.csv"),
    ("Señal 1", "signal1000.csv"),
    ("Señal 2", "signal2000.csv"),
    ("Señal 3", "signal3000.csv"),
    ("Señal 4", "signal4000.csv"),
    ("Señal 5", "signal5000.csv"),
]
OUTDIR = "graficos_multi_v1"
CHUNKSIZE = 100000                  # filas por bloque
MAX_JETS = None                     # entero (p. ej. 50000) para pruebas rápidas; None = todos
N_WORKERS = None                    # None = min(nº de muestras, nº de CPUs)
MIN_NEFF_OBJETIVO = 5000.0          # N_eff (Kish) objetivo para el clip de pesos (por muestra)
EXCLUIR_TOP_K_FONDO = 5             # jets de mayor peso que se excluyen en el FONDO (las señales: 0)
NEFF_MIN_BIN = 100.0                # umbral de aviso de N_eff por bin
NBINS_AUTO = 60                     # nº de bins de las variables con rango automático
N_MUESTRA_AUTO = 20000              # filas por archivo para estimar el rango automático

COLORES = ["black", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#7F7F7F"]
ESTILOS = ["-", "-", "--", "-.", ":", (0, (5, 1, 1, 1)), "-", "--"]

_PI = np.pi
VARS_CONFIG = {
    "JetPt":   {"cands": ["jetpt", "ak8pt", "fatjetpt"], "bins": np.linspace(15, 1000, 60), "xlabel": r"Jet $p_T$ [GeV]", "is_pf": False, "cat": False, "logy": False},
    "JetEta":  {"cands": ["jeteta", "ak8eta", "fatjeteta"], "bins": np.linspace(-3.0, 3.0, 60), "xlabel": r"Jet $\eta$", "is_pf": False, "cat": False, "logy": False},
    "JetPhi":  {"cands": ["jetphi", "ak8phi", "fatjetphi"], "bins": np.linspace(-_PI, _PI, 60), "xlabel": r"Jet $\phi$ [rad]", "is_pf": False, "cat": False, "logy": False},
    "JetMass": {"cands": ["jetmass", "ak8mass", "fatjetmass", "mjet"], "bins": np.linspace(0, 300, 60), "xlabel": "Jet Mass [GeV]", "is_pf": False, "cat": False, "logy": False},
    "SDMass":  {"cands": ["sdmass", "softdropmass", "msd", "jetsdmass"], "bins": np.linspace(0, 300, 60), "xlabel": "Soft-Drop Mass [GeV]", "is_pf": False, "cat": False, "logy": False},
    "Tau21":   {"cands": ["tau21", "tau2tau1", "jettau21"], "bins": np.linspace(0, 1, 50), "xlabel": r"$\tau_2 / \tau_1$", "is_pf": False, "cat": False, "logy": False},
    "PF_Pt":   {"cands": ["pfpt"], "bins": np.linspace(0, 100, 60), "xlabel": r"PF Constituent $p_T$ [GeV]", "is_pf": True, "cat": False, "logy": False},
    "PF_Eta":  {"cands": ["pfeta"], "bins": np.linspace(-3.0, 3.0, 60), "xlabel": r"PF Constituent $\eta$", "is_pf": True, "cat": False, "logy": False},
    "PF_Phi":  {"cands": ["pfphi"], "bins": np.linspace(-_PI, _PI, 60), "xlabel": r"PF Constituent $\phi$ [rad]", "is_pf": True, "cat": False, "logy": False},
    # ---- nuevas ----
    "PF_E":      {"cands": ["pfe", "pfenergy"], "bins": np.linspace(0, 200, 60), "xlabel": r"PF Constituent $E$ [GeV]", "is_pf": True, "cat": False, "logy": True},
    "PF_PdgId":  {"cands": ["pfpdgid", "pfpdg", "pdgid"], "bins": np.arange(-1000.5, 1001.5, 1.0), "xlabel": "PF PDG ID", "is_pf": True, "cat": True, "logy": True},
    "PF_Dxy":    {"cands": ["pfdxy"], "bins": "auto", "xlabel": r"PF $d_{xy}$", "is_pf": True, "cat": False, "logy": True},
    "PF_Dz":     {"cands": ["pfdz"], "bins": "auto", "xlabel": r"PF $d_{z}$", "is_pf": True, "cat": False, "logy": True},
}
# Las columnas PF_Dxy_* (PF_Dxy_Err, PF_Dxy_Sig, ...) se detectan solas y se añaden con rango automático.

_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf", re.I)
_SENT = "-9.87654321e+300"
_SENTV = float(_SENT)
_TR = str.maketrans({c: " " for c in "[](),;"})


def norm(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


# ---------------------------------------------------------------- utilidades
def parse_listas(items):
    """[R3] items: lista de strings tipo '[1.2, 3.4]' (uno por jet). Devuelve (valores planos, índice de jet)."""
    items = [x if isinstance(x, str) else "" for x in items]
    n = len(items)
    flat = None
    try:
        big = (" " + _SENT + " ").join(items).translate(_TR)
        flat = np.fromstring(big, sep=" ", dtype=np.float64)
        pos = np.flatnonzero(flat == _SENTV)
        if len(pos) != n - 1:
            flat = None
    except Exception:
        flat = None
    if flat is None:                                   # respaldo lento pero robusto
        vals, jid = [], []
        for j, s in enumerate(items):
            a = np.array(_FLOAT_RE.findall(s), dtype=np.float64)
            vals.append(a)
            jid.append(np.full(len(a), j, dtype=np.int64))
        if not vals:
            return np.array([]), np.array([], dtype=np.int64)
        return np.concatenate(vals), np.concatenate(jid)
    edges = np.concatenate([[-1], pos, [len(flat)]])
    lens = np.diff(edges) - 1
    vals = flat[flat != _SENTV]
    return vals, np.repeat(np.arange(n, dtype=np.int64), lens)


def bin_idx(v, lo, hi, nb):
    """[R4] Binning uniforme. v debe ser finito. Devuelve (índice, máscara de válidos)."""
    b = np.floor((v - lo) * (nb / (hi - lo))).astype(np.int64)
    b[v == hi] = nb - 1
    return b, (b >= 0) & (b < nb)


def sumsq_por_jet(bb, jj, wj, nb):
    """Σ_jets (n_jet,bin · w_jet)² por bin (varios constituyentes de un jet comparten peso). jj ascendente."""
    out = np.zeros(nb)
    if len(bb) == 0:
        return out
    JB = max(1, 4_000_000 // nb)
    for start in range(int(jj[0]), int(jj[-1]) + 1, JB):
        lo = np.searchsorted(jj, start)
        hi = np.searchsorted(jj, start + JB)
        if hi == lo:
            continue
        key = (jj[lo:hi] - start) * nb + bb[lo:hi]
        s = np.bincount(key, weights=wj[lo:hi], minlength=JB * nb).reshape(JB, nb)
        out += np.square(s).sum(axis=0)
    return out


def encontrar_wmax_para_neff(w_pos, target_neff):
    target = min(float(target_neff), len(w_pos) * 0.8)
    if (w_pos.sum() ** 2) / np.square(w_pos).sum() >= target:
        return float(w_pos.max())
    lo, hi = np.log10(max(w_pos.min(), 1e-300)), np.log10(w_pos.max())
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        wc = np.minimum(w_pos, 10.0 ** mid)
        if (wc.sum() ** 2) / np.square(wc).sum() < target:
            hi = mid
        else:
            lo = mid
    return float(10.0 ** (0.5 * (lo + hi)))


def inspeccionar(path):
    header = pd.read_csv(path, nrows=5)
    cols = list(header.columns)
    by = {}
    for c in cols:
        by.setdefault(norm(c), []).append(c)

    def first(n):
        return by[n][0] if n in by else None

    weight_col = next((first(n) for n in ("genweight", "weight", "eventweight") if n in by), None)
    if weight_col is None:
        raise RuntimeError(f"{path}: no se encontró GenWeight. Columnas: {cols}")
    nconst_col = first("nconstituents")
    pf_pt = first("pfpt")
    fmt = "filas" if (pf_pt and pd.api.types.is_numeric_dtype(header[pf_pt])) else "listas"
    if fmt == "filas" and not nconst_col:
        raise RuntimeError(f"{path}: formato 'filas' requiere la columna N_constituents.")

    dxy_cols = [c for c in cols if norm(c).startswith("pfdxy")]
    base_dxy = (next((c for c in dxy_cols if c.lower() == "pf_dxy"), None)
                or next((c for c in dxy_cols if norm(c) == "pfdxy"), None))
    active = {}
    for key, cfg in VARS_CONFIG.items():
        meta = {k: v for k, v in cfg.items() if k != "cands"}
        col = base_dxy if key == "PF_Dxy" else next((first(c) for c in cfg["cands"] if c in by), None)
        if col:
            active[key] = {"mode": "direct", "cols": [col], **meta}
        elif key == "Tau21" and "tau1" in by and "tau2" in by:
            active[key] = {"mode": "ratio", "cols": [first("tau1"), first("tau2")], **meta}
    for c in dxy_cols:                                    # PF_Dxy_*, PF_Dxy_Err, ...
        if c == base_dxy:
            continue
        suf = re.sub(r"^pf[\W_]*dxy[\W_]*", "", c, flags=re.I)
        active["PF_Dxy_" + suf] = {"mode": "direct", "cols": [c], "bins": "auto", "xlabel": f"PF {c}",
                                   "is_pf": True, "cat": False, "logy": True}
    return {"fmt": fmt, "weight_col": weight_col, "nconst_col": nconst_col, "active": active}


def muestrear(path):
    """[R5] Valores de una muestra pequeña para fijar el rango de las variables con bins='auto'."""
    info = inspeccionar(path)
    auto = {k: c for k, c in info["active"].items() if isinstance(c["bins"], str)}
    out = {}
    if auto:
        cols = sorted({c["cols"][0] for c in auto.values()})
        df = pd.read_csv(path, usecols=cols, nrows=N_MUESTRA_AUTO)
        for k, c in auto.items():
            col = c["cols"][0]
            if info["fmt"] == "listas":
                v, _ = parse_listas(df[col].tolist())
            else:
                v = df[col].to_numpy(dtype=np.float64)
            out[k] = v[np.isfinite(v)]
    return info, out


def escanear(path, fmt, wcol, ncol):
    """Paso 1: inicio, longitud y GenWeight de cada jet (lee solo 1-2 columnas)."""
    usecols = [wcol] + ([ncol] if fmt == "filas" else [])
    starts_l, lens_l, w_l = [], [], []
    curr, next_jet_row, n_jets, completo = 0, 0, 0, True
    for ch in pd.read_csv(path, usecols=usecols, chunksize=CHUNKSIZE):
        n_ch = len(ch)
        w_arr = ch[wcol].to_numpy(dtype=np.float64)
        if fmt == "filas":
            nc_np = np.maximum(ch[ncol].to_numpy(dtype=np.int64), 1)
            nc = nc_np.tolist()
            idx = next_jet_row - curr
            s = []
            while idx < n_ch:
                s.append(idx)
                idx += nc[idx]
                n_jets += 1
                if MAX_JETS and n_jets >= MAX_JETS:
                    break
            s = np.array(s, dtype=np.int64)
            starts_l.append(curr + s)
            lens_l.append(nc_np[s])
            w_l.append(w_arr[s])
            next_jet_row = curr + idx
        else:
            take = n_ch if not MAX_JETS else min(n_ch, MAX_JETS - n_jets)
            starts_l.append(curr + np.arange(take, dtype=np.int64))
            lens_l.append(np.ones(take, dtype=np.int64))
            w_l.append(w_arr[:take])
            n_jets += take
        curr += n_ch
        if MAX_JETS and n_jets >= MAX_JETS:
            completo = False
            break
    starts = np.concatenate(starts_l)
    lens = np.concatenate(lens_l)
    if completo:
        lens = np.minimum(lens, curr - starts)
    return starts, lens, np.concatenate(w_l)


def valores_jet(buf, c, ls):
    if c["mode"] == "direct":
        return buf[c["cols"][0]].to_numpy(dtype=np.float64)[ls]
    t1 = buf[c["cols"][0]].to_numpy(dtype=np.float64)[ls]
    t2 = buf[c["cols"][1]].to_numpy(dtype=np.float64)[ls]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(t1 > 0, t2 / t1, np.nan)


# ---------------------------------------------------------------- worker (una muestra)
def procesar_muestra(args):
    label, path, k_excl, bins_map = args
    t0 = time.time()
    info = inspeccionar(path)
    fmt = info["fmt"]
    active = {k: c for k, c in info["active"].items() if k in bins_map}

    # ---- Paso 1: pesos
    jet_starts, jet_lengths, w_signed = escanear(path, fmt, info["weight_col"], info["nconst_col"])
    jet_ends = jet_starts + jet_lengths
    n_jets = len(w_signed)
    all_w = np.abs(w_signed)
    all_w = np.where(np.isfinite(all_w) & (all_w > 0), all_w, 0.0)
    n_excl = 0
    if k_excl > 0 and n_jets > k_excl:
        all_w[np.argpartition(all_w, -k_excl)[-k_excl:]] = 0.0
        n_excl = int(k_excl)
    pos = all_w > 0
    n_pos = int(pos.sum())
    if n_pos == 0:
        raise RuntimeError(f"{label}: todos los jets tienen GenWeight == 0.")
    w_pos = all_w[pos]
    neff_raw = (w_pos.sum() ** 2) / np.square(w_pos).sum()
    w_max = encontrar_wmax_para_neff(w_pos, MIN_NEFF_OBJETIVO)
    w_cl = np.minimum(all_w, w_max)
    neff_clean = (w_cl.sum() ** 2) / np.square(w_cl).sum()
    n_topados = int((all_w > w_max).sum())
    del all_w, w_signed

    # ---- Paso 2: histogramas
    H = {k: np.zeros(len(bins_map[k]) - 1) for k in active}
    H2 = {k: np.zeros(len(bins_map[k]) - 1) for k in active}
    jet_keys = [k for k, c in active.items() if not c["is_pf"]]
    pf_keys = [k for k, c in active.items() if c["is_pf"]]
    lohinb = {k: (float(bins_map[k][0]), float(bins_map[k][-1]), len(bins_map[k]) - 1) for k in active}

    def bloque(buf, buf_start, jl, jh):
        ls = jet_starts[jl:jh] - buf_start
        ln = jet_lengths[jl:jh]
        w = w_cl[jl:jh]
        ok = w > 0
        nj = jh - jl
        for key in jet_keys:
            lo, hi, nb = lohinb[key]
            vals = valores_jet(buf, active[key], ls)
            m = ok & np.isfinite(vals)
            if not m.any():
                continue
            b, inb = bin_idx(vals[m], lo, hi, nb)
            bb, ww = b[inb], w[m][inb]
            H[key] += np.bincount(bb, weights=ww, minlength=nb)
            H2[key] += np.bincount(bb, weights=ww * ww, minlength=nb)
        if not pf_keys:
            return
        if fmt == "filas":
            offs = np.cumsum(ln) - ln
            tot = int(ln.sum())
            row_idx = np.arange(tot) - np.repeat(offs, ln) + np.repeat(ls, ln)
            jidx = np.repeat(np.arange(nj), ln)
            w_use = w
            pf = {k: (buf[active[k]["cols"][0]].to_numpy(dtype=np.float64)[row_idx], jidx) for k in pf_keys}
        else:
            sel = np.flatnonzero(ok)                       # solo se parsean los jets con peso > 0
            w_use = w[sel]
            ls_sel = ls[sel]
            pf = {k: parse_listas(buf[active[k]["cols"][0]].to_numpy()[ls_sel].tolist()) for k in pf_keys}
        pt = pf["PF_Pt"][0] if "PF_Pt" in pf else None
        for key in pf_keys:
            lo, hi, nb = lohinb[key]
            v, j = pf[key]
            m = np.isfinite(v)
            if pt is not None and len(pt) == len(v):
                m &= pt > 0
            wj = w_use[j]
            m &= wj > 0
            if not m.any():
                continue
            b, inb = bin_idx(v[m], lo, hi, nb)
            bb, ww, jj = b[inb], wj[m][inb], j[m][inb]
            H[key] += np.bincount(bb, weights=ww, minlength=nb)
            if active[key]["cat"]:
                H2[key] += np.bincount(bb, weights=ww * ww, minlength=nb)   # categórica: error por constituyente
            else:
                H2[key] += sumsq_por_jet(bb, jj, ww, nb)

    cols = sorted({col for c in active.values() for col in c["cols"]})
    buf, buf_start, n_proc = None, 0, 0
    for chunk in pd.read_csv(path, usecols=cols, chunksize=CHUNKSIZE):
        buf = chunk if buf is None else pd.concat([buf, chunk], ignore_index=True)
        buf_end = buf_start + len(buf)
        lo_j = int(np.searchsorted(jet_starts, buf_start, side="left"))
        hi_j = int(np.searchsorted(jet_ends, buf_end, side="right"))
        if hi_j > lo_j:
            bloque(buf, buf_start, lo_j, hi_j)
            n_proc += hi_j - lo_j
            new_start = min(int(jet_starts[hi_j]) if hi_j < n_jets else buf_end, buf_end)
            buf = buf.iloc[new_start - buf_start:]
            buf_start = new_start
        if hi_j >= n_jets:
            break

    return {"label": label, "hists": {k: (H[k], H2[k]) for k in active},
            "n_jets": n_jets, "n_pos": n_pos, "n_excl": n_excl, "w_max": w_max, "n_topados": n_topados,
            "neff_raw": float(neff_raw), "neff_clean": float(neff_clean), "n_proc": n_proc,
            "fmt": fmt, "secs": time.time() - t0}


# ---------------------------------------------------------------- main
def main():
    t_ini = time.time()
    samples = SAMPLES
    if len(sys.argv) > 1:
        samples = []
        for a in sys.argv[1:]:
            lab, sep, p = a.partition("=")
            samples.append((lab, p) if sep else (os.path.splitext(os.path.basename(a))[0], a))
    faltan = [p for _, p in samples if not os.path.exists(p)]
    if faltan:
        raise SystemExit("ERROR: no existen estos archivos (edita SAMPLES o pásalos por línea de comandos):\n  "
                         + "\n  ".join(faltan))
    os.makedirs(OUTDIR, exist_ok=True)
    nw = N_WORKERS or min(len(samples), os.cpu_count() or 1)
    print(f"{len(samples)} muestras, {nw} procesos en paralelo.", flush=True)

    with ProcessPoolExecutor(max_workers=nw) as ex:
        # 1) inspección + muestra pequeña para rangos automáticos
        pre = list(ex.map(muestrear, [p for _, p in samples]))
        bins_map = {k: c["bins"] for k, c in VARS_CONFIG.items() if not isinstance(c["bins"], str)}
        meta, auto_vals = {}, {}
        for (lab, _), (info, vals) in zip(samples, pre):
            for k, c in info["active"].items():
                meta.setdefault(k, c)
            for k, v in vals.items():
                auto_vals.setdefault(k, []).append(v)
            faltantes = [k for k in VARS_CONFIG if k not in info["active"]]
            print(f"  [{lab}] formato '{info['fmt']}'" + (f" | variables ausentes: {faltantes}" if faltantes else ""), flush=True)
        for k, lst in auto_vals.items():
            v = np.concatenate(lst)
            if len(v) == 0:
                continue
            lo_p, hi_p = np.percentile(v, 0.5), np.percentile(v, 99.5)
            if v.min() >= 0:
                lo, hi = 0.0, float(hi_p)
            else:
                m = float(max(abs(lo_p), abs(hi_p)))
                lo, hi = -m, m
            if not hi > lo:
                lo, hi = lo - 0.5, hi + 0.5
            bins_map[k] = np.linspace(lo, hi, NBINS_AUTO + 1)
        orden = [k for k in VARS_CONFIG if k in bins_map] + [k for k in meta if k not in VARS_CONFIG and k in bins_map]

        # 2) procesamiento paralelo de las muestras
        tareas = [(lab, p, EXCLUIR_TOP_K_FONDO if i == 0 else 0, bins_map) for i, (lab, p) in enumerate(samples)]
        resultados = list(ex.map(procesar_muestra, tareas))

    print("\nRESUMEN POR MUESTRA:", flush=True)
    print(f"  {'muestra':<16}{'jets':>10}{'peso>0':>10}{'excl.':>6}{'N_eff crudo':>13}{'N_eff limpio':>14}{'topados':>9}{'t [s]':>8}", flush=True)
    for r in resultados:
        print(f"  {r['label']:<16}{r['n_jets']:>10}{r['n_pos']:>10}{r['n_excl']:>6}{r['neff_raw']:>13.0f}"
              f"{r['neff_clean']:>14.0f}{r['n_topados']:>9}{r['secs']:>8.1f}", flush=True)

    # ---- gráficos
    print("\nGenerando gráficos...", flush=True)
    resumen = {}
    for key in orden:
        c = meta[key]
        bins = bins_map[key]
        widths = np.diff(bins)
        fig, ax = plt.subplots(figsize=(9, 6.2), constrained_layout=True)
        top, min_pos = 0.0, np.inf
        curvas = []
        for i, r in enumerate(resultados):
            if key not in r["hists"]:
                continue
            h, h2 = r["hists"][key]
            tot = h.sum()
            if tot <= 0:
                continue
            nb_eff = np.where(h2 > 0, h * h / np.where(h2 > 0, h2, 1.0), np.nan)
            resumen.setdefault(key, {})[r["label"]] = int(((h > 0) & (nb_eff < NEFF_MIN_BIN)).sum())
            curvas.append((i, r, h / tot, np.sqrt(h2) / tot))
        if not curvas:
            plt.close(fig)
            continue
        ns = len(curvas)
        etq = lambda r: (f"{r['label']} ($N_{{eff}}$={r['neff_clean']:.0f}" +
                         (f", {r['n_excl']} jets excl.)" if r["n_excl"] else ")"))
        if c["cat"]:
            M = np.array([cv[2] for cv in curvas])
            occ = np.flatnonzero((M >= 1e-4).any(axis=0))
            codes = np.round(0.5 * (bins[:-1] + bins[1:])).astype(int)[occ]
            x = np.arange(len(occ))
            wd = 0.85 / ns
            for n, (i, r, f, e) in enumerate(curvas):
                ax.bar(x + (n - (ns - 1) / 2) * wd, f[occ], wd, color=COLORES[i % len(COLORES)],
                       edgecolor="black" if i == 0 else "none", alpha=0.85, label=etq(r))
            ax.set_xticks(x)
            ax.set_xticklabels([str(cd) for cd in codes], rotation=60, ha="right")
            ax.set_ylabel("Fracción de constituyentes (ponderada)", fontsize=13)
            ax.set_yscale("log")
            ax.set_ylim(5e-5, max(float(M.max()) * 4, 1e-3))
        else:
            for i, r, f, e in curvas:
                dens, err = f / widths, e / widths
                col, ls_ = COLORES[i % len(COLORES)], ESTILOS[i % len(ESTILOS)]
                if i == 0:
                    ax.stairs(dens, bins, fill=True, color="0.80", alpha=0.7)
                ax.stairs(dens, bins, color=col, linestyle=ls_, linewidth=2.0 if i else 2.2, label=etq(r))
                ax.stairs(dens + err, bins, baseline=np.maximum(dens - err, 0), fill=True, color=col, alpha=0.18)
                top = max(top, float(np.max(dens + err)))
                if (dens > 0).any():
                    min_pos = min(min_pos, float(dens[dens > 0].min()))
            ax.set_xlim(bins[0], bins[-1])
            ax.set_ylabel("Densidad de probabilidad", fontsize=13)
            if c["logy"]:
                ax.set_yscale("log")
                ax.set_ylim(max(min_pos * 0.3, top * 1e-7), top * 4)
            else:
                ax.set_ylim(0, top * 1.35 if top > 0 else 1.0)
        ax.set_xlabel(c["xlabel"], fontsize=13)
        ax.set_title(f"{key} — Fondo vs señales (GenWeight limpio, área = 1)", fontsize=13)
        ax.legend(loc="upper right", fontsize=9, frameon=False)
        ax.tick_params(direction="in", top=True, right=True, length=5)
        safe = re.sub(r"[^A-Za-z0-9_]", "_", key)
        fig.savefig(os.path.join(OUTDIR, f"multi_{safe}.png"), dpi=170)
        plt.close(fig)
        print(f"  -> {OUTDIR}/multi_{safe}.png", flush=True)

    labs = [r["label"] for r in resultados]
    print(f"\nBINS CON N_eff < {NEFF_MIN_BIN:.0f} (errores poco fiables):", flush=True)
    print("  " + f"{'variable':<14}" + "".join(f"{l[:12]:>14}" for l in labs), flush=True)
    for key in orden:
        if key in resumen:
            print("  " + f"{key:<14}" + "".join(f"{resumen[key].get(l, '-'):>14}" for l in labs), flush=True)
    print(f"\n¡Listo en {time.time() - t_ini:.1f} s!", flush=True)


if __name__ == "__main__":
    main()
