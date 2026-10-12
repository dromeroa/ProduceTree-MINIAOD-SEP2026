#!/usr/bin/env python3
"""Espectros de las variables globales del jet (señal vs fondo) con estilo de figura de CMS,
sin el rótulo/logo de CMS. Lee UN solo CSV por muestra y solo los primeros N jets.

Variables (se buscan por nombre, sin distinguir mayúsculas, espacios, guiones ni guiones bajos;
las que no existan en el CSV se omiten con un aviso):
    JetPt, JetPt_raw, JetEta, JetPhi, JetMass, JetMass_raw, SDMass, Tau2/Tau1
Tau2/Tau1 (τ21) se calcula jet a jet con las columnas Tau2 y Tau1 (los jets con Tau1 = 0 o no
finito se excluyen de ese histograma y se avisa). Si tus columnas se llaman distinto:
    --map SDMass=SD_mass --map Tau1=Tau_1 --map Tau2=Tau_2
(claves: JetPt JetPt_raw JetEta JetPhi JetMass JetMass_raw SDMass Tau1 Tau2 Tau21).
Con --vars JetPt SDMass Tau21 eliges cuáles graficar.

Estilo: mplhep "CMS" (ticks hacia dentro, sin rejilla ni título; sin él se imita con rcParams),
histogramas en escalón, bandas estadísticas sqrt(N) por bin, ejes en inglés con unidades,
eje y en log para pT y lineal desde 0 para el resto (--logy-all / --linear-y para forzarlo),
panel inferior opcional señal/fondo (--ratio). No se dibuja el rótulo "CMS".

Rango de masa: JetMass, JetMass_raw y SDMass se grafican entre 0 y 300 GeV (--mass-max y --range
para cambiarlo). Con --cut-mass además se descartan, para TODAS las variables, los jets cuya
JetMass cae fuera de [0, 300] GeV (selección a nivel de jet); sin esa opción solo se acota el eje.
JetPt: además de la figura normalizada, se guarda jet_JetPt_counts (--counts-vars cambia cuáles):
el eje y es el número de jets por bin, sin normalizar (barras de error sqrt(N)).

Jets madre: JetPt_raw ya NO se grafica (solo si lo pides con --vars JetPt_raw). Cada leyenda indica
cuántos jets madre se usaron en señal y fondo, y se imprime y guarda el resumen en
--outdir/n_jets_utilizados.txt (leídos, usados, y los constituyentes PF).

Constituyentes PF (además de todo lo anterior): histogramas de N_constituents, PF_Pt, PF_Eta, PF_Phi,
PF_E, PF_PdgId, PF_Dxy, PF_Dz, PF_charge y PF_DeltaR_jet de los primeros N jets, en --outdir como
pf_<variable>.png/.pdf y, con --grid, pf_vars_grid. Se descartan en todas las variables los
constituyentes de relleno (PF_Pt <= 0). N_constituents se toma por jet y se compara con los
constituyentes realmente guardados (aviso si no coinciden: relleno, truncado u orden de filas).
PF_Pt y PF_E van en log-log, PF_Dxy y PF_Dz en y log con rango ±percentil 99.5 de |x| (unidades
en cm: ajusta la etiqueta si tu CSV usa otra), PdgId y carga como barras con la fracción por valor.
Opciones: --no-pf, --only-pf, --pf-vars, --pf-range VAR MIN MAX, --pf-linear, --pf-no-logx
y --map PF_Dz=MiColumna para columnas con otro nombre.

Figura extra (y también en el resumen --grid): JetPt de 0 a 2000 GeV (--pt-zoom MIN MAX), con el
eje y = número de jets por bin (sin normalizar) en escala LINEAL normal, desde 0 (--zoom-bins).

SIN pesos: todos los jets cuentan 1; no se lee GenWeight (el CSV puede no tenerlo). Ojo: si el
MC de fondo se generó plano en p̂T, el espectro sin pesos NO es el físico (el pT del fondo sale
plano hasta el límite cinemático). --norm area (defecto) normaliza cada muestra a 1 y compara
formas; --norm none dibuja el número de jets por unidad de la variable.

Salida: una figura por variable en --outdir (PNG y PDF) y, con --grid, un resumen con todas.

Formato del CSV (se detecta solo): "listas" (una fila por jet, PF_* como listas) o "filas"
(una fila por constituyente, N_constituents repetido; las variables globales se leen de la
primera fila de cada jet). No se leen los constituyentes, solo las variables del jet.

Uso:
    python plot_jet_vars_cms.py --signal senal.csv --background fondo.csv --n 3000
    python plot_jet_vars_cms.py --signal senal.csv --background fondo.csv --ratio --grid
    python plot_jet_vars_cms.py --signal senal.csv --background fondo.csv --vars JetPt SDMass Tau21
"""
import argparse
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullFormatter, ScalarFormatter

try:                                   # estilo oficial de CMS (pip install mplhep)
    import mplhep as hep
except ImportError:
    hep = None

JET_PT_NAMES = ["jetpt", "ak8pt", "fatjetpt", "ptjet", "jetptcorr", "jetptcorrected", "pt"]
JET_ETA_NAMES = ["jeteta", "ak8eta", "fatjeteta", "etajet", "eta"]
JET_TAGS = ("jet", "ak8", "fatjet")
_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|inf", re.I)


def norm(name):
    """'jet PT', 'Jet_pT' y 'jetPt' -> 'jetpt'."""
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def to_array(x):
    """Convierte una celda (texto '[1.2, 3.4]', lista o array) en un array float."""
    if isinstance(x, str):
        return np.array(_FLOAT_RE.findall(x), dtype=float)
    return np.atleast_1d(np.asarray(x, dtype=float))


class Columns:
    """Columnas de un CSV con búsqueda tolerante a mayúsculas, espacios y guiones."""

    def __init__(self, path):
        self.path = path
        self.names = list(pd.read_csv(path, nrows=0).columns)
        self.by_norm = {}
        for c in self.names:
            self.by_norm.setdefault(norm(c), c)

    def need(self, key, pretty):
        if key not in self.by_norm:
            raise SystemExit(f"{self.path}: no encontré la columna {pretty}.\n"
                             f"Columnas: {self.names}")
        return self.by_norm[key]

    def _explicit(self, explicit):
        if explicit in self.names:
            return explicit
        if norm(explicit) in self.by_norm:
            return self.by_norm[norm(explicit)]
        raise SystemExit(f"{self.path}: no existe la columna '{explicit}'.\n"
                         f"Columnas: {self.names}")

    def jet_col(self, explicit, candidates, key, what, flag):
        """Columna global del jet: la indicada, un nombre conocido o una única candidata."""
        if explicit:
            return self._explicit(explicit)
        for name in candidates:
            if name in self.by_norm:
                return self.by_norm[name]
        fuzzy = [c for c in self.names
                 if norm(c).endswith(key) and not norm(c).startswith("pf")
                 and any(t in norm(c) for t in JET_TAGS)]
        if len(fuzzy) == 1:
            return fuzzy[0]
        hint = f" Candidatas: {fuzzy}." if fuzzy else ""
        raise SystemExit(f"{self.path}: no pude identificar la columna del {what} del jet; "
                         f"indícala con {flag}.{hint}\nColumnas: {self.names}")


def detect_format(path, pt_col, eta_col):
    """'listas' si PF_Pt/PF_Eta son texto con varios valores; 'filas' si son numéricas."""
    head = pd.read_csv(path, nrows=5, usecols=[pt_col, eta_col])
    numeric = all(pd.api.types.is_numeric_dtype(head[c]) for c in head.columns)
    return "filas" if numeric else "listas"


def jet_starts(path, n, ncol, chunksize=500_000):
    """Fila donde empieza cada uno de los primeros n jets y total de filas que ocupan
    (formato 'filas'). Salta de jet en jet con N_constituents, sin cargar las demás
    columnas. Si el archivo tiene menos de n jets, devuelve los que haya."""
    parts = []
    starts = []
    for chunk in pd.read_csv(path, usecols=[ncol], chunksize=chunksize):
        parts.append(chunk[ncol].to_numpy())
        nc = np.concatenate(parts)
        starts, i = [], 0
        while len(starts) < n and i < len(nc):
            starts.append(i)
            i += max(int(nc[i]), 1)
        if len(starts) == n and i <= len(nc):
            return np.array(starts), i
    return np.array(starts), len(nc)


VARS = {
    "JetPt": dict(names=["jetpt", "ak8pt", "fatjetpt"], x=r"Jet $p_{T}$ [GeV]",
                  d=r"p_{T}", unit="GeV", kind="pt", logy=True),
    "JetPt_raw": dict(names=["jetptraw", "rawjetpt", "jetrawpt", "jetptuncorr"],
                      x=r"Jet $p_{T}^{\mathrm{raw}}$ [GeV]", d=r"p_{T}^{\mathrm{raw}}",
                      unit="GeV", kind="pt", logy=True),
    "JetEta": dict(names=["jeteta", "ak8eta", "fatjeteta"], x=r"Jet $\eta$",
                   d=r"\eta", unit="", kind="eta", logy=False),
    "JetPhi": dict(names=["jetphi", "ak8phi", "fatjetphi"], x=r"Jet $\phi$ [rad]",
                   d=r"\phi", unit="rad", kind="phi", logy=False),
    "JetMass": dict(names=["jetmass", "ak8mass", "fatjetmass", "mjet"], x=r"Jet mass [GeV]",
                    d=r"m", unit="GeV", kind="mass", logy=False),
    "JetMass_raw": dict(names=["jetmassraw", "rawjetmass", "jetrawmass", "jetmassuncorr"],
                        x=r"Jet $m^{\mathrm{raw}}$ [GeV]", d=r"m^{\mathrm{raw}}", unit="GeV",
                        kind="mass", logy=False),
    "SDMass": dict(names=["sdmass", "softdropmass", "msd", "jetsdmass", "jetmasssd", "sdmasscorr"],
                   x=r"Soft-drop mass [GeV]", d=r"m_{\mathrm{SD}}", unit="GeV", kind="mass",
                   logy=False),
    "Tau21": dict(names=["tau21", "tau2tau1", "jettau21"], x=r"$\tau_{2}/\tau_{1}$",
                  d=r"\tau_{2}/\tau_{1}", unit="", kind="tau", logy=False),
}
TAU_NAMES = {"Tau1": ["tau1", "jettau1"], "Tau2": ["tau2", "jettau2"]}


def find_col(cols, key, names, overrides):
    """Columna del CSV para `key`: la de --map o la primera coincidencia exacta normalizada."""
    if key in overrides:
        return cols._explicit(overrides[key])
    for nm in names:
        if nm in cols.by_norm:
            return cols.by_norm[nm]
    return None


def resolve_vars(cols, wanted, overrides):
    """Devuelve {variable: ('col', nombre) | ('ratio', col2, col1)} solo con las que se pueden hacer."""
    plan = {}
    for key in wanted:
        if key == "Tau21":
            t1 = find_col(cols, "Tau1", TAU_NAMES["Tau1"], overrides)
            t2 = find_col(cols, "Tau2", TAU_NAMES["Tau2"], overrides)
            if t1 and t2 and "Tau21" not in overrides:
                plan[key] = ("ratio", t2, t1)
                continue
        c = find_col(cols, key, VARS[key]["names"], overrides)
        if c:
            plan[key] = ("col", c)
    return plan


def load_jets(path, n, wanted, overrides):
    """Primeros n jets: variables globales pedidas (sin leer constituyentes ni GenWeight)."""
    cols = Columns(path)
    pf_pt = cols.need("pfpt", "PF_Pt")
    pf_eta = cols.need("pfeta", "PF_Eta")
    fmt = detect_format(path, pf_pt, pf_eta)
    plan = resolve_vars(cols, wanted, overrides)
    need = []
    for spec in plan.values():
        need += list(spec[1:])
    need = list(dict.fromkeys(need))
    if not need:
        raise SystemExit(f"{path}: no encontré ninguna de las variables pedidas.\nColumnas: {cols.names}")

    if fmt == "listas":
        df = pd.read_csv(path, nrows=n, usecols=need)
        njets = len(df)
        per_jet = lambda c: df[c].to_numpy(float)
    else:
        ncol = cols.need("nconstituents", "N_constituents")
        starts, nrows = jet_starts(path, n, ncol)
        njets = len(starts)
        df = pd.read_csv(path, nrows=nrows, usecols=need)
        per_jet = lambda c: df[c].to_numpy(float)[starts]

    out = {}
    for key, spec in plan.items():
        if spec[0] == "col":
            out[key] = per_jet(spec[1])
        else:
            num, den = per_jet(spec[1]), per_jet(spec[2])
            with np.errstate(divide="ignore", invalid="ignore"):
                out[key] = np.where(den > 0, num / den, np.nan)
    return {"v": out, "plan": plan, "jets": njets, "fmt": fmt}


def cms_style():
    """Estilo de figura de CMS; sin mplhep, imita lo esencial con rcParams."""
    if hep is not None:
        plt.style.use(hep.style.CMS)
        return
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans",
                                                          "DejaVu Sans"],
        "font.size": 14, "axes.labelsize": 16, "legend.fontsize": 13, "legend.frameon": False,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 8, "ytick.major.size": 8, "xtick.minor.size": 4, "ytick.minor.size": 4,
        "axes.grid": False, "figure.dpi": 100})


def histo(x, bins, norm_area):
    """Densidad por unidad de x y su incertidumbre estadística sqrt(N) en cada bin."""
    h = np.histogram(x, bins=bins)[0].astype(float)
    e = np.sqrt(h)
    scale = 1.0 / np.diff(bins)
    if norm_area:
        total = h.sum()
        scale = scale / total if total != 0 else scale
    return h * scale, e * scale


def draw(ax, bins, h, e, color, label, ls="-", lw=2.0, band=True):
    ax.stairs(h, bins, color=color, lw=lw, ls=ls, label=label)
    if band:
        ax.stairs(h + e, bins, baseline=np.maximum(h - e, 0), fill=True, color=color, alpha=0.25, lw=0)


def var_range(key, arrs, user, mass_max=300.0):
    """Rango común de x para señal y fondo."""
    if key in user:
        return user[key]
    allv = np.concatenate(arrs)
    kind = VARS[key]["kind"]
    if kind == "phi":
        return -np.pi, np.pi
    if kind == "tau":
        return 0.0, 1.0
    if kind == "eta":
        m = np.ceil(np.abs(allv).max() * 2) / 2
        return -m, m
    if kind == "mass":
        return 0.0, float(mass_max)
    return float(allv.min()), float(allv.max())


def plot_var(ax, rax, key, data, args, ranges, logy, colors):
    """Dibuja un histograma (señal y fondo) en `ax`; devuelve True si dibujó algo."""
    spec = VARS[key]
    lo, hi = ranges[key]
    use_log_x = args.logx and spec["kind"] == "pt" and lo > 0
    bins = (np.logspace(np.log10(lo), np.log10(hi), args.bins + 1) if use_log_x
            else np.linspace(lo, hi, args.bins + 1))
    norm_area = args.norm == "area"
    hs = {}
    for label, d in data.items():
        x = d["v"][key]
        xf = x[np.isfinite(x)]
        h, e = histo(xf, bins, norm_area)
        hs[label] = (h, e)
        n_in = int(np.histogram(xf, bins=bins)[0].sum())
        name = f"{label}: {len(xf)} jets" + ("" if n_in == len(xf) else f" ({n_in} in range)")
        draw(ax, bins, h, e, colors[label], name)

    unit = f" [1/{spec['unit']}]" if spec["unit"] else ""
    ax.set_ylabel((rf"$1/N\ \mathrm{{d}}N/\mathrm{{d}}{spec['d']}$" + unit) if norm_area
                  else "Jets" + (f" / {spec['unit']}" if spec["unit"] else ""))
    ax.set_xlim(lo, hi)
    if logy:
        ax.set_yscale("log")
        pos = [v for h, _ in hs.values() for v in h if v > 0]
        if pos:
            ax.set_ylim(min(pos) * 0.3, max(pos) * 10 ** 1.3)
    else:
        top = max(max(h + e) for h, e in hs.values())
        ax.set_ylim(0, top * 1.35)
    if use_log_x:
        ax.set_xscale("log")
        ticks = [t for t in (100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000) if lo <= t <= hi]
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(loc="upper right", title=args.text or None)

    if rax is not None:
        hsig, esig = hs["Signal"]
        hbkg, ebkg = hs["Background"]
        ok = (hbkg > 0) & (hsig > 0)
        sg, bk = np.where(hsig > 0, hsig, 1), np.where(hbkg > 0, hbkg, 1)
        ratio = np.where(ok, hsig / bk, np.nan)
        err = ratio * np.sqrt(np.where(ok, (esig / sg) ** 2 + (ebkg / bk) ** 2, np.nan))
        ctr = 0.5 * (bins[1:] + bins[:-1])
        rax.errorbar(ctr, ratio, yerr=err, fmt="o", color="black", ms=3.5, lw=1)
        rax.axhline(1, color="gray", lw=1, ls="--")
        if spec["kind"] in ("eta", "phi"):       # cocientes ~1: eje lineal
            rax.set_ylim(bottom=0)
        else:                                       # cocientes que varían órdenes de magnitud
            rax.set_yscale("log")
        rax.set_ylabel("Signal / Bkg.")
        rax.set_xlabel(spec["x"])
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(spec["x"])


def plot_counts(key, data, args, rng, colors, logy, nbins, ax=None):
    """Número de jets por bin (sin normalizar), señal y fondo, con error sqrt(N).
    rng = (min, max) del eje x. Si se pasa `ax` dibuja ahí; si no, crea y devuelve una figura."""
    spec = VARS[key]
    lo, hi = rng
    use_log_x = args.logx and spec["kind"] == "pt" and lo > 0
    bins = (np.logspace(np.log10(lo), np.log10(hi), nbins + 1) if use_log_x
            else np.linspace(lo, hi, nbins + 1))
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 9))
    tops = []
    for label, d in data.items():
        x = d["v"][key]
        x = x[np.isfinite(x)]
        h = np.histogram(x, bins=bins)[0].astype(float)
        name = f"{label}: {len(x)} jets" + ("" if int(h.sum()) == len(x) else f" ({int(h.sum())} in range)")
        draw(ax, bins, h, np.sqrt(h), colors[label], name)
        tops.append(h)
        out = np.mean((x < lo) | (x > hi)) * 100 if len(x) else 0
        if out > 1 and fig is not None:
            print(f"AVISO: {label}: {out:.1f}% de {key} queda fuera de [{lo:.4g}, {hi:.4g}] en la figura de conteos")
    width = (hi - lo) / nbins
    unit = f" {spec['unit']}" if spec["unit"] else ""
    ax.set_ylabel("Number of jets" + ("" if use_log_x else f" / {width:.3g}{unit}"))
    ax.set_xlabel(spec["x"])
    ax.set_xlim(lo, hi)
    if logy:
        ax.set_yscale("log")
        pos = [v for h in tops for v in h if v > 0]
        ax.set_ylim(0.5 if pos else 1, max(pos) * 10 ** 1.3 if pos else 10)
    else:                                              # escala lineal normal, y desde 0
        ax.set_ylim(0, max(h.max() for h in tops) * 1.35)
    if use_log_x:
        ax.set_xscale("log")
        ticks = [t for t in (100, 200, 300, 500, 700, 1000, 2000, 3000, 5000, 7000) if lo <= t <= hi]
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(loc="upper right", title=args.text or None)
    return fig


def write_summary(args, lines, mode):
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "n_jets_utilizados.txt"), mode, encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def run_jet_plots(args, overrides, user_range):
    wanted = args.vars or [k for k in VARS if k != "JetPt_raw"]   # JetPt_raw se omite (solo con --vars)
    data = {}
    for label, path in [("Signal", args.signal), ("Background", args.background)]:
        d = load_jets(path, args.n, list(dict.fromkeys(wanted + (["JetMass"] if args.cut_mass else []))),
                      overrides)
        if args.cut_mass:
            if "JetMass" not in d["v"]:
                raise SystemExit(f"{path}: --cut-mass necesita la columna JetMass (usa --map JetMass=COLUMNA).")
            m = d["v"]["JetMass"]
            keep = np.isfinite(m) & (m >= 0) & (m <= args.mass_max)
            print(f"{label}: --cut-mass deja {keep.sum()} de {len(m)} jets con JetMass en [0, {args.mass_max:g}] GeV")
            d["v"] = {k: v[keep] for k, v in d["v"].items()}
        d["n_used"] = len(next(iter(d["v"].values())))
        data[label] = d
        print(f"{label}: formato '{d['fmt']}'; {d['jets']} jets madre leídos (--n {args.n}), "
              f"{d['n_used']} usados")
        for key, spec in d["plan"].items():
            src = spec[1] if spec[0] == "col" else f"{spec[1]} / {spec[2]}"
            bad = np.sum(~np.isfinite(d["v"][key]))
            print(f"        {key} <- {src}" + (f"  ({bad} jets sin valor válido, excluidos)" if bad else ""))
        miss = [k for k in wanted if k not in d["plan"]]
        if miss:
            print(f"        AVISO: no encontré {miss} en este CSV; se omiten (usa --map CLAVE=COLUMNA).")

    keys = [k for k in wanted if all(k in d["plan"] for d in data.values())]
    if not keys:
        raise SystemExit("Ninguna variable está disponible en ambos CSV.")
    ranges = {}
    for key in keys:
        ranges[key] = var_range(key, [np.concatenate([d["v"][key][np.isfinite(d["v"][key])]
                                                      for d in data.values()])], user_range,
                              args.mass_max)
        lo, hi = ranges[key]
        for label, d in data.items():
            x = d["v"][key][np.isfinite(d["v"][key])]
            out = np.mean((x < lo) | (x > hi)) * 100 if len(x) else 0
            if out > 1:
                print(f"AVISO: {label}: {out:.1f}% de {key} queda fuera de [{lo:.4g}, {hi:.4g}]")

    os.makedirs(args.outdir, exist_ok=True)
    summary = [f"Jets madre utilizados (primeros --n {args.n} de cada CSV)"]
    for label, d in data.items():
        summary.append(f"  {label}: {d['jets']} leídos, {d['n_used']} usados"
                       + (f" (tras --cut-mass 0-{args.mass_max:g} GeV)" if args.cut_mass else ""))
    print("\n" + "\n".join(summary))
    write_summary(args, summary, "w")
    cms_style()
    colors = {"Signal": "#0072B2", "Background": "#D55E00"}   # validados: CVD y contraste >= 3:1
    logy_of = lambda k: False if args.linear_y else (True if args.logy_all else VARS[k]["logy"])

    for key in keys:
        if args.ratio:
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 11), sharex=True,
                                          gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        else:
            fig, ax = plt.subplots(figsize=(10, 9))
            rax = None
        plot_var(ax, rax, key, data, args, ranges, logy_of(key), colors)
        base = os.path.join(args.outdir, f"jet_{key}")
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"Figura guardada: {base}.png (+ .pdf)")

    for key in args.counts_vars:
        if key not in keys:
            print(f"AVISO: {key} no está disponible; no se hace la figura de número de jets.")
            continue
        fig = plot_counts(key, data, args, ranges[key], colors, logy_of(key), args.bins)
        base = os.path.join(args.outdir, f"jet_{key}_counts")
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"Figura guardada: {base}.png (+ .pdf)  [número de jets, sin normalizar]")

    zoom = "JetPt" in keys
    if zoom:
        zr = (args.pt_zoom[0], args.pt_zoom[1])
        fig = plot_counts("JetPt", data, args, zr, colors, False, args.zoom_bins)
        base = os.path.join(args.outdir, f"jet_JetPt_counts_{zr[0]:g}-{zr[1]:g}")
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"Figura guardada: {base}.png (+ .pdf)  [número de jets, eje y lineal]")
    else:
        print("AVISO: JetPt no está disponible; no se hace la figura de JetPt 0-2000 GeV.")

    if args.grid:
        npanel = len(keys) + (1 if zoom else 0)
        ncol = (4 if npanel <= 8 else 3) if npanel > 3 else npanel
        nrow = int(np.ceil(npanel / ncol))
        plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "legend.fontsize": 11,
                             "xtick.labelsize": 12, "ytick.labelsize": 12})
        fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 4.8 * nrow), squeeze=False)
        for ax, key in zip(axes.flat, keys):
            plot_var(ax, None, key, data, args, ranges, logy_of(key), colors)
        if zoom:
            plot_counts("JetPt", data, args, zr, colors, False, args.zoom_bins, ax=axes.flat[len(keys)])
        for ax in list(axes.flat)[npanel:]:
            ax.axis("off")
        fig.tight_layout()
        base = os.path.join(args.outdir, "jet_vars_grid")
        fig.savefig(base + ".png", dpi=170, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        print(f"Figura resumen: {base}.png (+ .pdf)")


PFVARS = {
    "N_constituents": dict(names=["nconstituents"], x="Number of constituents per jet",
                           d=r"N_{\mathrm{const}}", unit="", kind="count", logy=False),
    "PF_Pt": dict(names=["pfpt"], x=r"PF constituent $p_{T}$ [GeV]", d=r"p_{T}", unit="GeV",
                  kind="logpos", logy=True),
    "PF_Eta": dict(names=["pfeta"], x=r"PF constituent $\eta$", d=r"\eta", unit="",
                   kind="eta", logy=False),
    "PF_Phi": dict(names=["pfphi"], x=r"PF constituent $\phi$ [rad]", d=r"\phi", unit="rad",
                   kind="phi", logy=False),
    "PF_E": dict(names=["pfe", "pfenergy"], x=r"PF constituent $E$ [GeV]", d=r"E", unit="GeV",
                 kind="logpos", logy=True),
    "PF_PdgId": dict(names=["pfpdgid", "pfpdgidabs"], x="PF PDG ID", d="", unit="",
                     kind="discrete", logy=False),
    "PF_Dxy": dict(names=["pfdxy"], x=r"PF $d_{xy}$ [cm]", d=r"d_{xy}", unit="cm",
                   kind="tail", logy=True),
    "PF_Dz": dict(names=["pfdz"], x=r"PF $d_{z}$ [cm]", d=r"d_{z}", unit="cm",
                  kind="tail", logy=True),
    "PF_charge": dict(names=["pfcharge"], x="PF charge", d="", unit="", kind="discrete",
                      logy=False),
    "PF_DeltaR_jet": dict(names=["pfdeltarjet", "pfdrjet", "pfdeltar"],
                          x=r"PF $\Delta R$ to the jet axis", d=r"\Delta R", unit="",
                          kind="dr", logy=False),
}


def load_pf(path, n, wanted, overrides):
    """Constituyentes PF de los primeros n jets. Los constituyentes con PF_Pt <= 0 (relleno) se
    descartan en TODAS las variables. N_constituents se toma por jet. También comprueba que
    N_constituents coincida con el número de constituyentes reales de cada jet."""
    cols = Columns(path)
    pf_pt = cols.need("pfpt", "PF_Pt")
    pf_eta = cols.need("pfeta", "PF_Eta")
    fmt = detect_format(path, pf_pt, pf_eta)
    found = {}
    for key in wanted:
        c = find_col(cols, key, PFVARS[key]["names"], overrides)
        if c:
            found[key] = c
    found["PF_Pt"] = pf_pt
    percol = [k for k in found if k != "N_constituents"]
    need = list(dict.fromkeys(list(found.values())))

    if fmt == "listas":
        df = pd.read_csv(path, nrows=n, usecols=need)
        njets = len(df)
        cells = {k: [to_array(x) for x in df[found[k]]] for k in percol}
        lengths = np.array([len(a) for a in cells["PF_Pt"]], dtype=int)
        for k in percol:
            lk = np.array([len(a) for a in cells[k]], dtype=int)
            if (lk != lengths).any():
                raise SystemExit(f"{path}: {found[k]} y {pf_pt} tienen distinta longitud en "
                                 f"{(lk != lengths).sum()} jets (posible error al llenar las ramas).")
        arrs = {k: (np.concatenate(cells[k]) if njets else np.array([])) for k in percol}
        jet_id = np.repeat(np.arange(njets), lengths)
        nconst = df[found["N_constituents"]].to_numpy(float) if "N_constituents" in found else None
        list_len = lengths
    else:
        ncol = cols.need("nconstituents", "N_constituents")
        found.setdefault("N_constituents", ncol)
        need = list(dict.fromkeys(need + [ncol]))
        starts, nrows = jet_starts(path, n, ncol)
        njets = len(starts)
        df = pd.read_csv(path, nrows=nrows, usecols=need)
        arrs = {k: df[found[k]].to_numpy(float) for k in percol}
        blocks = np.diff(np.append(starts, nrows))
        jet_id = np.repeat(np.arange(njets), blocks)
        nconst = df[ncol].to_numpy(float)[starts]
        list_len = blocks

    valid = (arrs["PF_Pt"] > 0) & np.isfinite(arrs["PF_Pt"])
    padded = 1 - valid.mean() if len(valid) else 0.0
    arrs = {k: v[valid] for k, v in arrs.items()}
    jet_id = jet_id[valid]
    n_valid = np.bincount(jet_id, minlength=njets)
    mism = float(np.mean(nconst != n_valid) * 100) if nconst is not None and njets else None
    mism_len = float(np.mean(nconst != list_len) * 100) if nconst is not None and njets else None
    return {"pf": arrs, "nconst": nconst, "found": found, "jets": njets, "fmt": fmt,
            "padded": padded * 100, "mism_valid": mism, "mism_len": mism_len,
            "n_valid": n_valid}


def pf_values(d, key):
    x = d["nconst"] if key == "N_constituents" else d["pf"][key]
    return x[np.isfinite(x)]


def pf_range(key, data, user):
    if key in user:
        return user[key]
    allv = np.concatenate([pf_values(d, key) for d in data.values()])
    kind = PFVARS[key]["kind"]
    if kind == "count":
        return 0.0, float(allv.max())
    if kind == "logpos":
        pos = allv[allv > 0]
        return float(pos.min()), float(pos.max())
    if kind == "eta":
        m = np.ceil(np.abs(allv).max() * 2) / 2
        return -m, m
    if kind == "phi":
        return -np.pi, np.pi
    if kind == "tail":
        m = float(np.percentile(np.abs(allv), 99.5))
        m = m if m > 0 else 1.0
        return -m, m
    if kind == "dr":
        return 0.0, float(np.ceil(allv.max() * 10) / 10)
    return float(allv.min()), float(allv.max())


def pf_bins(key, rng, args):
    lo, hi = rng
    kind = PFVARS[key]["kind"]
    if kind == "count" and hi - lo <= 150:
        return np.arange(np.floor(lo) - 0.5, np.ceil(hi) + 1.5, 1.0)
    if kind == "logpos" and not args.pf_no_logx and lo > 0:
        return np.logspace(np.log10(lo), np.log10(hi), args.bins + 1)
    return np.linspace(lo, hi, args.bins + 1)


def plot_discrete(ax, key, data, args, colors):
    """Variables discretas (PDG ID, carga): fracción de constituyentes por valor."""
    spec = PFVARS[key]
    vals = {lab: np.round(pf_values(d, key)).astype(int) for lab, d in data.items()}
    fr = {}
    for lab, v in vals.items():
        u, c = np.unique(v, return_counts=True)
        fr[lab] = dict(zip(u, c / max(len(v), 1)))
        fr[lab]["_n"] = len(v)
    cats = {c for lab in fr for c in fr[lab] if c != "_n" and fr[lab][c] >= 0.001}
    cats = sorted(cats, key=(lambda c: c) if key == "PF_charge" else (lambda c: (abs(c), c)))[:14]
    pos = np.arange(len(cats))
    w = 0.38
    for i, (lab, d) in enumerate(data.items()):
        h = np.array([fr[lab].get(c, 0.0) for c in cats])
        e = np.sqrt(h * fr[lab]["_n"]) / max(fr[lab]["_n"], 1)
        off = (i - 0.5) * w
        ax.bar(pos + off, h, width=w, color=colors[lab], alpha=0.9,
               label=f"{lab}: {d['jets']} jets, {fr[lab]['_n']} PF")
        ax.errorbar(pos + off, h, yerr=e, fmt="none", ecolor="black", lw=1, capsize=0)
    ax.set_xticks(pos)
    ax.set_xticklabels([str(c) for c in cats], rotation=0 if len(cats) <= 8 else 45)
    ax.set_xlabel(spec["x"])
    ax.set_ylabel("Fraction of constituents")
    ax.set_ylim(0, max(max(fr[l].get(c, 0.0) for c in cats) for l in fr) * 1.35)
    ax.legend(loc="upper right", title=args.text or None)
    shown = {lab: sum(fr[lab].get(c, 0.0) for c in cats) * 100 for lab in fr}
    return cats, shown


def plot_pf_var(ax, rax, key, data, args, ranges, colors):
    """Histograma de una variable PF (señal y fondo) en `ax`; ratio opcional en `rax`."""
    spec = PFVARS[key]
    if spec["kind"] == "discrete":
        plot_discrete(ax, key, data, args, colors)
        return
    lo, hi = ranges[key]
    bins = pf_bins(key, (lo, hi), args)
    logx = spec["kind"] == "logpos" and not args.pf_no_logx and lo > 0
    norm_area = args.norm == "area"
    hs = {}
    for label, d in data.items():
        x = pf_values(d, key)
        h, e = histo(x, bins, norm_area)
        hs[label] = (h, e)
        n_in = int(np.histogram(x, bins=bins)[0].sum())
        name = (f"{label}: {n_in} jets" if key == "N_constituents"
                else f"{label}: {d['jets']} jets, {n_in} PF")
        draw(ax, bins, h, e, colors[label], name)
    unit = f" [1/{spec['unit']}]" if spec["unit"] else ""
    if norm_area:
        ax.set_ylabel(rf"$1/N\ \mathrm{{d}}N/\mathrm{{d}}{spec['d']}$" + unit)
    else:
        ax.set_ylabel("Number of jets" if key == "N_constituents" else "Number of PF constituents")
    ax.set_xlim(lo, hi)
    logy = spec["logy"] and not args.pf_linear
    if logy:
        ax.set_yscale("log")
        pos = [v for h, _ in hs.values() for v in h if v > 0]
        if pos:
            ax.set_ylim(min(pos) * 0.3, max(pos) * 10 ** 1.3)
    else:
        ax.set_ylim(0, max(max(h + e) for h, e in hs.values()) * 1.35)
    if logx:
        ax.set_xscale("log")
    ax.legend(loc="upper right", title=args.text or None)
    if rax is not None:
        hsig, esig = hs["Signal"]
        hbkg, ebkg = hs["Background"]
        ok = (hbkg > 0) & (hsig > 0)
        sg, bk = np.where(hsig > 0, hsig, 1), np.where(hbkg > 0, hbkg, 1)
        ratio = np.where(ok, hsig / bk, np.nan)
        err = ratio * np.sqrt(np.where(ok, (esig / sg) ** 2 + (ebkg / bk) ** 2, np.nan))
        ctr = np.sqrt(bins[1:] * bins[:-1]) if logx else 0.5 * (bins[1:] + bins[:-1])
        rax.errorbar(ctr, ratio, yerr=err, fmt="o", color="black", ms=3.5, lw=1)
        rax.axhline(1, color="gray", lw=1, ls="--")
        if spec["kind"] in ("eta", "phi"):
            rax.set_ylim(bottom=0)
        else:
            rax.set_yscale("log")
        rax.set_ylabel("Signal / Bkg.")
        rax.set_xlabel(spec["x"])
        ax.tick_params(labelbottom=False)
    else:
        ax.set_xlabel(spec["x"])


def run_pf_plots(args, overrides):
    """Figuras de los constituyentes PF (una por variable, en --outdir) y resumen con --grid."""
    pf_user = {}
    for key, lo, hi in args.pf_range:
        if key not in PFVARS:
            raise SystemExit(f"--pf-range: variable desconocida '{key}'. Opciones: {list(PFVARS)}")
        pf_user[key] = (float(lo), float(hi))
    wanted = args.pf_vars or list(PFVARS)
    print("\n== Constituyentes PF ==")
    data = {}
    for label, path in [("Signal", args.signal), ("Background", args.background)]:
        d = load_pf(path, args.n, wanted, overrides)
        data[label] = d
        print(f"{label}: formato '{d['fmt']}', {d['jets']} jets, {len(d['pf']['PF_Pt'])} constituyentes "
              f"(se descartó {d['padded']:.1f}% de relleno con PF_Pt <= 0)")
        for key in wanted:
            print(f"        {key} <- {d['found'].get(key, 'NO ENCONTRADA')}")
        if d["mism_valid"] is not None:
            print(f"        N_constituents != nº de constituyentes con PF_Pt>0 en {d['mism_valid']:.1f}% de los jets; "
                  f"!= largo de la lista/bloque en {d['mism_len']:.1f}%")
            if d["mism_valid"] > 1:
                print("        AVISO: N_constituents no coincide con los constituyentes guardados; "
                      "revisa el llenado del tree/CSV (relleno, truncado a un máximo, orden de filas).")
    pf_sum = ["Constituyentes PF utilizados"]
    for label, d in data.items():
        pf_sum.append(f"  {label}: {d['jets']} jets madre, {len(d['pf']['PF_Pt'])} constituyentes")
    print("\n".join(pf_sum))
    write_summary(args, pf_sum, "w" if args.only_pf else "a")
    keys = [k for k in wanted if all(k in d["found"] and (k == "N_constituents" or k in d["pf"])
                                     for d in data.values())]
    miss = [k for k in wanted if k not in keys]
    if miss:
        print(f"AVISO: variables PF no disponibles en ambos CSV, se omiten: {miss} (usa --map CLAVE=COLUMNA).")
    if not keys:
        raise SystemExit("Ninguna variable PF está disponible en ambos CSV.")

    ranges = {}
    for key in keys:
        if PFVARS[key]["kind"] == "discrete":
            continue
        ranges[key] = pf_range(key, data, pf_user)
        lo, hi = ranges[key]
        for label, d in data.items():
            x = pf_values(d, key)
            out = np.mean((x < lo) | (x > hi)) * 100 if len(x) else 0
            if out > 1:
                print(f"AVISO: {label}: {out:.1f}% de {key} queda fuera de [{lo:.4g}, {hi:.4g}]")

    os.makedirs(args.outdir, exist_ok=True)
    cms_style()
    colors = {"Signal": "#0072B2", "Background": "#D55E00"}
    for key in keys:
        discrete = PFVARS[key]["kind"] == "discrete"
        if args.ratio and not discrete:
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(10, 11), sharex=True,
                                          gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
        else:
            fig, ax = plt.subplots(figsize=(10, 9))
            rax = None
        plot_pf_var(ax, rax, key, data, args, ranges, colors)
        base = os.path.join(args.outdir, f"pf_{key}")
        fig.savefig(base + ".png", dpi=200, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        plt.close(fig)
        print(f"Figura guardada: {base}.png (+ .pdf)")

    if args.grid:
        npanel = len(keys)
        ncol = npanel if npanel <= 3 else (4 if npanel <= 8 else (3 if npanel == 9 else 5))
        nrow = int(np.ceil(npanel / ncol))
        plt.rcParams.update({"font.size": 13, "axes.labelsize": 15, "legend.fontsize": 11,
                             "xtick.labelsize": 12, "ytick.labelsize": 12})
        fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 4.8 * nrow), squeeze=False)
        for ax, key in zip(axes.flat, keys):
            plot_pf_var(ax, None, key, data, args, ranges, colors)
        for ax in list(axes.flat)[npanel:]:
            ax.axis("off")
        fig.tight_layout()
        base = os.path.join(args.outdir, "pf_vars_grid")
        fig.savefig(base + ".png", dpi=170, bbox_inches="tight")
        fig.savefig(base + ".pdf", bbox_inches="tight")
        print(f"Figura resumen PF: {base}.png (+ .pdf)")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal", required=True, help="CSV de la señal")
    ap.add_argument("--background", required=True, help="CSV del fondo")
    ap.add_argument("--n", type=int, default=3000, help="jets por muestra")
    ap.add_argument("--vars", nargs="+", choices=list(VARS), help="variables del jet a graficar (defecto: todas)")
    ap.add_argument("--map", action="append", default=[], metavar="CLAVE=COLUMNA",
                    help="nombre real de una columna, p. ej. --map SDMass=SD_mass (repetible)")
    ap.add_argument("--norm", choices=["area", "none"], default="area",
                    help="area: cada muestra a 1 (defecto); none: número de jets por unidad")
    ap.add_argument("--bins", type=int, default=50)
    ap.add_argument("--range", action="append", nargs=3, default=[], metavar=("VAR", "MIN", "MAX"),
                    help="rango de una variable, p. ej. --range SDMass 0 300 (repetible)")
    ap.add_argument("--logx", action="store_true", help="eje x log en JetPt y JetPt_raw")
    ap.add_argument("--logy-all", action="store_true", help="eje y log en todas las variables")
    ap.add_argument("--linear-y", action="store_true", help="eje y lineal en todas (también pT)")
    ap.add_argument("--mass-max", type=float, default=300.0,
                    help="límite superior (GeV) del eje de las masas (defecto 300; mínimo 0)")
    ap.add_argument("--cut-mass", action="store_true",
                    help="descarta los jets con JetMass fuera de [0, mass-max] en TODAS las variables")
    ap.add_argument("--counts-vars", nargs="*", default=["JetPt"], choices=list(VARS),
                    help="variables que además se grafican como número de jets sin normalizar "
                         "(defecto: JetPt; vacío para ninguna)")
    ap.add_argument("--pt-zoom", nargs=2, type=float, default=[0.0, 2000.0], metavar=("MIN", "MAX"),
                    help="rango (GeV) de la figura extra de JetPt en número de jets, eje y lineal "
                         "(defecto 0 2000)")
    ap.add_argument("--zoom-bins", type=int, default=50, help="bins de esa figura extra (defecto 50)")
    ap.add_argument("--no-pf", action="store_true", help="no graficar los constituyentes PF")
    ap.add_argument("--only-pf", action="store_true", help="graficar solo los constituyentes PF")
    ap.add_argument("--pf-vars", nargs="+", choices=list(PFVARS), help="variables PF a graficar (defecto: todas)")
    ap.add_argument("--pf-range", action="append", nargs=3, default=[], metavar=("VAR", "MIN", "MAX"),
                    help="rango de una variable PF, p. ej. --pf-range PF_Dxy -1 1 (repetible)")
    ap.add_argument("--pf-linear", action="store_true", help="eje y lineal en TODAS las variables PF")
    ap.add_argument("--pf-no-logx", action="store_true", help="eje x lineal también en PF_Pt y PF_E")
    ap.add_argument("--ratio", action="store_true", help="panel inferior con señal/fondo")
    ap.add_argument("--text", default="", help="título de la leyenda (p. ej. 'AK8 jets')")
    ap.add_argument("--grid", action="store_true", help="además, una figura resumen con todas")
    ap.add_argument("--outdir", default="cms_jet_vars")
    args = ap.parse_args()

    overrides = {}
    for m in args.map:
        if "=" not in m:
            raise SystemExit(f"--map espera CLAVE=COLUMNA, recibí '{m}'")
        k, v = m.split("=", 1)
        overrides[k.strip()] = v.strip()
    user_range = {}
    for key, lo, hi in args.range:
        if key not in VARS:
            raise SystemExit(f"--range: variable desconocida '{key}'. Opciones: {list(VARS)}")
        user_range[key] = (float(lo), float(hi))

    if not args.only_pf:
        run_jet_plots(args, overrides, user_range)
    if not args.no_pf:
        run_pf_plots(args, overrides)


if __name__ == "__main__":
    main()
