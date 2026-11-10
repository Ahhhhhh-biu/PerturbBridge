"""Shared helpers for the figures in 15_manuscript/figures: fonts/style, and complex classes of the dictionary RBPs derived from Complex Portal
(06_features/rbp_functional/complexportal/9606.tsv, UniProt IDs from rbp_id_map.tsv). Classes are assigned in the priority order of CLASS_DEFS."""
import re, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import font_manager
PB = "/well/blagden/users/eic789/PerturbBridge"; FIG = f"{PB}/15_manuscript/figures"
for f in ("Arial.ttf", "Arial-Bold.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "savefig.dpi": 300, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
# (class label, Complex Portal accessions, name pattern or None); first match wins
CLASS_DEFS = [("U2 snRNP (SF3A/SF3B)", ["CPX-2539"], None), ("U4/U6.U5 tri-snRNP", ["CPX-2391", "CPX-26421"], None), ("Exon junction complex", ["CPX-1941"], None),
              ("Catalytic spliceosome (B-act/C)", ["CPX-26451", "CPX-26463"], None), ("Other spliceosomal", [], r"spliceosom|small nuclear ribonucleoprotein|U2AF"),
              ("Ribosomal proteins", ["CPX-5223", "CPX-5183"], None), ("eIF3", ["CPX-6036"], None), ("SSU processome", ["CPX-2511"], None)]
CLASS_COL = {"U2 snRNP (SF3A/SF3B)": "#c0392b", "U4/U6.U5 tri-snRNP": "#e67e22", "Exon junction complex": "#8e44ad", "Catalytic spliceosome (B-act/C)": "#d4ac0d", "Other spliceosomal": "#f1948a",
             "Ribosomal proteins": "#2471a3", "eIF3": "#17a589", "SSU processome": "#7f8c8d", "other": "#d9d9d9"}
REPRO_COL = {"high": "#c0392b", "mid": "#e59866", "low": "#bfc5cc", "na": "#bfc5cc"}
def complexes():
    ids = pd.read_csv(f"{PB}/06_features/rbp_functional/rbp_id_map.tsv", sep="\t"); u2t = dict(zip(ids.uniprot, ids.target)); cp = pd.read_csv(f"{PB}/06_features/rbp_functional/complexportal/9606.tsv", sep="\t"); out = {}
    for ac, name, parts in zip(cp["#Complex ac"], cp["Recommended name"], cp["Expanded participant list"].fillna("")):
        m = {u2t[u] for u in {re.split(r"[-(]", p)[0] for p in parts.split("|")} if u in u2t}
        if m: out[ac] = (name, m)
    return out
def rbp_classes():
    cx = complexes(); cls = {}
    for label, acs, pat in CLASS_DEFS:
        mem = set().union(*[cx[a][1] for a in acs if a in cx]) if acs else set()
        if pat: mem |= set().union(*[m for n, m in cx.values() if re.search(pat, n, re.I)])
        for t in mem: cls.setdefault(t, label)
    return cls
def shared_complex():
    cx = complexes(); mem = {}
    for ac, (n, m) in cx.items():
        for t in m: mem.setdefault(t, set()).add(ac)
    return lambda a, b: len(mem.get(a, set()) & mem.get(b, set())) > 0
def save(fig, name): fig.savefig(f"{FIG}/{name}.png"); fig.savefig(f"{FIG}/{name}.pdf"); print("saved", f"{FIG}/{name}.png")
def class_legend(ax, present, **kw):
    from matplotlib.patches import Patch
    return ax.legend(handles=[Patch(color=CLASS_COL[c], label=c) for c in CLASS_COL if c in present], **kw)
def repel_labels(ax, xs, ys, texts, fontsize=5.2, dx=4, pad=1.15, **kw):
    """Greedy label placement: labels go to the right of their point (left if they would leave the axes) and are pushed up until they do not
    overlap a previously placed label; a thin connector is drawn when a label was moved."""
    fig = ax.figure; fig.canvas.draw(); r = fig.canvas.get_renderer(); tr = ax.transData; inv = tr.inverted(); h = fontsize * fig.dpi / 72 * pad; placed = []; axb = ax.get_window_extent(r)
    pts = sorted(zip(xs, ys, texts), key=lambda p: tr.transform((p[0], p[1]))[1])
    for x, y, t in pts:
        px, py = tr.transform((x, y)); w = len(t) * fontsize * fig.dpi / 72 * 0.58; left = px + dx * fig.dpi / 72 + w > axb.x1; x0 = px - dx * fig.dpi / 72 - w if left else px + dx * fig.dpi / 72; y0 = py - h / 2; moved = False
        for _ in range(60):
            hit = [b for b in placed if not (x0 + w < b[0] or b[0] + b[2] < x0 or y0 + h < b[1] or b[1] + h < y0)]
            if not hit: break
            y0 = max(b[1] for b in hit) + h; moved = True
        placed.append((x0, y0, w)); tx, ty = inv.transform((x0 + (w if left else 0), y0 + h / 2))
        ax.annotate(t, (x, y), xytext=(tx, ty), textcoords="data", fontsize=fontsize, ha="right" if left else "left", va="center", arrowprops=dict(arrowstyle="-", lw=0.3, color="0.4", shrinkA=0, shrinkB=1.5) if moved else None, zorder=5, **kw)
