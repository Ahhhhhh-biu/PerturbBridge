"""Workflow schematic of the PerturbBridge project (five stages, top to bottom).
All coordinates are in millimetres on a 180 mm wide canvas; the y axis points down.
After drawing, every label is checked against its box (prints FIT WARNING if a label leaves its box).
Usage: envs/pb_gpu_venv/bin/python 15_manuscript/fig_workflow.py -> 15_manuscript/fig_workflow.png (300 DPI), fig_workflow.pdf"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyBboxPatch

for f in ("Arial.ttf", "Arial-Bold.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 6.2, "savefig.dpi": 300, "pdf.fonttype": 42})
O = "/well/blagden/users/eic789/PerturbBridge/15_manuscript"

# palette: (edge, tint) per stage
C = {"data": ("#5B7083", "#EEF2F6"), "art": ("#C98A1B", "#FCF3E0"), "fwd": ("#2F6DB5", "#E8F0FA"),
     "prog": ("#3B8F5A", "#EAF5EE"), "inv": ("#7A4FA3", "#F1EAF7"), "out": ("#8A8A8A", "#F6F6F6")}
INK = "#222222"; ARROW = "#444444"
W = 180.0; X0 = 11.0; BW = 166.0; PAD = 3.0          # band left edge, band width, inner padding
HEAD, BODY, PITCH = 6.8, 6.2, 2.95                   # font sizes (pt) and body line pitch (mm)
GAP = 6.0                                            # vertical gap between bands

def box_h(n_lines, extra=0.0): return 2.2 + 3.3 + n_lines * PITCH + 1.6 + extra

checks = []                                          # (text artist, x0, x1, y0, y1) for the fit check
def box(ax, x, y, w, h, head, lines, key, fc="white", ls="-", head_color=None, reserve=0.0):
    """Rounded box with a bold heading and centred body lines; the text block is centred vertically in the box
    (minus `reserve` mm kept free at the bottom for a result tag). A line may be (text, colour, weight)."""
    ec = C[key][0]; lines = [ln for ln in lines if ln]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.4", fc=fc, ec=ec, lw=0.7, ls=ls, zorder=2))
    yt = y + (h - reserve - (3.3 + len(lines) * PITCH - 0.6)) / 2
    t = ax.text(x + w / 2, yt, head, ha="center", va="top", fontsize=HEAD, fontweight="bold", color=head_color or INK, zorder=3)
    checks.append((t, x, x + w, y, y + h))
    for i, ln in enumerate(lines):
        txt, col, wt = (ln if isinstance(ln, tuple) else (ln, INK, "normal"))
        t = ax.text(x + w / 2, yt + 3.3 + i * PITCH, txt, ha="center", va="top", fontsize=BODY, color=col, fontweight=wt, zorder=3)
        checks.append((t, x, x + w, y, y + h))

def band(ax, y, h, num, title, key, dashed=False, tag=None):
    ec, tint = C[key]
    ax.add_patch(FancyBboxPatch((X0, y), BW, h, boxstyle="round,pad=0,rounding_size=2.0", fc=tint, ec=ec, lw=0.8 if dashed else 0.5,
                                ls=(0, (4, 2.5)) if dashed else "-", zorder=1))
    ax.text(X0 + PAD + 1.6, y + 4.3, str(num), ha="center", va="center", fontsize=7, fontweight="bold", color="white",
            bbox=dict(boxstyle="circle,pad=0.22", fc=ec, ec="none"), zorder=3)
    t = ax.text(X0 + PAD + 5.2, y + 4.3, title, ha="left", va="center", fontsize=8, fontweight="bold", color=ec, zorder=3)
    if tag:
        ax.figure.canvas.draw(); bb = t.get_window_extent().transformed(ax.transData.inverted())
        ax.text(max(bb.x0, bb.x1) + 3.4, y + 4.3, tag, ha="left", va="center", fontsize=6, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc=ec, ec="none"), zorder=3)

def arrow(ax, p0, p1, ls="-", color=ARROW):
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="-|>", lw=0.75, color=color, ls=ls, shrinkA=0, shrinkB=0, mutation_scale=6.5), zorder=4)

def row(ax, y, widths, gap, specs, key, n_lines, arrows=False, **kw):
    """Lay out boxes left to right inside the band; returns [(x, w)] and the common box height."""
    h = box_h(n_lines); x = X0 + PAD; pos = []
    for w, (head, lines) in zip(widths, specs):
        box(ax, x, y, w, h, head, lines, key, **kw); pos.append((x, w)); x += w + gap
    if arrows:
        for (xa, wa), (xb, _) in zip(pos[:-1], pos[1:]): arrow(ax, (xa + wa + 0.5, y + h / 2), (xb - 0.5, y + h / 2))
    return pos, h

# ---------------------------------------------------------------- geometry
TITLE = 8.2                                           # space taken by the band title row
h1, h2, h3, h4 = box_h(2), box_h(3), box_h(3), box_h(4)
h5a, h5b, BUS = box_h(2), box_h(4, extra=4.2), 8.0
bands_h = [TITLE + h + PAD for h in (h1, h2, h3, h4)] + [TITLE + h5a + BUS + h5b + PAD]
H = 2.5 + sum(bands_h) + 4 * GAP + 2.5
fig = plt.figure(figsize=(W / 25.4, H / 25.4)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.axis("off")
ys = [2.5]
for bh in bands_h[:-1]: ys.append(ys[-1] + bh + GAP)
IN = BW - 2 * PAD                                     # usable width inside a band (160 mm)
XC = X0 + BW / 2

# ---------------------------------------------------------------- 1 public data
band(ax, ys[0], bands_h[0], 1, "Public data", "data")
row(ax, ys[0] + TITLE, [(IN - 9) / 4] * 4, 3, [
    ("ENCODE knockdown RNA-seq", ["876 knockdown/control pairs", "478 RBPs, HepG2 + K562"]),
    ("ENCODE eCLIP", ["252 experiments", "168 RBPs"]),
    ("RBP descriptors", ["ESM2 protein embedding,", "GenePT text embedding, motifs"]),
    ("Event descriptors", ["sequence, splice-site strength,", "conservation, 5-mer spectra"])], "data", 2)

# ---------------------------------------------------------------- 2 reprocessing and artefact removal
band(ax, ys[1], bands_h[1], 2, "Junction-level reprocessing and artefact removal", "art")
fix = "#8A5A00"
pos2, _ = row(ax, ys[1] + TITLE, [(IN - 12) / 4] * 4, 4, [
    ("Junction-level PSI", ["rMATS-turbo on all BAM pairs", "546,917 events", "keyed by coordinates"]),
    ("Shared-control artefact", ["unrelated RBPs, same control:", "r 0.49–0.51", ("fix: control-centring", fix, "bold")]),
    ("Regression to the mean", ["control-deviation rule alone:", "direction accuracy 0.83 / 0.79", ("fix: clean target", fix, "bold")]),
    ("Measurement ceiling", ["same RBP, independent runs:", "direction agreement 0.56", "two thirds not reproducible"])], "art", 3, arrows=True)

# ---------------------------------------------------------------- 3 forward benchmark
band(ax, ys[2], bands_h[2], 3, "Forward benchmark on held-out RBPs", "fwd")
y3 = ys[2] + TITLE; w3 = [31, 31, 31, IN - 93 - 12]
pos3, _ = row(ax, y3, w3[:3], 4, [
    ("Held-out RBP split", ["paralog-grouped 5-fold", "BLASTP ≥ 40 % identity", "kept in one fold"]),
    ("Models", ["ridge regression", "deep bilinear (v3)", "1D-CNN (v4)"]),
    ("Controls", ["shuffled RBP identity", "event-only model", ""])], "fwd", 3, arrows=True)
xr = pos3[-1][0] + pos3[-1][1] + 4
box(ax, xr, y3, w3[3], h3, "Result", ["direction accuracy 0.51–0.52; linear = deep", "eCLIP / motif add < 0.7 points",
                                      "core spliceosome is the one predictable class"], "out", fc=C["out"][1])
arrow(ax, (xr - 3.5, y3 + h3 / 2), (xr - 0.5, y3 + h3 / 2))

# ---------------------------------------------------------------- 4 programme-level upgrade (in progress)
band(ax, ys[3], bands_h[3], 4, "Programme-level upgrade", "prog", dashed=True, tag="in progress")
row(ax, ys[3] + TITLE, [60, 52, IN - 112 - 6], 3, [
    ("a  Programme-level targets", ["10 response programmes (PCA on training RBPs)", "instead of ~50,000 single events",
                                    "cross-cell-line reproducibility: r 0.52", "(programme 1) vs 0.05 (single events)"]),
    ("b  Functional RBP descriptors", ["genome-wide Perturb-seq (K562)", "STRING network / protein complexes", "DepMap co-dependency", ""]),
    ("c  Observational prior", ["RBP expression × event PSI", "co-variation across", "TCGA tumours", ""])], "prog", 4, ls=(0, (3, 2)))

# ---------------------------------------------------------------- 5 inverse mapping
band(ax, ys[4], bands_h[4], 5, "Inverse mapping: from a splicing shift to the perturbed RBP", "inv")
y5 = ys[4] + TITLE
pos5, _ = row(ax, y5, [(IN - 12) / 3] * 3, 6, [
    ("Dictionary", ["measured control-centred responses", "of 478 RBPs"]),
    ("Query", ["splicing shift of a sample group", "(case vs control)"]),
    ("Scoring", ["Pearson r ranking over shared events", "+ label-permutation null"])], "inv", 2, arrows=True)
y5b = y5 + h5a + BUS; wb = (IN - 9) / 4
outs = [("ENCODE retrieval", ["across cell lines:", "top-5 25 % vs 1 % by chance", "high-reproducibility RBPs:", "median rank 1"], "retrieved", C["inv"][0]),
        ("HeyA8 LARP1 knockout", ["in-house, post-freeze", "class-level hit:", "PCBP2, SF3B1, SF3A1", ""], "class-level hit", C["inv"][0]),
        ("TCGA mutant tumours", ["52 tests, negative overall", "UVM SF3B1 hotspot retrieves", "the core-spliceosome class", ""], "one class-level hit", C["inv"][0]),
        ("Premalignant cohorts", ["fallopian tube, breast, sorted", "mammary epithelium: negative", "after composition / batch", "adjustment"], "negative", C["out"][0])]
ybus = y5 + h5a + BUS / 2; xs = []
for i, (head, lines, tag, tc) in enumerate(outs):
    x = X0 + PAD + i * (wb + 3); xs.append(x + wb / 2)
    box(ax, x, y5b, wb, h5b, head, lines, "out", fc=C["out"][1], reserve=4.6)
    ax.text(x + wb / 2, y5b + h5b - 3.1, tag, ha="center", va="center", fontsize=6, fontweight="bold", color="white",
            bbox=dict(boxstyle="round,pad=0.45,rounding_size=0.9", fc=tc, ec="none"), zorder=3)
    arrow(ax, (x + wb / 2, ybus), (x + wb / 2, y5b - 0.5))
xs_sc = pos5[2][0] + pos5[2][1] / 2                   # centre of the scoring box
ax.plot([xs[0], xs[-1]], [ybus, ybus], color=ARROW, lw=0.75, zorder=4, solid_capstyle="round")
ax.plot([xs_sc, xs_sc], [y5 + h5a + 0.5, ybus], color=ARROW, lw=0.75, zorder=4)

# ---------------------------------------------------------------- arrows between stages
for i in range(3): arrow(ax, (XC, ys[i] + bands_h[i] + 0.6), (XC, ys[i + 1] - 0.6))
ax.text(XC + 2.5, ys[2] + bands_h[2] + GAP / 2, "bottleneck: noisy event-level targets, uninformative RBP descriptors", ha="left", va="center", fontsize=6, style="italic", color="#555555")
arrow(ax, (XC, ys[3] + bands_h[3] + 0.6), (XC, ys[4] - 0.6), ls=(0, (3, 2)), color=C["prog"][0])
ax.text(XC + 2.5, ys[3] + bands_h[3] + GAP / 2, "programme activity scores (planned)", ha="left", va="center", fontsize=6, style="italic", color=C["prog"][0])
# bypass: the dictionary is built from the measured, control-centred responses of stage 2, not from model predictions
ya = ys[1] + TITLE + h2 / 2; yb = y5 + h5a / 2; xl = 5.5
ax.plot([X0, xl, xl], [ya, ya, yb], color=C["inv"][0], lw=0.75, zorder=4, solid_joinstyle="round")
arrow(ax, (xl, yb), (X0 + PAD - 0.5, yb), color=C["inv"][0])
ax.text(xl - 2.6, (ya + yb) / 2, "measured, control-centred responses", rotation=90, ha="center", va="center", fontsize=6, style="italic", color=C["inv"][0])

# ---------------------------------------------------------------- fit check, save
fig.canvas.draw(); inv = ax.transData.inverted(); bad = 0
for t, x0, x1, y0, y1 in checks:
    if not t.get_text(): continue
    bb = t.get_window_extent().transformed(inv); tx0, tx1 = sorted((bb.x0, bb.x1)); ty0, ty1 = sorted((bb.y0, bb.y1))
    if tx0 < x0 + 0.8 or tx1 > x1 - 0.8 or ty0 < y0 + 0.5 or ty1 > y1 - 0.5:
        bad += 1; print(f"FIT WARNING: '{t.get_text()}' text x {tx0:.1f}-{tx1:.1f} box x {x0:.1f}-{x1:.1f}; text y {ty0:.1f}-{ty1:.1f} box y {y0:.1f}-{y1:.1f}")
print(f"canvas {W:.0f} x {H:.1f} mm; labels checked {len(checks)}; fit warnings {bad}")
fig.savefig(f"{O}/fig_workflow.png"); fig.savefig(f"{O}/fig_workflow.pdf"); print("fig", f"{O}/fig_workflow.png")
