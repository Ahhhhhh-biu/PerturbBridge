"""Shared style for the manuscript figures in this folder (Arial, 300 DPI, editable PDF text)."""
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; from matplotlib import font_manager
for f in ("Arial.ttf", "Arial-Bold.ttf"): font_manager.fontManager.addfont(f"/well/blagden/users/eic789/.fonts/{f}")
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                     "savefig.dpi": 300, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/15_manuscript/figures"; CELLS = ("HepG2", "K562"); COL = {"HepG2": "#2f5d8a", "K562": "#8fb4d9"}
def panel(ax, letter, title): ax.set_title(f"{letter}  {title}", loc="left", fontweight="bold", fontsize=7.5)
def save(fig, name): fig.savefig(f"{OUT}/{name}.png"); fig.savefig(f"{OUT}/{name}.pdf"); print("saved", f"{OUT}/{name}.png")
