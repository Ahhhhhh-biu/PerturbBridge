"""Forward model (design 3): event encoder x perturbation encoder with a low-rank bilinear interaction, trained on
control-centred dPSI, evaluated on RBP-held-out folds.
  E(e)   = MLP([event features, baseline PSI, projected 5-mer spectrum])                       -> d
  P(r,e) = MLP([proj ESM2, proj GenePT, motif density(3) + flag, eCLIP(6) + flag, assay])       -> d
  h      = E * (W P) + cell/assay bias ; dPSI = tanh(MLP(h))
Tiers are input masks applied to the perturbation encoder (A: all; B: no eCLIP; C: no eCLIP, no motif), at train and test time.
Usage: python 01_train.py <cell_line> <fold> <tier A|B|C> [--epochs N] [--device auto|cpu|cuda] [--test]
  training options added 2026-11-07 (all default to the first-round behaviour): --sel weighted (model selection by n_sig-weighted
  RBP-mean r_sig), --sigw W (loss weight of significant rows), --rpca K (PCA-compress the RBP embedding), --rdrop p, --cos
Outputs: 13_model/runs/<cell_line>_f<fold>_<tier>/ {per_run.tsv, summary.tsv, pred_heldout.npz, model.pt, log.txt}"""
import sys, os, time, argparse, numpy as np, pandas as pd, torch, torch.nn as nn
ap = argparse.ArgumentParser(); ap.add_argument("cell_line"); ap.add_argument("fold", type=int); ap.add_argument("tier"); ap.add_argument("--epochs", type=int, default=12)
ap.add_argument("--device", default="auto"); ap.add_argument("--test", action="store_true"); ap.add_argument("--d", type=int, default=128); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--batch", type=int, default=8192); ap.add_argument("--lr", type=float, default=1e-3); ap.add_argument("--tag", default="")
ap.add_argument("--sel", default="rbp", choices=["rbp", "weighted"], help="model selection: RBP-mean r_sig, or n_sig-weighted RBP mean")
ap.add_argument("--sigw", type=float, default=1.0, help="loss weight of significant rows (FDR<0.05 & |dPSI|>0.1) relative to the rest")
ap.add_argument("--rpca", type=int, default=0, help="PCA-compress ESM2+GenePT to this many dims (fitted on all targets, unsupervised); 0 = raw 2816")
ap.add_argument("--rdrop", type=float, default=0.0, help="dropout on the projected RBP embedding")
ap.add_argument("--cos", action="store_true", help="cosine learning-rate decay to lr/10 over all steps")
ap.add_argument("--wd", type=float, default=1e-4, help="AdamW weight decay")
ap.add_argument("--target", default="fdr", choices=["fdr", "clean"], help="which rows count as significant for loss weighting, model selection and metrics: fdr = FDR<0.05 & |dPSI|>0.1 vs own control; clean = replicate-consistent deviation from the control-group consensus (13_model/clean_flags_<cell>.npz)")
ap.add_argument("--no_psi", action="store_true", help="replace the run-specific control PSI input by the event's mean control PSI over all runs (removes the control-deviation shortcut)")
ap.add_argument("--rbp_mode", default="normal", choices=["normal", "none", "shuffle"], help="control: none = zero all RBP inputs (event-only model); shuffle = each RBP gets another RBP's features (fixed permutation)")
a = ap.parse_args(); PB = "/well/blagden/users/eic789/PerturbBridge"; torch.manual_seed(a.seed); np.random.seed(a.seed)
dev = torch.device("cuda" if a.device == "auto" and torch.cuda.is_available() else ("cpu" if a.device == "auto" else a.device))
O = f"{PB}/13_model/runs/{a.cell_line}_f{a.fold}_{a.tier}{a.tag}" + ("_test" if a.test else ""); os.makedirs(O, exist_ok=True); log = open(f"{O}/log.txt", "a")
def P(*x): print(*x, flush=True); print(*x, file=log, flush=True)
z = np.load(f"{PB}/13_model/tensors_{a.cell_line}.npz", allow_pickle=True)
Xe = torch.tensor(z["Xe"]); K = torch.tensor(z["K"]).reshape(len(Xe), -1); EC = torch.tensor(z["EC"]); tei = torch.tensor(z["target_ecidx"])
R = torch.tensor(np.hstack([z["R_esm"], z["R_gp"]])); R = (R - R.mean(0)) / (R.std(0) + 1e-6)
if a.rpca: U_, S_, V_ = torch.linalg.svd(R, full_matrices=False); R = R @ V_[:a.rpca].T; R = (R - R.mean(0)) / (R.std(0) + 1e-6)
Rmot = torch.tensor(z["R_mot"]); Rhas = torch.tensor(z["R_has_mot"])
run_target = torch.tensor(z["run_target"]); run_assay = torch.tensor(z["run_assay"]); run_fold = z["run_fold"]; run_tier = z["run_tier"]
tperm = torch.randperm(len(z["targets"]), generator=torch.Generator().manual_seed(1000 + a.seed)) if a.rbp_mode == "shuffle" else torch.arange(len(z["targets"]))
rows = {k: z[k] for k in ("row_e", "row_r", "row_psi", "row_y", "row_yraw", "row_sig")}
if a.target == "clean": rows["row_sig"] = np.load(f"{PB}/13_model/clean_flags_{a.cell_line}.npz")["row_clean"]
if a.no_psi: _pm = np.bincount(rows["row_e"], weights=rows["row_psi"]) / np.maximum(np.bincount(rows["row_e"]), 1); rows["row_psi"] = _pm[rows["row_e"]].astype(np.float32)
if a.test:
    keep = rows["row_e"] < 2000; rows = {k: v[keep] for k, v in rows.items()}; a.epochs = min(a.epochs, 2)
# standardise continuous event features (columns not one-hot / flags)
cont = torch.tensor([not (n in ("SE", "A3SS", "A5SS", "MXE", "RI") or n.endswith("_missing")) for n in z["Xe_names"]]); mu, sd = Xe[:, cont].mean(0), Xe[:, cont].std(0) + 1e-6; Xe[:, cont] = (Xe[:, cont] - mu) / sd
Kn = (K - K.mean(0)) / (K.std(0) + 1e-6)
# motif density per (event, target) computed on the fly: Kraw[e] (3x1024) . Rmot[t] (1024) -> 3
Kraw = torch.tensor(z["K"]); 
# split: held-out fold = test runs; inner validation = runs of fold (fold+1)%5 among training; rest = train
rf = run_fold[rows["row_r"]]; te = rf == a.fold; va = rf == (a.fold + 1) % 5; tr = ~te & ~va
P(f"{a.cell_line} fold {a.fold} tier {a.tier} device {dev} | rows train {tr.sum():,} val {va.sum():,} test {te.sum():,} | events {len(Xe):,} runs {len(run_target)}")
Xe, Kn, Kraw, EC, R, Rmot, Rhas, run_target, run_assay, tei, tperm = [x.to(dev) for x in (Xe, Kn, Kraw, EC, R, Rmot, Rhas, run_target, run_assay, tei, tperm)]
class Net(nn.Module):
    def __init__(s, d):
        super().__init__(); s.kproj = nn.Linear(3 * 1024, 64); s.E = nn.Sequential(nn.Linear(Xe.shape[1] + 1 + 64, 256), nn.GELU(), nn.Linear(256, d))
        s.rproj = nn.Sequential(nn.Linear(R.shape[1], 64), nn.Dropout(a.rdrop)); s.Pn = nn.Sequential(nn.Linear(64 + 3 + 1 + 6 + 1 + 1, 256), nn.GELU(), nn.Linear(256, d)); s.W = nn.Linear(d, d, bias=False)
        s.head = nn.Sequential(nn.Linear(d, 128), nn.GELU(), nn.Linear(128, 1)); s.abias = nn.Embedding(2, d)
    def forward(s, e, r, psi):
        t = tperm[run_target[r]]; E = s.E(torch.cat([Xe[e], psi[:, None], s.kproj(Kn[e])], 1))
        mot = torch.einsum("bkn,bn->bk", Kraw[e], Rmot[t]) * Rhas[t][:, None]; hasm = Rhas[t][:, None]
        ei = tei[t]; ec = torch.where((ei >= 0)[:, None], EC[ei.clamp(min=0), e], torch.zeros(len(e), 6, device=dev)); hase = (ei >= 0).float()[:, None]
        if a.tier in ("B", "C"): ec = ec * 0; hase = hase * 0
        if a.tier == "C": mot = mot * 0; hasm = hasm * 0
        rp = s.rproj(R[t])
        if a.rbp_mode == "none": rp, mot, hasm, ec, hase = rp * 0, mot * 0, hasm * 0, ec * 0, hase * 0
        Pv = s.Pn(torch.cat([rp, mot, hasm, ec, hase, run_assay[r][:, None]], 1))
        h = E * s.W(Pv) + s.abias(run_assay[r].long()); return torch.tanh(s.head(h))[:, 0]
net = Net(a.d).to(dev); opt = torch.optim.AdamW(net.parameters(), lr=a.lr, weight_decay=a.wd); loss_fn = nn.HuberLoss(delta=0.1, reduction="none")
nsteps = a.epochs * ((tr.sum() + a.batch - 1) // a.batch); sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=int(nsteps), eta_min=a.lr / 10) if a.cos else None
def tens(idx): return (torch.tensor(rows["row_e"][idx]).to(dev), torch.tensor(rows["row_r"][idx]).to(dev), torch.tensor(rows["row_psi"][idx]).to(dev), torch.tensor(rows["row_y"][idx]).to(dev))
def wts(idx): return torch.tensor(1.0 + (a.sigw - 1.0) * rows["row_sig"][idx].astype(np.float32)).to(dev)
def predict(mask):
    net.eval(); idx = np.where(mask)[0]; out = np.zeros(len(idx), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(idx), 65536): e, r, p, _ = tens(idx[i:i + 65536]); out[i:i + 65536] = net(e, r, p).cpu().numpy()
    net.train(); return idx, out
def run_metrics(idx, pred):
    d = pd.DataFrame({"r": rows["row_r"][idx], "y": rows["row_y"][idx], "p": pred, "sig": rows["row_sig"][idx]}); out = []
    for r, g in d.groupby("r"):
        s = g[g.sig]; m = dict(run=z["runs"][r], target=z["targets"][run_target[r].item()], tier_data="A" if run_tier[r] else "B", n=len(g), n_sig=len(s))
        m["r_all"] = np.corrcoef(g.y, g.p)[0, 1] if g.p.std() > 0 else np.nan; m["r_sig"] = np.corrcoef(s.y, s.p)[0, 1] if len(s) > 2 and s.p.std() > 0 else np.nan
        m["sign_sig"] = (np.sign(s.y) == np.sign(s.p)).mean() if len(s) else np.nan; m["mse_all"] = float(((g.y - g.p) ** 2).mean())
        for lo, hi in [(0.05, 0.1), (0.1, 0.2), (0.2, 1.01)]:
            b = g[(g.y.abs() >= lo) & (g.y.abs() < hi)]; m[f"sign_{lo}_{hi}"] = (np.sign(b.y) == np.sign(b.p)).mean() if len(b) else np.nan
        out.append(m)
    return pd.DataFrame(out)
tr_idx = np.where(tr)[0]; best = (-1e9, None); t0 = time.time()
for ep in range(a.epochs):
    perm = np.random.permutation(tr_idx); tot = 0.0
    for i in range(0, len(perm), a.batch):
        e, r, p, y = tens(perm[i:i + a.batch]); w = wts(perm[i:i + a.batch]); loss = (loss_fn(net(e, r, p), y) * w).sum() / w.sum(); opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item() * len(e)
        if sched: sched.step()
    vi, vp = predict(va); vm = run_metrics(vi, vp); vb = vm.groupby("target").agg(r_sig=("r_sig", "mean"), n_sig=("n_sig", "sum")).dropna()
    score = np.average(vb.r_sig, weights=vb.n_sig) if a.sel == "weighted" else np.nanmean(vm.groupby("target").r_sig.mean())
    P(f"epoch {ep} train huber {tot / len(perm):.5f} | val r_sig(by RBP) {score:.4f} sign_sig {np.nanmean(vm.groupby('target').sign_sig.mean()):.4f} r_all {np.nanmean(vm.groupby('target').r_all.mean()):.4f} | {time.time() - t0:.0f}s")
    if score > best[0]: best = (score, {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}, ep)
net.load_state_dict(best[1]); P(f"best epoch {best[2]} val r_sig {best[0]:.4f}")
ti, tp = predict(te); tm = run_metrics(ti, tp); tm["fold"] = a.fold; tm["tier_model"] = a.tier; tm.to_csv(f"{O}/per_run.tsv", sep="\t", index=False)
np.savez_compressed(f"{O}/pred_heldout.npz", idx=ti, pred=tp); torch.save(best[1], f"{O}/model.pt")
by = tm.groupby("target").mean(numeric_only=True); S = by.agg(["mean", "median", "count"]).T; S.to_csv(f"{O}/summary.tsv", sep="\t")
P("held-out (by RBP):", {c: round(float(by[c].mean()), 4) for c in ("r_all", "r_sig", "sign_sig", "sign_0.1_0.2", "sign_0.2_1.01")}, "| tier-A RBPs:", {c: round(float(by[tm.groupby('target').tier_data.first() == 'A'][c].mean()), 4) for c in ("r_sig", "sign_sig")})
P("DONE", O)
