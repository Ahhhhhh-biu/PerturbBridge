"""Gene-expression output layer of the forward model: gene encoder x perturbation encoder with a low-rank bilinear interaction,
trained on control-centred delta log2 TPM (13_model/tensors_expr_<cell>.npz), evaluated on RBP-held-out folds (same folds as splicing).
  Gn(g)  = MLP([gene features, proj GenePT(64), learned gene embedding(32)])                       -> d
  P(r)   = MLP([proj ESM2+GenePT of the RBP (PCA 32 -> 64, dropout), gene-level eCLIP(2) + flag, assay])  -> d
  h      = Gn * (W P) + assay bias ; y_hat = MLP(h)   (linear output, log2 units)
Tiers: A = all RBP inputs; B = no eCLIP. --rbp_mode none (gene-only model) / shuffle (fixed RBP permutation) are the controls.
Baselines reported per run without training: mean baseline (mean Y of the training runs per gene = Ahlmann-Eltze style linear baseline).
Usage: python 11_train_expr.py <cell_line> <fold> <tier A|B> [--epochs 20 --batch 4096 --lr 3e-4 --sigw 10 --rpca 32 --rdrop 0.3 --cos --seed 0 --tag _x --rbp_mode normal]
Outputs: 13_model/runs_expr/<cell_line>_f<fold>_<tier><tag>/ {per_run.tsv, pred_heldout.npz, model.pt, log.txt}"""
import sys, os, time, argparse, numpy as np, pandas as pd, torch, torch.nn as nn
ap = argparse.ArgumentParser(); ap.add_argument("cell_line"); ap.add_argument("fold", type=int); ap.add_argument("tier"); ap.add_argument("--epochs", type=int, default=20)
ap.add_argument("--device", default="auto"); ap.add_argument("--d", type=int, default=128); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--batch", type=int, default=4096)
ap.add_argument("--lr", type=float, default=3e-4); ap.add_argument("--tag", default=""); ap.add_argument("--sigw", type=float, default=10.0); ap.add_argument("--rpca", type=int, default=32)
ap.add_argument("--rdrop", type=float, default=0.3); ap.add_argument("--cos", action="store_true"); ap.add_argument("--wd", type=float, default=1e-4); ap.add_argument("--huber", type=float, default=0.5)
ap.add_argument("--rbp_mode", default="normal", choices=["normal", "none", "shuffle"]); ap.add_argument("--test", action="store_true")
a = ap.parse_args(); PB = "/well/blagden/users/eic789/PerturbBridge"; torch.manual_seed(a.seed); np.random.seed(a.seed)
dev = torch.device("cuda" if a.device == "auto" and torch.cuda.is_available() else ("cpu" if a.device == "auto" else a.device))
O = f"{PB}/13_model/runs_expr/{a.cell_line}_f{a.fold}_{a.tier}{a.tag}" + ("_test" if a.test else ""); os.makedirs(O, exist_ok=True); log = open(f"{O}/log.txt", "a")
def P(*x): print(*x, flush=True); print(*x, file=log, flush=True)
z = np.load(f"{PB}/13_model/tensors_expr_{a.cell_line}.npz", allow_pickle=True)
Y = z["Y"]; SIG = z["SIG"]; IS_T = z["IS_T"]; nR, nG = Y.shape
Xg = torch.tensor(z["Xg"]); Xg[:, :2] = (Xg[:, :2] - Xg[:, :2].mean(0)) / (Xg[:, :2].std(0) + 1e-6); Ggp = torch.tensor(z["G_gp"]); Ggp = (Ggp - Ggp.mean(0)) / (Ggp.std(0) + 1e-6)
R = torch.tensor(np.hstack([z["R_esm"], z["R_gp"]])); R = (R - R.mean(0)) / (R.std(0) + 1e-6)
if a.rpca: U_, S_, V_ = torch.linalg.svd(R, full_matrices=False); R = R @ V_[:a.rpca].T; R = (R - R.mean(0)) / (R.std(0) + 1e-6)
EC = torch.tensor(z["EC"]); tei = torch.tensor(z["target_ecidx"]); run_target = torch.tensor(z["run_target"]); run_assay = torch.tensor(z["run_assay"]); run_fold = z["run_fold"]; run_tier = z["run_tier"]
tperm = torch.randperm(len(z["targets"]), generator=torch.Generator().manual_seed(1000 + a.seed)) if a.rbp_mode == "shuffle" else torch.arange(len(z["targets"]))
rr, gg = np.meshgrid(np.arange(nR), np.arange(nG), indexing="ij"); keep = ~IS_T.ravel()
if a.test: keep &= gg.ravel() < 1500; a.epochs = min(a.epochs, 2)
row_r = rr.ravel()[keep].astype(np.int32); row_g = gg.ravel()[keep].astype(np.int32); row_y = Y.ravel()[keep]; row_sig = SIG.ravel()[keep]
rf = run_fold[row_r]; te = rf == a.fold; va = rf == (a.fold + 1) % 5; tr = ~te & ~va
train_runs = np.unique(row_r[tr]); mean_base = np.nanmean(np.where(IS_T[train_runs], np.nan, Y[train_runs]), 0)   # per-gene mean response of the training runs
P(f"{a.cell_line} fold {a.fold} tier {a.tier} mode {a.rbp_mode} device {dev} | rows train {tr.sum():,} val {va.sum():,} test {te.sum():,} | genes {nG:,} runs {nR} | sig frac {row_sig.mean():.4f}")
Xg, Ggp, R, EC, run_target, run_assay, tei, tperm = [x.to(dev) for x in (Xg, Ggp, R, EC, run_target, run_assay, tei, tperm)]
class Net(nn.Module):
    def __init__(s, d):
        super().__init__(); s.gproj = nn.Linear(1536, 64); s.gemb = nn.Embedding(nG, 32); s.G = nn.Sequential(nn.Linear(Xg.shape[1] + 64 + 32, 256), nn.GELU(), nn.Linear(256, d))
        s.rproj = nn.Sequential(nn.Linear(R.shape[1], 64), nn.Dropout(a.rdrop)); s.Pn = nn.Sequential(nn.Linear(64 + 2 + 1 + 1, 256), nn.GELU(), nn.Linear(256, d)); s.W = nn.Linear(d, d, bias=False)
        s.head = nn.Sequential(nn.Linear(d, 128), nn.GELU(), nn.Linear(128, 1)); s.abias = nn.Embedding(2, d)
    def forward(s, g, r):
        t = tperm[run_target[r]]; G = s.G(torch.cat([Xg[g], s.gproj(Ggp[g]), s.gemb(g)], 1))
        ei = tei[t]; ec = torch.where((ei >= 0)[:, None], EC[ei.clamp(min=0), g], torch.zeros(len(g), 2, device=dev)); hase = (ei >= 0).float()[:, None]
        if a.tier == "B": ec = ec * 0; hase = hase * 0
        rp = s.rproj(R[t])
        if a.rbp_mode == "none": rp, ec, hase = rp * 0, ec * 0, hase * 0
        Pv = s.Pn(torch.cat([rp, ec, hase, run_assay[r][:, None]], 1)); h = G * s.W(Pv) + s.abias(run_assay[r].long()); return s.head(h)[:, 0]
net = Net(a.d).to(dev); opt = torch.optim.AdamW(net.parameters(), lr=a.lr, weight_decay=a.wd); loss_fn = nn.HuberLoss(delta=a.huber, reduction="none")
nsteps = a.epochs * ((tr.sum() + a.batch - 1) // a.batch); sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=int(nsteps), eta_min=a.lr / 10) if a.cos else None
def tens(idx): return torch.tensor(row_g[idx]).to(dev), torch.tensor(row_r[idx]).to(dev), torch.tensor(row_y[idx]).to(dev)
def wts(idx): return torch.tensor(1.0 + (a.sigw - 1.0) * row_sig[idx].astype(np.float32)).to(dev)
def predict(mask):
    net.eval(); idx = np.where(mask)[0]; out = np.zeros(len(idx), dtype=np.float32)
    with torch.no_grad():
        for i in range(0, len(idx), 65536): g, r, _ = tens(idx[i:i + 65536]); out[i:i + 65536] = net(g, r).cpu().numpy()
    net.train(); return idx, out
def corr(x, y): return np.corrcoef(x, y)[0, 1] if len(x) > 2 and np.std(x) > 0 and np.std(y) > 0 else np.nan
def run_metrics(idx, pred):
    d = pd.DataFrame({"r": row_r[idx], "g": row_g[idx], "y": row_y[idx], "p": pred, "sig": row_sig[idx]}); out = []
    for r, g in d.groupby("r"):
        s = g[g.sig]; mb = mean_base[g.g.values]; mbs = mean_base[s.g.values]
        out.append(dict(run=z["runs"][r], target=z["targets"][run_target[r].item()], tier_data="A" if run_tier[r] else "B", n=len(g), n_sig=len(s),
                        r_all=corr(g.y, g.p), r_sig=corr(s.y, s.p), sign_sig=(np.sign(s.y) == np.sign(s.p)).mean() if len(s) else np.nan, mse_all=float(((g.y - g.p) ** 2).mean()),
                        r_all_meanbase=corr(g.y, mb), r_sig_meanbase=corr(s.y, mbs), sign_sig_meanbase=(np.sign(s.y) == np.sign(mbs)).mean() if len(s) else np.nan,
                        r_all_resid=corr(g.y - mb, g.p - mb), mean_abs_pred=float(np.abs(g.p).mean()), mean_abs_y_sig=float(np.abs(s.y).mean()) if len(s) else np.nan))
    return pd.DataFrame(out)
tr_idx = np.where(tr)[0]; best = (-1e9, None, -1); t0 = time.time()
for ep in range(a.epochs):
    perm = np.random.permutation(tr_idx); tot = 0.0
    for i in range(0, len(perm), a.batch):
        g, r, y = tens(perm[i:i + a.batch]); w = wts(perm[i:i + a.batch]); loss = (loss_fn(net(g, r), y) * w).sum() / w.sum(); opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item() * len(g)
        if sched: sched.step()
    vi, vp = predict(va); vm = run_metrics(vi, vp); vb = vm.groupby("target").mean(numeric_only=True); score = float(np.nanmean(vb.r_all))
    P(f"epoch {ep} train huber {tot / len(perm):.5f} | val by RBP: r_all {score:.4f} r_sig {np.nanmean(vb.r_sig):.4f} sign_sig {np.nanmean(vb.sign_sig):.4f} | meanbase r_all {np.nanmean(vb.r_all_meanbase):.4f} | {time.time() - t0:.0f}s")
    if score > best[0]: best = (score, {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}, ep)
net.load_state_dict(best[1]); P(f"best epoch {best[2]} val r_all {best[0]:.4f}")
ti, tp = predict(te); tm = run_metrics(ti, tp); tm["fold"] = a.fold; tm["tier_model"] = a.tier; tm["rbp_mode"] = a.rbp_mode; tm.to_csv(f"{O}/per_run.tsv", sep="\t", index=False)
np.savez_compressed(f"{O}/pred_heldout.npz", idx=ti, pred=tp); torch.save(best[1], f"{O}/model.pt")
by = tm.groupby("target").mean(numeric_only=True)
P("held-out (by RBP):", {c: round(float(np.nanmean(by[c])), 4) for c in ("r_all", "r_sig", "sign_sig", "r_all_meanbase", "r_sig_meanbase", "sign_sig_meanbase", "r_all_resid")}); P("DONE", O)
