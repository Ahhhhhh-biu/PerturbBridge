"""Linear (ridge) baseline on exactly the same folds, rows, target and loss weighting as v3: y = w . [event block, RBP block, event x RBP
interaction block] with event block = standardised event features (23) + mean control PSI + 32 PCs of the 5-mer spectrum, RBP block = 32 PCs of
ESM2+GenePT + motif density (3) + flag + eCLIP (6) + flag + assay, interaction block = outer product of the first 8 event PCs and first 8 RBP PCs (64).
Closed-form weighted ridge (weights 100 on clean rows, as in v3), lambda chosen on the inner validation fold by n_sig-weighted RBP-mean r_sig.
Usage: python 02_ridge_baseline.py <cell_line> <fold> [--rbp_mode normal|shuffle] -> 13_model/runs/<cell>_f<fold>_A_ridge_s0/ {per_run.tsv, pred_heldout.npz, log.txt}"""
import sys, os, argparse, numpy as np, pandas as pd, torch
ap = argparse.ArgumentParser(); ap.add_argument("cell_line"); ap.add_argument("fold", type=int); ap.add_argument("--rbp_mode", default="normal", choices=["normal", "shuffle"]); ap.add_argument("--sigw", type=float, default=100.0); a = ap.parse_args()
PB = "/well/blagden/users/eic789/PerturbBridge"; dev = torch.device("cuda" if torch.cuda.is_available() else "cpu"); tag = "_ridge_s0" if a.rbp_mode == "normal" else "_ridge_ctrl_shuf_s0"
O = f"{PB}/13_model/runs/{a.cell_line}_f{a.fold}_A{tag}"; os.makedirs(O, exist_ok=True); log = open(f"{O}/log.txt", "w")
def P(*x): print(*x, flush=True); print(*x, file=log, flush=True)
z = np.load(f"{PB}/13_model/tensors_{a.cell_line}.npz", allow_pickle=True)
Xe = torch.tensor(z["Xe"]); names = z["Xe_names"]; cont = torch.tensor([not (n in ("SE", "A3SS", "A5SS", "MXE", "RI") or n.endswith("_missing")) for n in names]); Xe[:, cont] = (Xe[:, cont] - Xe[:, cont].mean(0)) / (Xe[:, cont].std(0) + 1e-6)
K = torch.tensor(z["K"]).reshape(len(Xe), -1); Kn = (K - K.mean(0)) / (K.std(0) + 1e-6); U_, S_, V_ = torch.linalg.svd(Kn, full_matrices=False); Kp = Kn @ V_[:32].T; Kp = (Kp - Kp.mean(0)) / (Kp.std(0) + 1e-6)
R = torch.tensor(np.hstack([z["R_esm"], z["R_gp"]])); R = (R - R.mean(0)) / (R.std(0) + 1e-6); U_, S_, V_ = torch.linalg.svd(R, full_matrices=False); Rp = R @ V_[:32].T; Rp = (Rp - Rp.mean(0)) / (Rp.std(0) + 1e-6)
Rmot = torch.tensor(z["R_mot"]); Rhas = torch.tensor(z["R_has_mot"]); EC = torch.tensor(z["EC"]); tei = torch.tensor(z["target_ecidx"]); Kraw = torch.tensor(z["K"])
run_target = torch.tensor(z["run_target"]); run_assay = torch.tensor(z["run_assay"]); run_fold = z["run_fold"]; run_tier = z["run_tier"]
tperm = torch.randperm(len(z["targets"]), generator=torch.Generator().manual_seed(1000)) if a.rbp_mode == "shuffle" else torch.arange(len(z["targets"]))
rows = {k: z[k] for k in ("row_e", "row_r", "row_psi", "row_y")}; rows["row_sig"] = np.load(f"{PB}/13_model/clean_flags_{a.cell_line}.npz")["row_clean"]
pm = np.bincount(rows["row_e"], weights=rows["row_psi"]) / np.maximum(np.bincount(rows["row_e"]), 1); rows["row_psi"] = pm[rows["row_e"]].astype(np.float32)
rf = run_fold[rows["row_r"]]; te = rf == a.fold; va = rf == (a.fold + 1) % 5; tr = ~te & ~va
Xe, Kp, Kraw, Rp, Rmot, Rhas, EC, tei, run_target, run_assay, tperm = [x.to(dev) for x in (Xe, Kp, Kraw, Rp, Rmot, Rhas, EC, tei, run_target, run_assay, tperm)]
def feats(idx):
    e = torch.tensor(rows["row_e"][idx]).to(dev); r = torch.tensor(rows["row_r"][idx]).to(dev); t = tperm[run_target[r]]; psi = torch.tensor(rows["row_psi"][idx]).to(dev)[:, None]
    mot = torch.einsum("bkn,bn->bk", Kraw[e], Rmot[t]) * Rhas[t][:, None]; ei = tei[t]; ec = torch.where((ei >= 0)[:, None], EC[ei.clamp(min=0), e], torch.zeros(len(e), 6, device=dev)); hase = (ei >= 0).float()[:, None]
    Eb = torch.cat([Xe[e], psi, Kp[e]], 1); Rb = torch.cat([Rp[t], mot, Rhas[t][:, None], ec, hase, run_assay[r][:, None]], 1); inter = (Kp[e][:, :8, None] * Rp[t][:, None, :8]).reshape(len(e), 64)
    return torch.cat([torch.ones(len(e), 1, device=dev), Eb, Rb, inter], 1)
def wts(idx): return torch.tensor(1.0 + (a.sigw - 1.0) * rows["row_sig"][idx].astype(np.float32)).to(dev)
tr_idx = np.where(tr)[0]; D = feats(tr_idx[:10]).shape[1]; XtX = torch.zeros(D, D, device=dev, dtype=torch.float64); Xty = torch.zeros(D, device=dev, dtype=torch.float64)
for i in range(0, len(tr_idx), 200000):
    idx = tr_idx[i:i + 200000]; X = feats(idx).double(); w = wts(idx).double(); y = torch.tensor(rows["row_y"][idx]).to(dev).double(); XtX += X.T @ (X * w[:, None]); Xty += X.T @ (w * y)
P(f"{a.cell_line} fold {a.fold} ridge {a.rbp_mode}: train rows {len(tr_idx):,}, features {D}")
def predict(beta, mask):
    idx = np.where(mask)[0]; out = np.zeros(len(idx), dtype=np.float32)
    for i in range(0, len(idx), 200000): out[i:i + 200000] = (feats(idx[i:i + 200000]).double() @ beta).float().cpu().numpy()
    return idx, out
def run_metrics(idx, pred):
    d = pd.DataFrame({"r": rows["row_r"][idx], "y": rows["row_y"][idx], "p": pred, "sig": rows["row_sig"][idx]}); out = []
    for r, g in d.groupby("r"):
        s = g[g.sig]; m = dict(run=z["runs"][r], target=z["targets"][run_target[r].item()], tier_data="A" if run_tier[r] else "B", n=len(g), n_sig=len(s))
        m["r_all"] = np.corrcoef(g.y, g.p)[0, 1] if g.p.std() > 0 else np.nan; m["r_sig"] = np.corrcoef(s.y, s.p)[0, 1] if len(s) > 2 and s.p.std() > 0 else np.nan
        m["sign_sig"] = (np.sign(s.y) == np.sign(s.p)).mean() if len(s) else np.nan; m["mse_all"] = float(((g.y - g.p) ** 2).mean()); out.append(m)
    return pd.DataFrame(out)
best = (-1e9, None, None); I = torch.eye(D, device=dev, dtype=torch.float64); I[0, 0] = 0
for lam in (1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0):
    beta = torch.linalg.solve(XtX + lam * I, Xty); vi, vp = predict(beta, va); vm = run_metrics(vi, vp); vb = vm.groupby("target").agg(r_sig=("r_sig", "mean"), n_sig=("n_sig", "sum")).dropna(); score = np.average(vb.r_sig, weights=vb.n_sig)
    P(f"lambda {lam:g}: val weighted r_sig {score:.4f} sign_sig {np.nanmean(vm.groupby('target').sign_sig.mean()):.4f}")
    if score > best[0]: best = (score, lam, beta)
P(f"best lambda {best[1]:g} val r_sig {best[0]:.4f}"); ti, tp = predict(best[2], te); tm = run_metrics(ti, tp); tm["fold"] = a.fold; tm["tier_model"] = "A"; tm.to_csv(f"{O}/per_run.tsv", sep="\t", index=False)
np.savez_compressed(f"{O}/pred_heldout.npz", idx=ti, pred=tp); np.save(f"{O}/beta.npy", best[2].cpu().numpy()); by = tm.groupby("target").mean(numeric_only=True)
P("held-out (by RBP):", {c: round(float(np.nanmean(by[c])), 4) for c in ("r_all", "r_sig", "sign_sig")}); P("DONE", O)
