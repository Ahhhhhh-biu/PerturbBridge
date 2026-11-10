"""ESM2-650M mean embeddings for atlas proteins outside the ENCODE set (same model, pooling and 2,000-residue truncation as 06_features/rbp_annotation/02_esm2_embeddings.py).
Output: 06_features/rbp_functional/atlas/atlas_new_esm2.npy and atlas_new_esm2_index.tsv"""
import os, torch, numpy as np, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; A = f"{PB}/06_features/rbp_functional/atlas"; os.environ.setdefault("HF_HOME", f"{PB}/envs/hf_cache")
from transformers import AutoTokenizer, EsmModel; torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
seqs = []
for line in open(f"{A}/atlas_new_sequences.fa"):
    if line.startswith(">"): seqs.append(line[1:].strip().split("|") + [""])
    else: seqs[-1][2] += line.strip()
dev = "cuda" if torch.cuda.is_available() else "cpu"; tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D"); model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D").eval().to(dev); E = []
with torch.no_grad():
    for i, (g, acc, s) in enumerate(seqs):
        x = tok(s[:2000], return_tensors="pt").to(dev); E.append(model(**x).last_hidden_state[0, 1:-1].mean(0).float().cpu().numpy())
        if i % 100 == 0: print(i, g, len(s), flush=True)
np.save(f"{A}/atlas_new_esm2.npy", np.stack(E)); pd.DataFrame({"gene": [s[0] for s in seqs], "uniprot": [s[1] for s in seqs]}).to_csv(f"{A}/atlas_new_esm2_index.tsv", sep="\t", index=False); print("saved", len(E), dev)
