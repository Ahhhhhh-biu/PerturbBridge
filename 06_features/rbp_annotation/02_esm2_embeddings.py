"""ESM2 (facebook/esm2_t33_650M_UR50D) mean-pooled per-protein embeddings for all RBP sequences (CPU). Output: parquet with target, uniprot, 1280-d vector."""
import os, sys, torch, numpy as np, pandas as pd
from transformers import AutoTokenizer, EsmModel
PB = "/well/blagden/users/eic789/PerturbBridge"; D = f"{PB}/06_features/rbp_annotation"
os.environ.setdefault("HF_HOME", f"{PB}/envs/hf_cache")
seqs = []
for line in open(f"{D}/rbp_sequences.fa"):
    if line.startswith(">"): name, acc = line[1:].split()[0].split("|")[:2]; seqs.append([name, acc, ""])
    else: seqs[-1][2] += line.strip()
print("proteins:", len(seqs), flush=True)
tok = AutoTokenizer.from_pretrained("facebook/esm2_t33_650M_UR50D"); model = EsmModel.from_pretrained("facebook/esm2_t33_650M_UR50D").eval()
torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
rows = []
with torch.no_grad():
    for i, (name, acc, s) in enumerate(seqs):
        s = s[:2000]  # truncate very long proteins to bound memory
        x = tok(s, return_tensors="pt")
        h = model(**x).last_hidden_state[0, 1:-1]  # drop BOS/EOS
        rows.append((name, acc, len(s), h.mean(0).numpy().astype(np.float32)))
        if i % 25 == 0: print(i, name, len(s), flush=True)
df = pd.DataFrame({"target": [r[0] for r in rows], "uniprot": [r[1] for r in rows], "length_used": [r[2] for r in rows]})
E = np.stack([r[3] for r in rows]); np.save(f"{D}/rbp_esm2_650M_mean.npy", E); df.to_csv(f"{D}/rbp_esm2_650M_index.tsv", sep="\t", index=False)
print("saved", E.shape)
