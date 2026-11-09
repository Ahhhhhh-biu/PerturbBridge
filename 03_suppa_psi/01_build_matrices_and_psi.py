#!/usr/bin/env python
"""Build transcript-TPM and gene-TPM matrices per cell line from ENCODE RSEM quantifications,
then run SUPPA2 psiPerEvent. Excludes 'whole organism' experiments. Sample id = <experiment>_rep<n>."""
import pandas as pd, numpy as np, os, subprocess, sys
PB = "/well/blagden/users/eic789/PerturbBridge"
PY = f"{PB}/envs/pb_venv/bin/python"; SUPPA = f"{PB}/envs/SUPPA/suppa.py"
man = pd.read_csv(f"{PB}/01_encode_metadata/download_manifest.tsv", sep="\t", dtype=str)
man = man[man.cell_line.isin(["K562", "HepG2"]) & man.assay.isin(["shRNA RNA-seq", "CRISPR RNA-seq"])]
OUT = f"{PB}/03_suppa_psi"; os.makedirs(f"{OUT}/matrices", exist_ok=True); os.makedirs(f"{OUT}/psi", exist_ok=True)

def load_quant(p, kind):
    """ENCODE ships two transcript-quant formats: RSEM (transcript_id, TPM) and kallisto
    (target_id = 'ENST|ENSG|...|', tpm). Gene quant is RSEM. Returns (Series TPM indexed by id, tool)."""
    hdr = open(p).readline().rstrip("\n").split("\t")
    if kind == "transcript" and hdr[0] == "target_id":
        d = pd.read_csv(p, sep="\t", usecols=["target_id", "tpm"])
        d["id"] = d.target_id.str.split("|").str[0]
        return d.set_index("id")["tpm"], "kallisto"
    idcol = "transcript_id" if kind == "transcript" else "gene_id"
    d = pd.read_csv(p, sep="\t", usecols=[idcol, "TPM"], index_col=0)
    d = d[d.index.astype(str).str.startswith("ENS")]          # drop spike-ins
    return d["TPM"], "rsem"

for cl in ["K562", "HepG2"]:
    for kind, idcol, valcol in [("transcript quantifications", "transcript_id", "TPM"), ("gene quantifications", "gene_id", "TPM")]:
        sub = man[(man.cell_line == cl) & (man.output_type == kind)]
        tag = "transcript" if kind.startswith("transcript") else "gene"
        if os.path.exists(f"{OUT}/matrices/{cl}_{tag}_TPM.tsv") and os.path.exists(f"{OUT}/matrices/{cl}_{tag}_samples.tsv"):
            print(cl, tag, "matrix exists, skipping rebuild", file=sys.stderr)
        else:
          # preallocate float32 matrix on a fixed feature index (union of ids seen in the first RSEM + kallisto file)
          paths = [(r, f"{PB}/02_encode_download/{r.dest}") for r in sub.itertuples()]
          paths = [(r, p) for r, p in paths if os.path.exists(p)]
          missing = len(sub) - len(paths)
          idx = set()
          for r, p in paths[:60]:                       # union of ids over the first 60 files covers both quant formats
              idx |= set(load_quant(p, tag)[0].index)
          idx = pd.Index(sorted(idx))
          X = np.full((len(idx), len(paths)), np.nan, dtype=np.float32)
          sids, sheet = [], []
          for j, (r, p) in enumerate(paths):
              v, tool = load_quant(p, tag)
              X[:, j] = v.reindex(idx).to_numpy(dtype=np.float32)
              sid = f"{r.experiment}_rep{r.bio_reps.replace(';', '_')}"
              if sid in sids: sid += f"_{r.file}"
              sids.append(sid)
              sheet.append({"sample": sid, "experiment": r.experiment, "assay": r.assay, "target": r.target,
                            "cell_line": cl, "bio_rep": r.bio_reps, "file": r.file, "quant_tool": tool})
          M = pd.DataFrame(X, index=idx, columns=sids)
          M = M[M.notna().any(axis=1)]
          M.to_csv(f"{OUT}/matrices/{cl}_{tag}_TPM.tsv", sep="\t", float_format="%.4f")
          pd.DataFrame(sheet).to_csv(f"{OUT}/matrices/{cl}_{tag}_samples.tsv", sep="\t", index=False)
          print(cl, tag, M.shape, "missing files:", missing, "tools:", pd.DataFrame(sheet).quant_tool.value_counts().to_dict(), file=sys.stderr)
          del M, X
    # SUPPA input: RSEM columns only, header WITHOUT the leading empty id cell (SUPPA counts header fields = n samples)
    sheet = pd.read_csv(f"{OUT}/matrices/{cl}_transcript_samples.tsv", sep="\t", dtype=str)
    rsem = sheet[sheet.quant_tool == "rsem"]["sample"].tolist()
    T = pd.read_csv(f"{OUT}/matrices/{cl}_transcript_TPM.tsv", sep="\t", index_col=0, usecols=["Unnamed: 0"] + rsem, dtype={c: np.float32 for c in rsem})
    T = T.fillna(0.0)
    sup = f"{OUT}/matrices/{cl}_transcript_TPM_rsem_suppa.tsv"
    with open(sup, "w") as fh:
        fh.write("\t".join(T.columns) + "\n")
        T.to_csv(fh, sep="\t", header=False, float_format="%.4f")
    print(cl, "SUPPA input", T.shape, file=sys.stderr)
    subprocess.run([PY, SUPPA, "psiPerEvent", "-i", f"{OUT}/events/gencode_v29_all_strict.ioe",
                    "-e", sup, "-o", f"{OUT}/psi/{cl}_events"], check=True)
    print(cl, "psi done", file=sys.stderr)
