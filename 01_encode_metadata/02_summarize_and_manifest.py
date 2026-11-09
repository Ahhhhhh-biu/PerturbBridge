#!/usr/bin/env python
"""Summarise ENCODE metadata and write the download manifest (quantification TSVs + eCLIP peaks)."""
import pandas as pd, os
D = os.path.dirname(os.path.abspath(__file__))
ex = pd.read_csv(f"{D}/experiments.tsv", sep="\t", dtype=str).fillna("")
fi = pd.read_csv(f"{D}/files.tsv", sep="\t", dtype=str).fillna("")
fi["file_size"] = pd.to_numeric(fi["file_size"], errors="coerce")

kd = ex[ex.assay.isin(["shRNA RNA-seq", "CRISPR RNA-seq"])]
ctrl = kd[kd.control_type != ""]
kd_t = kd[kd.control_type == ""]
print("== experiments by assay / cell line (targeted only) ==")
print(kd_t.groupby(["assay", "cell_line"]).size().to_string())
print("\n== control experiments by assay / cell line ==")
print(ctrl.groupby(["assay", "cell_line", "control_type"]).size().to_string())
print("\n== unique RBP targets ==")
sh = set(kd_t[kd_t.assay == "shRNA RNA-seq"].target); cr = set(kd_t[kd_t.assay == "CRISPR RNA-seq"].target)
ec = set(ex[ex.assay == "eCLIP"].target)
print(f"shRNA targets {len(sh)} | CRISPR targets {len(cr)} | union KD {len(sh|cr)} | shRNA∩CRISPR {len(sh&cr)}")
print(f"eCLIP targets {len(ec)} | KD∩eCLIP {len((sh|cr)&ec)} | shRNA∩eCLIP {len(sh&ec)} | CRISPR∩eCLIP {len(cr&ec)}")
# same RBP, same cell line, both technologies
both = (kd_t[kd_t.assay=="shRNA RNA-seq"][["target","cell_line"]].drop_duplicates()
        .merge(kd_t[kd_t.assay=="CRISPR RNA-seq"][["target","cell_line"]].drop_duplicates()))
print(f"same RBP+cell line with BOTH shRNA and CRISPR (positive-control pairs): {len(both)}")
print("\n== eCLIP by cell line ==")
print(ex[ex.assay=="eCLIP"].groupby("cell_line").size().to_string())
print("\n== quantification files: annotation versions ==")
q = fi[fi.output_type.isin(["transcript quantifications","gene quantifications"])]
print(q.groupby(["assay","output_type","genome_annotation"]).size().to_string())
print("\n== missing controls ==")
print(f"targeted KD experiments without possible_controls: {(kd_t.controls=='').sum()}")

# manifest: quant TSVs for all KD + control experiments (GRCh38, prefer V29); eCLIP peaks (bed narrowPeak)
q = q[q.assembly=="GRCh38"]
ann = q.genome_annotation.value_counts().index[0]
q = q[q.genome_annotation==ann]
peaks = fi[(fi.assay=="eCLIP") & (fi.output_type.isin(["peaks","replicated peaks"])) & (fi.file_format=="bed") & (fi.assembly=="GRCh38")]
man = pd.concat([q, peaks])
man["dest"] = man.apply(lambda r: {"transcript quantifications":"transcript_quant","gene quantifications":"gene_quant"}.get(r.output_type,"eclip")
                        + "/" + r.file + ("." + r.file_format + ".gz" if r.file_format=="bed" else ".tsv"), axis=1)
man.to_csv(f"{D}/download_manifest.tsv", sep="\t", index=False)
print(f"\nmanifest rows {len(man)}; quant files {len(q)} ({ann}); eCLIP peak files {len(peaks)}; total {man.file_size.sum()/1e9:.1f} GB")
