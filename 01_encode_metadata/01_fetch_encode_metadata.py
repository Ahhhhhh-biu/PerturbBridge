#!/usr/bin/env python
"""Fetch ENCODE metadata for shRNA RNA-seq, CRISPR RNA-seq and eCLIP experiments.

Outputs (tab-separated):
  experiments.tsv  one row per experiment
  files.tsv        one row per released file of interest
"""
import json, sys, time, urllib.request, urllib.parse, csv, os

BASE = "https://www.encodeproject.org"
OUT = os.path.dirname(os.path.abspath(__file__))
ASSAYS = ["shRNA RNA-seq", "CRISPR RNA-seq", "eCLIP"]
EXP_FIELDS = ["accession", "assay_title", "target.label", "target.genes.symbol",
              "biosample_ontology.term_name", "possible_controls.accession",
              "control_type", "lab.title", "date_released", "status",
              "replicates.biological_replicate_number",
              "replicates.library.nucleic_acid_term_name",
              "replicates.library.strand_specificity",
              "files.accession", "files.output_type", "files.file_format",
              "files.file_format_type", "files.assembly", "files.genome_annotation",
              "files.biological_replicates", "files.file_size", "files.href",
              "files.md5sum", "files.status", "files.preferred_default"]
KEEP_OUTPUT = {"transcript quantifications", "gene quantifications", "alignments",
               "peaks", "replicated peaks"}

def get(url, tries=6):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.load(r)
        except Exception as e:
            print(f"retry {i+1}: {e}", file=sys.stderr); time.sleep(5 * (i + 1))
    raise RuntimeError(url)

exp_rows, file_rows = [], []
for assay in ASSAYS:
    q = {"type": "Experiment", "assay_title": assay, "status": "released", "limit": "all"}
    url = BASE + "/search/?" + urllib.parse.urlencode(q) + "".join(f"&field={f}" for f in EXP_FIELDS)
    d = get(url)
    print(assay, len(d["@graph"]), file=sys.stderr)
    for e in d["@graph"]:
        tgt = e.get("target") or {}
        reps = e.get("replicates") or []
        exp_rows.append({
            "accession": e["accession"], "assay": assay,
            "target": tgt.get("label", ""),
            "target_genes": ";".join(g.get("symbol", "") for g in tgt.get("genes", [])),
            "cell_line": e.get("biosample_ontology", {}).get("term_name", ""),
            "controls": ";".join(c["accession"] for c in e.get("possible_controls", [])),
            "control_type": e.get("control_type", ""),
            "lab": e.get("lab", {}).get("title", ""),
            "date_released": e.get("date_released", ""),
            "n_bio_reps": len({r.get("biological_replicate_number") for r in reps}),
            "nucleic_acid": ";".join(sorted({(r.get("library") or {}).get("nucleic_acid_term_name", "") for r in reps})),
            "strand_specificity": ";".join(sorted({(r.get("library") or {}).get("strand_specificity", "") for r in reps})),
        })
        for f in e.get("files", []):
            if f.get("status") != "released" or f.get("output_type") not in KEEP_OUTPUT:
                continue
            if f.get("assembly") not in (None, "GRCh38"):
                continue
            file_rows.append({
                "experiment": e["accession"], "assay": assay, "target": tgt.get("label", ""),
                "cell_line": e.get("biosample_ontology", {}).get("term_name", ""),
                "file": f["accession"], "output_type": f["output_type"],
                "file_format": f["file_format"], "file_format_type": f.get("file_format_type", ""),
                "assembly": f.get("assembly", ""), "genome_annotation": f.get("genome_annotation", ""),
                "bio_reps": ";".join(map(str, f.get("biological_replicates", []))),
                "file_size": f.get("file_size", ""), "md5sum": f.get("md5sum", ""),
                "preferred_default": f.get("preferred_default", ""),
                "url": BASE + f["href"],
            })

for name, rows in (("experiments.tsv", exp_rows), ("files.tsv", file_rows)):
    with open(os.path.join(OUT, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print(name, len(rows), file=sys.stderr)
