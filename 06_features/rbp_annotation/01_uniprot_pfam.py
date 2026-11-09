#!/usr/bin/env python
"""Annotate every ENCODE KD/eCLIP target with UniProt accession, protein length, Pfam domains and sequence
(for ESM2 later). UniProt REST, human reviewed entries, gene-name exact match; retries for the flaky proxy."""
import pandas as pd, json, time, urllib.request, urllib.parse, sys, os
PB = "/well/blagden/users/eic789/PerturbBridge"; OUT = f"{PB}/06_features/rbp_annotation"
ex = pd.read_csv(f"{PB}/01_encode_metadata/experiments.tsv", sep="\t", dtype=str).fillna("")
targets = sorted({t for t in ex[ex.cell_line.isin(["K562", "HepG2"])].target if t})
print("targets:", len(targets), file=sys.stderr)
def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}), timeout=120) as r:
                return json.load(r)
        except Exception as e:
            time.sleep(3 * (i + 1))
    return None
rows, fasta = [], []
for i, t in enumerate(targets):
    q = urllib.parse.quote(f'gene_exact:{t} AND organism_id:9606 AND reviewed:true')
    d = get(f"https://rest.uniprot.org/uniprotkb/search?query={q}&fields=accession,gene_primary,protein_name,length,xref_pfam,sequence&size=5")
    hit = None
    if d and d.get("results"):
        prim = [r for r in d["results"] if (r.get("genes") or [{}])[0].get("geneName", {}).get("value") == t]
        hit = (prim or d["results"])[0]
    if not hit:
        rows.append({"target": t, "uniprot": "", "protein_name": "", "length": "", "pfam": "", "note": "no reviewed human hit"}); continue
    pf = ";".join(x["id"] + "|" + (x.get("properties") or [{}])[0].get("value", "") for x in hit.get("uniProtKBCrossReferences", []) if x["database"] == "Pfam")
    rows.append({"target": t, "uniprot": hit["primaryAccession"],
                 "protein_name": hit.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", ""),
                 "length": hit.get("sequence", {}).get("length", ""), "pfam": pf, "note": ""})
    fasta.append(f">{t}|{hit['primaryAccession']}\n{hit['sequence']['value']}\n")
    if i % 50 == 0: print(i, t, file=sys.stderr)
pd.DataFrame(rows).to_csv(f"{OUT}/rbp_uniprot_pfam.tsv", sep="\t", index=False)
open(f"{OUT}/rbp_sequences.fa", "w").write("".join(fasta))
print("done; no-hit:", sum(1 for r in rows if r["note"]), file=sys.stderr)
