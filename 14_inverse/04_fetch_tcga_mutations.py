"""Fetch somatic mutations of all ENCODE knockdown targets (plus RBM10, SRSF2, ZRSR2) in six TCGA PanCanAtlas studies from the
cBioPortal public API. Output 14_inverse/tcga_mutations.tsv: study, sampleId (TCGA-XX-XXXX-01), gene, proteinChange, mutationType.
Runs on the login node (needs internet)."""
import json, time, urllib.request, pandas as pd
PB = "/well/blagden/users/eic789/PerturbBridge"; API = "https://www.cbioportal.org/api"
STUDIES = ["uvm_tcga_pan_can_atlas_2018", "luad_tcga_pan_can_atlas_2018", "brca_tcga_pan_can_atlas_2018", "ov_tcga_pan_can_atlas_2018", "laml_tcga_pan_can_atlas_2018", "skcm_tcga_pan_can_atlas_2018"]
genes = sorted(set(pd.read_csv(f"{PB}/11_splits/run_tiers.tsv", sep="\t").target) | {"RBM10", "SRSF2", "ZRSR2"})
def post(path, body):
    for k in range(5):
        try:
            req = urllib.request.Request(f"{API}{path}", data=json.dumps(body).encode(), headers={"Content-Type": "application/json", "Accept": "application/json"}); return json.load(urllib.request.urlopen(req, timeout=300))
        except Exception as e: err = e; time.sleep(5 * (k + 1))
    raise RuntimeError(f"{path}: {err}")
gm = post("/genes/fetch?geneIdType=HUGO_GENE_SYMBOL", genes); ent = {g["hugoGeneSymbol"]: g["entrezGeneId"] for g in gm}; print("genes mapped", len(ent), "/", len(genes), "| missing:", sorted(set(genes) - set(ent))[:20])
rows = []
for st in STUDIES:
    m = post(f"/molecular-profiles/{st}_mutations/mutations/fetch?projection=SUMMARY", {"sampleListId": f"{st}_all", "entrezGeneIds": list(ent.values())})
    inv = {v: k for k, v in ent.items()}
    for x in m: rows.append(dict(study=st, sampleId=x["sampleId"], patientId=x["patientId"], gene=inv.get(x["entrezGeneId"]), proteinChange=x.get("proteinChange"), mutationType=x.get("mutationType"), variantType=x.get("variantType")))
    sl = post(f"/sample-lists/{st}_all", {}) if False else None
    print(st, "mutations", len(m))
D = pd.DataFrame(rows); D.to_csv(f"{PB}/14_inverse/tcga_mutations.tsv", sep="\t", index=False)
trunc = D[D.mutationType.isin(["Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site"])]
print("\n== truncating mutations per (study, gene), top 25 ==\n", trunc.groupby(["study", "gene"]).sampleId.nunique().sort_values(ascending=False).head(25).to_string())
print("\n== hotspot genes ==\n", D[D.gene.isin(["SF3B1", "U2AF1", "SRSF2"])].groupby(["study", "gene", "proteinChange"]).size().sort_values(ascending=False).head(20).to_string())
# sample lists (all profiled samples) so wild-type sets are defined
S = []
for st in STUDIES:
    sl = json.load(urllib.request.urlopen(f"{API}/sample-lists/{st}_sequenced?projection=DETAILED", timeout=300)); S += [dict(study=st, sampleId=s) for s in sl["sampleIds"]]
pd.DataFrame(S).to_csv(f"{PB}/14_inverse/tcga_sequenced_samples.tsv", sep="\t", index=False); print("sequenced samples", len(S))
