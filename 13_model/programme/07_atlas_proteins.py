"""Protein list for the atlas: reviewed human UniProt entries annotated with GO:0003723 (RNA binding). Writes
06_features/rbp_functional/atlas/atlas_proteins.tsv (uniprot, gene, length, string_id, in_encode, target) and atlas_new_sequences.fa (proteins not in the ENCODE set)."""
import urllib.request, io, pandas as pd, time
PB = "/well/blagden/users/eic789/PerturbBridge"; A = f"{PB}/06_features/rbp_functional/atlas"
u = "https://rest.uniprot.org/uniprotkb/stream?query=(reviewed:true)%20AND%20(organism_id:9606)%20AND%20(go:0003723)&format=tsv&fields=accession,gene_primary,protein_name,length,xref_string,sequence"
for t in range(5):
    try: txt = urllib.request.urlopen(u, timeout=300).read().decode(); break
    except Exception as e: print("retry", e); time.sleep(10)
d = pd.read_csv(io.StringIO(txt), sep="\t"); d.columns = ["uniprot", "gene", "protein_name", "length", "string_id", "sequence"]; d["string_id"] = d.string_id.fillna("").str.split(";").str[0]
ids = pd.read_csv(f"{PB}/06_features/rbp_functional/rbp_id_map.tsv", sep="\t"); u2t = dict(zip(ids.uniprot, ids.target)); d["target"] = d.uniprot.map(u2t); d["in_encode"] = d.target.notna()
enc = ids[~ids.uniprot.isin(d.uniprot)]; print("GO:0003723 reviewed human", len(d), "| of these in ENCODE set", int(d.in_encode.sum()), "| ENCODE targets without the GO term (added with their existing descriptors)", len(enc))
d[d.gene.notna()].drop(columns="sequence").to_csv(f"{A}/atlas_proteins.tsv", sep="\t", index=False)
with open(f"{A}/atlas_new_sequences.fa", "w") as f:
    for r in d[~d.in_encode & d.gene.notna()].itertuples(): f.write(f">{r.gene}|{r.uniprot}\n{r.sequence}\n")
print("new proteins to embed", int((~d.in_encode & d.gene.notna()).sum()), "| LARP1 present", bool((d.gene == "LARP1").any()))
