"""Build sample_table.tsv and runinfo.csv for a GEO series from GEO (sample characteristics, SRX relation) and SRA runinfo.
Usage: python make_cohort_tables.py <GSE> [--write]   (without --write only a summary is printed; nothing is created)
Output (with --write): 05_precancer_cohorts/<GSE>/{sample_table.tsv (gsm, title, source, characteristics..., srx, runs, n_runs), runinfo.csv}"""
import sys, re, io, os, time, urllib.request, pandas as pd
gse = sys.argv[1]; write = "--write" in sys.argv; PB = "/well/blagden/users/eic789/PerturbBridge"
def get(u):
    for t in range(5):
        try: return urllib.request.urlopen(u, timeout=120).read().decode("utf-8", "replace")
        except Exception as e: time.sleep(3 + 3 * t)
    raise RuntimeError(u)
txt = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=gsm&form=text&view=brief"); rows = []; cur = None
for line in txt.splitlines():
    if line.startswith("^SAMPLE"): cur = {"gsm": line.split(" = ")[1].strip()}; rows.append(cur)
    elif cur is not None and " = " in line:
        k, v = line.split(" = ", 1)
        if k == "!Sample_title": cur["title"] = v
        elif k == "!Sample_source_name_ch1": cur["source"] = v
        elif k == "!Sample_library_strategy": cur["strategy"] = v
        elif k == "!Sample_organism_ch1": cur["organism"] = v
        elif k == "!Sample_characteristics_ch1" and ": " in v: a, b = v.split(": ", 1); cur[re.sub(r"\W+", "_", a.strip().lower())] = b.strip()
        elif k == "!Sample_relation" and "SRX" in v: cur["srx"] = re.search(r"SRX\d+", v).group(0)
S = pd.DataFrame(rows); print(gse, "GEO samples", len(S), "| with SRX", int(S.get("srx", pd.Series(dtype=str)).notna().sum()))
ser = get(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}&targ=self&form=text&view=brief")
for l in ser.splitlines():
    if l.startswith(("!Series_title", "!Series_overall_design", "!Series_pubmed_id", "!Series_relation = BioProject", "!Series_relation = SubSeries")): print("  " + l[:400])
for c in S.columns:
    if c not in ("gsm", "title", "srx") and S[c].nunique() <= 12: print(f"  {c}: " + ", ".join(f"{k} {v}" for k, v in S[c].value_counts(dropna=False).head(12).items()))
    elif c not in ("gsm", "title", "srx"): print(f"  {c}: {S[c].nunique()} distinct, e.g. {list(S[c].dropna().unique()[:4])}")
print("  titles e.g.", list(S.title.head(4)))
if "srx" not in S or S.srx.notna().sum() == 0: print("  NO SRA RAW DATA LINKED"); sys.exit(0)
srx = S.srx.dropna().unique(); R = []
for i in range(0, len(srx), 150):
    term = "+OR+".join(srx[i:i + 150]); ids = re.findall(r"<Id>(\d+)</Id>", get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmax=5000&term={term}")); time.sleep(0.5)
    for j in range(0, len(ids), 200): R.append(pd.read_csv(io.StringIO(get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&rettype=runinfo&retmode=text&id={','.join(ids[j:j + 200])}")))); time.sleep(0.5)
R = pd.concat(R).drop_duplicates("Run"); R = R[R.Experiment.isin(srx)]
print(f"  SRA runs {len(R)} | layout {dict(R.LibraryLayout.value_counts())} | selection {dict(R.LibrarySelection.value_counts())} | platform {dict(R.Model.value_counts())} | avgLength median {R.avgLength.median():.0f} | spots median {R.spots.median() / 1e6:.1f} M (min {R.spots.min() / 1e6:.1f}) | total size {R.size_MB.sum() / 1e3:.0f} GB (sra)")
runs = R.groupby("Experiment").Run.apply(lambda x: ";".join(sorted(x))); S["runs"] = S.srx.map(runs); S["n_runs"] = S.runs.fillna("").apply(lambda x: len([y for y in x.split(";") if y]))
if write:
    d = f"{PB}/05_precancer_cohorts/{gse}"; os.makedirs(d, exist_ok=True); S[S.n_runs > 0].to_csv(f"{d}/sample_table.tsv", sep="\t", index=False); R.to_csv(f"{d}/runinfo.csv", index=False); print("  written", d, "samples", int((S.n_runs > 0).sum()))
