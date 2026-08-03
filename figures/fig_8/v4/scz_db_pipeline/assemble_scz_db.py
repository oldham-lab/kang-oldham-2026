#!/usr/bin/env python3
"""Assemble scz_db_summary_table_dfc.csv + scz_db.R in v4/ from the merged Sonnet-5
classification + prefilter DB. DFC, BOTH tracks (literature + curated-DB prefilter).
Matches the v1 9-column schema so panel e / scz_db.R consumption is unchanged.
Adapted from fig_7/v8/ad_db_pipeline/assemble_addb.py."""
import os
import json, csv, os, re
HERE=os.path.dirname(os.path.abspath(__file__))
V4=os.path.dirname(HERE)
PREFILTER_CSV=os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/scz_gene_databases_prefilter.csv")

CATS={1:"Dopaminergic signaling dysfunction",2:"Glutamatergic / NMDA receptor hypofunction",
3:"GABAergic interneuron dysfunction",4:"Synaptic vesicle cycling / neurotransmitter release",
5:"Dendritic spine morphology or density",6:"Neurodevelopmental processes",
7:"Myelination / oligodendrocyte function",8:"Mitochondrial dysfunction",
9:"Neuroinflammation / complement / immune dysregulation",10:"Epigenetic regulation",
11:"RNA processing / splicing / post-transcriptional regulation",12:"Schizophrenia GWAS (common-variant)",
13:"Calcium-channel signaling",14:"Rare coding / CNV variant burden"}
COLS=["Gene","Region","N_papers","SCZ_mechanism","References","References_verified",
      "SCZ_involvement_verified","Agent3_resolution","Pubmed_total_hits"]
WEB_PRIORITY=[("In_OMIM","omim.org"),("In_ClinVar","ncbi.nlm.nih.gov/clinvar"),
              ("In_OpenTargets_Platform","platform.opentargets.org")]

def load(p): return json.load(open(p))

def main():
    merged=load(os.path.join(HERE,"scz_classification_final.json"))
    cache=load(os.path.join(HERE,"ncbi_cache.json"))
    alias=load(os.path.join(HERE,"alias_cache.json")) if os.path.exists(os.path.join(HERE,"alias_cache.json")) else {}
    pcount=load(os.path.join(HERE,"prefilter_counts.json"))
    curator=load(os.path.join(HERE,"curator_genes.json")) if os.path.exists(os.path.join(HERE,"curator_genes.json")) else {}
    overlap=[g.strip() for g in open(os.path.join(HERE,"agent_artifacts","dfc_candidates.txt")) if g.strip()]
    pref={r["Gene"]:r for r in csv.DictReader(open(PREFILTER_CSV))}
    prefilter_set=set(pref)

    # fold curator-preserved genes into merged
    for g,d in curator.items():
        merged[g]={"scz_linked":True,"category":d["_category"],"mechanism":d["SCZ_mechanism"],
                   "n_papers":d["_n"],"refs":[str(x) for x in d["_refs"]],"curator":True}

    def total_hits(g):
        c=cache.get(g,{}).get("count")
        a=alias.get(g,{}).get("alias_count")
        vals=[v for v in (c,a) if isinstance(v,int) and v>=0]
        return max(vals) if vals else -1
    def npapers(g):
        r=merged[g]; n=r["n_papers"]; cnt=total_hits(g)
        if isinstance(n,str): return n
        return "5+" if (n>=5 and cnt>5) else str(n)

    def lit_row(g):
        r=merged[g]; cur=r.get("curator")
        res=("Curator-preserved: SCZ evidence documented under a protein-family name unreachable "
             "by symbol/alias PubMed query; PMIDs from prior manual curation") if cur else \
            "Confirmed (Sonnet-5 review; PMIDs pre-fetched from NCBI)"
        return {"Gene":g,"Region":"DFC","N_papers":npapers(g),"SCZ_mechanism":r["mechanism"],
            "References":"|".join(map(str,r["refs"])),
            "References_verified":"Y (curator)" if cur else "Y",
            "SCZ_involvement_verified":"Y (curator)" if cur else "Y",
            "Agent3_resolution":res,"Pubmed_total_hits":total_hits(g)}
    def pref_row(g):
        d=pref.get(g,{}); dbs=d.get("Databases",""); name=d.get("Gene_name",g)
        srcs=d.get("OpenTargets_sources","")
        mech=f"{name}: identified as SCZ-associated gene in {dbs}; mechanism not curated by literature review pipeline"
        if srcs and "OpenTargets" in dbs: mech+=f" (OpenTargets sources: {srcs})"
        mech+="."
        return {"Gene":g,"Region":"DFC","N_papers":"NA","SCZ_mechanism":mech,
            "References":dbs,"References_verified":f"Prefilter DB ({dbs})",
            "SCZ_involvement_verified":"Y (prefilter DB)",
            "Agent3_resolution":f"Included via prefilter database entry ({dbs})",
            "Pubmed_total_hits":pcount.get(g,-1)}

    included=set(g for g,r in merged.items() if r.get("scz_linked"))
    # priority: curator (lit) > prefilter DB > literature include
    rows=[]; tracks=[]
    for g in sorted(overlap):
        if merged.get(g,{}).get("curator"):
            rows.append(lit_row(g)); tracks.append((g,"lit"))
        elif g in prefilter_set:
            rows.append(pref_row(g)); tracks.append((g,"pref"))
        elif g in included:
            rows.append(lit_row(g)); tracks.append((g,"lit"))
    with open(os.path.join(V4,"scz_db_summary_table_dfc.csv"),"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=COLS); w.writeheader(); w.writerows(rows)
    n_lit=sum(1 for _,t in tracks if t=="lit"); n_pref=sum(1 for _,t in tracks if t=="pref")
    print(f"scz_db_summary_table_dfc.csv: {len(rows)} rows ({n_lit} literature/curator + {n_pref} prefilter)")

    # ---- scz_db.R (match v1 format: n quoted, ref=PMID int for lit, ref=url for prefilter) ----
    def web_ref(g):
        d=pref.get(g,{})
        for col,url in WEB_PRIORITY:
            if str(d.get(col,"")).strip().upper()=="Y": return url
        return "platform.opentargets.org"
    def comment(txt):
        txt=re.sub(r"\s+"," ",txt).strip(); return txt[:88]+("..." if len(txt)>88 else "")
    lines=["scz_db_dfc <- list("]
    from collections import defaultdict
    bycat=defaultdict(list)
    for g,t in tracks:
        if t=="lit": bycat[merged[g]["category"]].append(g)
    for cat in sorted(bycat):
        lines.append(f"\n  # {CATS.get(cat,'Other')}")
        for g in sorted(bycat[cat]):
            r=merged[g]; ref=r["refs"][0] if r["refs"] else 0
            tag="[curator] " if r.get("curator") else ""
            lines.append(f"  {g} = list(n = '{npapers(g)}', ref = {ref}),  # {tag}{comment(r['mechanism'])}")
    prefg=[g for g,t in tracks if t=="pref"]
    if prefg:
        lines.append("\n  # Prefilter database entries")
        for g in sorted(prefg):
            d=pref.get(g,{}); dbs=d.get("Databases","")
            lines.append(f'  {g} = list(n = 1, ref = "{web_ref(g)}"),  # {d.get("Gene_name",g)}; SCZ-associated in {dbs}')
    out="\n".join(lines)+"\n)\n"
    out=re.sub(r"\),(\s*#[^\n]*)?\n\)", r")\1\n)", out)  # drop trailing comma on last elem
    hdr=("# scz_db: Schizophrenia gene database for DFC brain region\n"
         "# Rebuilt 2026-07-22 from dfc_overlaps.csv (HGNC-corrected) via pre-fetched-PMID Sonnet-5\n"
         "#   literature review (deterministic NCBI prefetch + alias augmentation) + curated prefilter DB.\n"
         "# n = number of SCZ-linked papers (up to 5; '5+' if more); n = 1 for curated-database genes\n"
         "# ref = PMID of most relevant reference (literature genes); ref = website for database-sourced genes\n\n")
    with open(os.path.join(V4,"scz_db.R"),"w") as f:
        f.write(hdr); f.write(out)
    print(f"scz_db.R: scz_db_dfc = {len(tracks)} entries ({n_lit} lit/curator + {n_pref} prefilter)")
    print("includes not in overlap:", [g for g,_ in tracks if g not in set(overlap)])

if __name__=="__main__": main()
