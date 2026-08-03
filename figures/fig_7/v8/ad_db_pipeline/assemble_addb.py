#!/usr/bin/env python3
"""Phase 4b: assemble ad_db_summary_table_{dfc,mtg}.csv + ad_db.R in v8/ from the merged
Sonnet-5 classification + OpenTargets scores + prefilter DB. Minimal-file output.
Rules per prompt_ad_db_v3.txt: DFC = literature track only; MTG = literature + prefilter track."""
import os
import json, csv, os, re
SC="/tmp/claude-1001/-home-gugene/7c5dcd16-a0b7-4476-940f-e565da2c0449/scratchpad"
V8=os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_7/v8")

CATS={1:"APP processing / trafficking",2:"Amyloid-beta generation / aggregation / clearance",
3:"Tau phosphorylation / clearance",4:"Autophagy-lysosome pathway",5:"Ubiquitin-proteasome system",
6:"Mitochondrial dysfunction",7:"Endosomal-lysosomal trafficking",8:"ER stress / UPR / ERAD",
9:"Axonal transport",10:"Synaptic failure / loss",11:"Neuroinflammation",12:"Ca2+ dysregulation",
13:"AD genome-wide association (GWAS)"}
COLS=["Gene","Region","N_papers","AD_mechanism","References","References_verified",
      "AD_involvement_verified","Agent3_resolution","Pubmed_total_hits","OpenTargets_AD_score"]

def load_overlap(region):
    import csv as c
    p=f"{V8}/{region}_overlaps.csv"
    return set(r[1] for i,r in enumerate(c.reader(open(p))) if i>0)

def main():
    merged=json.load(open(f"{SC}/agent1_merged_final.json"))
    cmap=json.load(open(f"{SC}/count_map.json"))
    cache=json.load(open(f"{SC}/ncbi_cache.json"))
    ot=json.load(open(f"{SC}/ot_scores.json"))
    pcount=json.load(open(f"{SC}/prefilter_counts.json"))
    # full curated prefilter DB (each region intersects it with its own overlap -> symmetric tracks)
    prefilter_set=set(r["Gene"] for r in csv.DictReader(open(f"{V8}/ad_gene_databases_prefilter.csv")))
    dfc_ov=load_overlap("dfc"); mtg_ov=load_overlap("mtg")
    # curator-preserved genes: real AD genes whose literature uses a protein-family name
    # (calcineurin=PPP3CB, spectrin=SPTAN1, ERAD=EDEM3) unreachable by symbol/alias PubMed query.
    CURATOR=json.load(open(f"{SC}/curator_genes.json"))
    for g,d in CURATOR.items():
        merged[g]={"ad_linked":True,"category":d["_category"],"mechanism":d["AD_mechanism"],
                   "n_papers":d["_n"],"refs":d["References"].split("|"),"curator":True}
        ph=str(d["Pubmed_total_hits"]); cmap[g]=int(ph) if ph.lstrip("-").isdigit() else -1
        sc=d["OpenTargets_AD_score"]; ot[g]={"ensembl":None,"score":(float(sc) if sc not in ("NA","") else None)}
    # prefilter DB details
    pref={}
    for r in csv.DictReader(open(f"{V8}/ad_gene_databases_prefilter.csv")):
        pref[r["Gene"]]=r

    def ot_str(g):
        v=ot.get(g);
        return f"{v['score']:.4f}" if v and v.get("score") is not None else "NA"
    def npapers(g):
        r=merged[g]; n=r["n_papers"]; cnt=cmap.get(g,0) or 0
        if isinstance(n,str): return n
        return "5+" if (n>=5 and cnt>5) else str(n)

    def lit_row(g,region):
        r=merged[g]; cur=r.get("curator")
        res=("Curator-preserved: AD evidence documented under a protein-family name "
             "(e.g. calcineurin, spectrin) unreachable by symbol/alias PubMed query; PMIDs from prior manual curation") \
            if cur else "Confirmed (Sonnet-5 review; PMIDs pre-fetched from NCBI)"
        return {"Gene":g,"Region":region,"N_papers":npapers(g),"AD_mechanism":r["mechanism"],
            "References":"|".join(map(str,r["refs"])),
            "References_verified":"Y (curator)" if cur else "Y",
            "AD_involvement_verified":"Y (curator)" if cur else "Y",
            "Agent3_resolution":res,
            "Pubmed_total_hits":cmap.get(g,-1),"OpenTargets_AD_score":ot_str(g)}
    def pref_row(g,region):
        d=pref.get(g,{}); dbs=d.get("Databases","").replace("|","; "); name=d.get("Gene_name",g)
        srcs=d.get("OpenTargets_sources","")
        mech=f"{name}: identified as AD-associated gene in {dbs}; mechanism not curated by literature review pipeline"
        if srcs: mech+=f" (OpenTargets sources: {srcs})"
        return {"Gene":g,"Region":region,"N_papers":"NA","AD_mechanism":mech,
            "References":dbs,"References_verified":f"Prefilter DB ({dbs})",
            "AD_involvement_verified":"Y (prefilter DB)",
            "Agent3_resolution":f"Included via prefilter database entry ({dbs})",
            "Pubmed_total_hits":pcount.get(g, cmap.get(g,-1)),"OpenTargets_AD_score":ot_str(g)}

    included=set(g for g,r in merged.items() if r.get("ad_linked"))
    # BOTH regions apply BOTH tracks identically (only input genes differ):
    # curator-preserved > prefilter-track (curated DB) > literature-track include.
    def build(region, ov):
        rows=[]; tracks=[]
        for g in sorted(ov):
            if merged.get(g,{}).get("curator"):
                rows.append(lit_row(g,region)); tracks.append((g,"lit"))
            elif g in prefilter_set:
                rows.append(pref_row(g,region)); tracks.append((g,"pref"))
            elif g in included:
                rows.append(lit_row(g,region)); tracks.append((g,"lit"))
        return rows, tracks
    dfc_rows,dfc_tracks=build("DFC",dfc_ov)
    mtg_rows,mtg_final=build("MTG",mtg_ov)
    # write CSVs
    for region,rows in (("dfc",dfc_rows),("mtg",mtg_rows)):
        with open(f"{V8}/ad_db_summary_table_{region}.csv","w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=COLS); w.writeheader(); w.writerows(rows)
        print(f"{region}: {len(rows)} rows -> ad_db_summary_table_{region}.csv")

    # ---- ad_db.R ----
    WEB_PRIORITY=[("In_OMIM","omim.org"),("In_ClinVar","ncbi.nlm.nih.gov/clinvar"),
                  ("In_AlzForum","alzforum.org"),("In_OpenTargets_Platform","platform.opentargets.org")]
    def web_ref(g):
        d=pref.get(g,{})
        for col,url in WEB_PRIORITY:
            if str(d.get(col,"")).lower() in ("true","1","yes"): return url
        return "platform.opentargets.org"
    def comment(txt):
        txt=re.sub(r"\s+"," ",txt).strip()
        return txt[:90]+("..." if len(txt)>90 else "")
    def emit_list(name, genes_tracks):
        lines=[f"{name} <- list("]
        # group literature genes by category, then a prefilter section
        lit=[g for g,t in genes_tracks if t=="lit"]
        prefg=[g for g,t in genes_tracks if t=="pref"]
        from collections import defaultdict
        bycat=defaultdict(list)
        for g in lit: bycat[merged[g]["category"]].append(g)
        for cat in sorted(bycat):
            lines.append(f"\n  # {CATS.get(cat,'Other')}")
            for g in sorted(bycat[cat]):
                r=merged[g]; ref=r["refs"][0] if r["refs"] else 0
                tag="[curator] " if r.get("curator") else ""
                lines.append(f"  {g} = list(n = {npapers(g)}, ref = {ref}),  # {tag}{comment(r['mechanism'])}")
        if prefg:
            lines.append("\n  # Prefilter database entries")
            for g in sorted(prefg):
                d=pref.get(g,{}); dbs=d.get("Databases","").replace("|","; ")
                lines.append(f'  {g} = list(n = 1, ref = "{web_ref(g)}"),  # {d.get("Gene_name",g)}; AD-associated in {dbs}')
        # fix n for 5+ (needs quotes)
        out="\n".join(lines)
        out=re.sub(r'n = (\d+)\+', r'n = "\1+"', out)
        out=out+"\n)\n"
        # remove trailing comma on the last list element (R rejects it)
        out=re.sub(r"\),(\s*#[^\n]*)?\n\)", r")\1\n)", out)
        return out

    hdr=('# ad_db: Alzheimer\'s disease gene database for DFC and MTG brain regions\n'
         '# Rebuilt 2026-07-22 from dfc_overlaps.csv / mtg_overlaps.csv via pre-fetched-PMID Sonnet-5 literature review + prefilter DB\n'
         '# n = number of AD-linked papers (up to 5; "5+" if more); n = 1 for curated-database genes\n'
         '# ref = PMID of most relevant reference (literature genes); ref = "website" for database-sourced genes\n\n')
    with open(f"{V8}/ad_db.R","w") as f:
        f.write(hdr)
        f.write(emit_list("ad_db_dfc",dfc_tracks)); f.write("\n")
        f.write(emit_list("ad_db_mtg",mtg_final))
    print(f"ad_db.R: ad_db_dfc={len(dfc_tracks)}, ad_db_mtg={len(mtg_final)}")
    # quick consistency
    print("DFC includes not in overlap:", [g for g,_ in dfc_tracks if g not in dfc_ov])
    print("MTG rows not in overlap:", [g for g,_ in mtg_final if g not in mtg_ov])

if __name__=="__main__": main()
