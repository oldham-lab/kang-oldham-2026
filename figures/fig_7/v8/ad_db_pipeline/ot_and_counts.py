#!/usr/bin/env python3
"""Phase 4a (deterministic): OpenTargets AD-association scores for all final-table genes,
and PubMed total hit counts for the 96 MTG prefilter-track genes (lit genes already have
counts in ncbi_cache). Paced + retry. Outputs ot_scores.json, prefilter_counts.json."""
import json, os, time, urllib.request, urllib.parse
SC="/tmp/claude-1001/-home-gugene/7c5dcd16-a0b7-4476-940f-e565da2c0449/scratchpad"
OT="https://api.platform.opentargets.org/api/v4/graphql"
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMMON={"tool":"consensus_analysis_fig7","email":"gugene.kang@ucsf.edu"}

def post(url, payload):
    for wait in (0,10,30):
        if wait: time.sleep(wait)
        try:
            req=urllib.request.Request(url, data=json.dumps(payload).encode(),
                headers={"Content-Type":"application/json","User-Agent":"consensus_analysis_fig7/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r: return json.loads(r.read())
        except Exception as e: print("  OT retry",e)
    return None

def get(url):
    for wait in (0,10,30):
        if wait: time.sleep(wait)
        try:
            req=urllib.request.Request(url, headers={"User-Agent":"consensus_analysis_fig7/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r: return r.read().decode("utf-8","replace")
        except Exception as e: print("  NCBI retry",e)
    return None

def ot_score(gene):
    q={"query":'{search(queryString:"%s",entityNames:["target"]){hits{id name object{... on Target{approvedSymbol}}}}}'%gene}
    d=post(OT,q); time.sleep(0.35)
    if not d: return None
    ens=None
    for h in d.get("data",{}).get("search",{}).get("hits",[]):
        sym=(h.get("object") or {}).get("approvedSymbol","")
        if sym.upper()==gene.upper(): ens=h["id"]; break
    if not ens: return {"ensembl":None,"score":None}
    q2={"query":'{target(ensemblId:"%s"){associatedDiseases(BFilter:"MONDO_0004975",page:{size:1,index:0}){rows{score}}}}'%ens}
    d2=post(OT,q2); time.sleep(0.35)
    rows=(((d2 or {}).get("data",{}) or {}).get("target",{}) or {}).get("associatedDiseases",{}).get("rows",[])
    score=rows[0]["score"] if rows else None
    return {"ensembl":ens,"score":score}

def pubmed_count(gene):
    term=f'"{gene}"[TIAB] AND "Alzheimer"[TIAB]'
    q=urllib.parse.urlencode({**COMMON,"db":"pubmed","term":term,"retmax":0,"retmode":"json",
                              "mindate":2000,"maxdate":2026,"datetype":"pdat"})
    txt=get(f"{EUT}/esearch.fcgi?{q}"); time.sleep(0.4)
    try: return int(json.loads(txt)["esearchresult"]["count"])
    except Exception: return -1

def main():
    merged=json.load(open(f"{SC}/agent1_merged.json"))
    included=sorted(g for g,r in merged.items() if r.get("ad_linked"))
    prefilter=[g.strip() for g in open(f"{SC}/genes_mtg_prefilter.txt") if g.strip()]
    ot_genes=sorted(set(included)|set(prefilter))
    print(f"OT scores for {len(ot_genes)} genes; PubMed counts for {len(prefilter)} prefilter genes")
    ot=json.load(open(f"{SC}/ot_scores.json")) if os.path.exists(f"{SC}/ot_scores.json") else {}
    for i,g in enumerate(ot_genes,1):
        if g in ot: continue
        ot[g]=ot_score(g)
        if i%20==0 or i==len(ot_genes): json.dump(ot,open(f"{SC}/ot_scores.json","w")); print(f"  OT [{i}/{len(ot_genes)}]")
    json.dump(ot,open(f"{SC}/ot_scores.json","w"))
    pc=json.load(open(f"{SC}/prefilter_counts.json")) if os.path.exists(f"{SC}/prefilter_counts.json") else {}
    for i,g in enumerate(prefilter,1):
        if g in pc: continue
        pc[g]=pubmed_count(g)
        if i%20==0 or i==len(prefilter): json.dump(pc,open(f"{SC}/prefilter_counts.json","w")); print(f"  PMcount [{i}/{len(prefilter)}]")
    json.dump(pc,open(f"{SC}/prefilter_counts.json","w"))
    scored=sum(1 for g in ot if ot[g] and ot[g].get("score"))
    print(f"DONE. OT scored {scored}/{len(ot)}; prefilter counts {len(pc)}")

if __name__=="__main__": main()
