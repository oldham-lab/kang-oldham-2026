#!/usr/bin/env python3
"""Alias-augmented re-fetch (deterministic, NO model tokens). For each lit-track gene:
  1. NCBI Gene esummary -> OtherAliases (e.g. IARS1->IARS, MAPK8->JNK)
  2. PubMed esearch ("SYM"[TIAB] OR "ALIAS"[TIAB]...) AND "Schizophrenia"[TIAB], top5 relevance
  3. new_pmids = alias-query PMIDs NOT already in symbol-only cache
  4. efetch abstracts for new PMIDs
Output alias_cache.json + reclassify_set.json (not-yet-included genes that gained >=1 new abstract).
Adapted from fig_7/v8/ad_db_pipeline/alias_fetch.py (SCZ; DFC)."""
import json, os, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET
HERE=os.path.dirname(os.path.abspath(__file__))
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMMON={"tool":"consensus_analysis_fig8","email":"gugene.kang@ucsf.edu"}
DELAY=0.40
def get(url):
    for w in (0,10,30):
        if w: time.sleep(w)
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"consensus_analysis_fig8/1.0"})
            with urllib.request.urlopen(req,timeout=30) as r: return r.read().decode("utf-8","replace")
        except Exception as e: print("  retry",e)
    return None
def gene_aliases(sym):
    q=urllib.parse.urlencode({**COMMON,"db":"gene","term":f"{sym}[sym] AND Homo sapiens[orgn]","retmode":"json","retmax":1})
    txt=get(f"{EUT}/esearch.fcgi?{q}"); time.sleep(DELAY)
    try: uid=json.loads(txt)["esearchresult"]["idlist"][0]
    except Exception: return []
    q2=urllib.parse.urlencode({**COMMON,"db":"gene","id":uid,"retmode":"json"})
    txt=get(f"{EUT}/esummary.fcgi?{q2}"); time.sleep(DELAY)
    try:
        rec=json.loads(txt)["result"][uid]; al=[a.strip() for a in rec.get("otheraliases","").split(",") if a.strip()]
    except Exception: al=[]
    seen={sym.upper()}; out=[]
    for a in al:
        if len(a)>=3 and a.upper() not in seen: out.append(a); seen.add(a.upper())
        if len(out)>=6: break
    return out
def pubmed_alias(sym, aliases):
    terms=" OR ".join(f'"{t}"[TIAB]' for t in [sym]+aliases)
    term=f"({terms}) AND \"Schizophrenia\"[TIAB]"
    q=urllib.parse.urlencode({**COMMON,"db":"pubmed","term":term,"retmax":5,"retmode":"json",
                              "sort":"relevance","mindate":2000,"maxdate":2026,"datetype":"pdat"})
    txt=get(f"{EUT}/esearch.fcgi?{q}"); time.sleep(DELAY)
    try: res=json.loads(txt)["esearchresult"]; return int(res.get("count",0)), res.get("idlist",[])
    except Exception: return None,[]
def efetch(pmids):
    if not pmids: return []
    q=urllib.parse.urlencode({**COMMON,"db":"pubmed","id":",".join(pmids),"retmode":"xml"})
    txt=get(f"{EUT}/efetch.fcgi?{q}"); time.sleep(DELAY)
    if not txt: return []
    try: root=ET.fromstring(txt)
    except Exception: return []
    recs=[]
    for art in root.findall(".//PubmedArticle"):
        recs.append({"pmid":art.findtext(".//PMID",""),
            "title":("".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else "").strip(),
            "journal":art.findtext(".//Journal/Title","").strip(),
            "year":(art.findtext(".//JournalIssue/PubDate/Year","") or art.findtext(".//JournalIssue/PubDate/MedlineDate","")).strip(),
            "abstract":" ".join("".join(a.itertext()) for a in art.findall(".//Abstract/AbstractText")).strip()})
    return recs
def main():
    genes=[g.strip() for g in open(os.path.join(HERE,"literature_track.txt")) if g.strip()]
    orig=json.load(open(os.path.join(HERE,"ncbi_cache.json")))
    merged=json.load(open(os.path.join(HERE,"scz_classification.json")))
    included=set(g for g,r in merged.items() if r.get("scz_linked"))
    ac=json.load(open(os.path.join(HERE,"alias_cache.json"))) if os.path.exists(os.path.join(HERE,"alias_cache.json")) else {}
    todo=[g for g in genes if g not in ac and g not in included]  # skip already-included genes
    print(f"total {len(genes)} | included {len(included)} | cached {len(ac)} | to do {len(todo)}",flush=True)
    for i,g in enumerate(todo,1):
        aliases=gene_aliases(g)
        count,pmids=pubmed_alias(g,aliases)
        orig_pmids={r["pmid"] for r in orig.get(g,{}).get("records",[])}
        new_pmids=[p for p in pmids if p not in orig_pmids]
        new_recs=efetch(new_pmids) if new_pmids else []
        ac[g]={"aliases":aliases,"alias_count":count,"alias_pmids":pmids,"new_pmids":new_pmids,"new_records":new_recs}
        if i%10==0 or i==len(todo):
            json.dump(ac,open(os.path.join(HERE,"alias_cache.json"),"w"))
            print(f"[{i}/{len(todo)}] {g} aliases={aliases} new={len(new_recs)}",flush=True)
    json.dump(ac,open(os.path.join(HERE,"alias_cache.json"),"w"))
    reclass=sorted(g for g in ac if g not in included and any(r["abstract"] for r in ac[g]["new_records"]))
    json.dump(reclass,open(os.path.join(HERE,"reclassify_set.json"),"w"))
    tot_new=sum(1 for g in ac if any(r["abstract"] for r in ac[g]["new_records"]))
    print(f"DONE. genes with >=1 new abstract: {tot_new} | RECLASSIFY: {len(reclass)} -> {reclass}",flush=True)
if __name__=="__main__": main()
