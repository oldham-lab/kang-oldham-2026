#!/usr/bin/env python3
"""Count-only esearch of SCZ PubMed total hits for the prefilter-track genes (populate
Pubmed_total_hits). '"GENE"[TIAB] AND "Schizophrenia"[TIAB]', 2000-present. ~2.5 req/s."""
import json, os, time, urllib.request, urllib.parse
HERE=os.path.dirname(os.path.abspath(__file__))
EUT="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMMON={"tool":"consensus_analysis_fig8","email":"gugene.kang@ucsf.edu"}
def get(url):
    for w in (0,10,30):
        if w: time.sleep(w)
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"consensus_analysis_fig8/1.0"})
            with urllib.request.urlopen(req,timeout=30) as r: return r.read().decode("utf-8","replace")
        except Exception as e: print("  retry",e)
    return None
def count(g):
    term=f'"{g}"[TIAB] AND "Schizophrenia"[TIAB]'
    q=urllib.parse.urlencode({**COMMON,"db":"pubmed","term":term,"retmax":0,"retmode":"json",
                              "mindate":2000,"maxdate":2026,"datetype":"pdat"})
    txt=get(f"{EUT}/esearch.fcgi?{q}")
    try: return int(json.loads(txt)["esearchresult"]["count"])
    except Exception: return -1
genes=[g.strip() for g in open(os.path.join(HERE,"agent_artifacts","prefilter_track.txt")) if g.strip()]
out={}
for g in genes:
    out[g]=count(g); time.sleep(0.4)
json.dump(out,open(os.path.join(HERE,"prefilter_counts.json"),"w"))
print("prefilter_counts.json:",len(out),"genes; nonzero:",sum(1 for v in out.values() if v>0))
