#!/usr/bin/env python3
"""Phase 1 (deterministic, no model): fetch NCBI PubMed data for the fig_8 scz lit-review genes.
Query '"GENE"[TIAB] AND "Schizophrenia"[TIAB]', 2000-present, top 5 by relevance; also capture
total hit count. One esearch + one batched efetch per gene => real PMIDs by construction (no
hallucination possible). Paced at ~2.5 req/s (no API key), retry/backoff, incremental checkpoint.
Adapted from fig_7/v8/ad_db_pipeline/ncbi_prefetch.py (SCZ instead of AD; DFC literature track)."""
import json, os, sys, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET

HERE  = os.path.dirname(os.path.abspath(__file__))
GENES = os.path.join(HERE, "literature_track.txt")
CACHE = os.path.join(HERE, "ncbi_cache.json")
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
COMMON = {"tool": "consensus_analysis_fig8", "email": "gugene.kang@ucsf.edu"}
DELAY = 0.40  # ~2.5 req/s, safely under the 3/s unkeyed limit

def get(url):
    for wait in (0, 10, 30):
        if wait: time.sleep(wait)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "consensus_analysis_fig8/1.0"})
            with urllib.request.urlopen(req, timeout=30) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            sys.stderr.write(f"  retry ({e})\n")
    return None

def esearch(gene):
    term = f'"{gene}"[TIAB] AND "Schizophrenia"[TIAB]'
    q = urllib.parse.urlencode({**COMMON, "db": "pubmed", "term": term, "retmax": 5,
                                "retmode": "json", "sort": "relevance",
                                "mindate": 2000, "maxdate": 2026, "datetype": "pdat"})
    txt = get(f"{BASE}/esearch.fcgi?{q}")
    if not txt: return None, []
    try:
        res = json.loads(txt)["esearchresult"]
        return int(res.get("count", 0)), res.get("idlist", [])
    except Exception:
        return None, []

def efetch(pmids):
    if not pmids: return []
    q = urllib.parse.urlencode({**COMMON, "db": "pubmed", "id": ",".join(pmids), "retmode": "xml"})
    txt = get(f"{BASE}/efetch.fcgi?{q}")
    if not txt: return []
    recs = []
    try:
        root = ET.fromstring(txt)
    except Exception:
        return []
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID", "")
        title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
        journal = art.findtext(".//Journal/Title", "")
        year = art.findtext(".//JournalIssue/PubDate/Year", "") or art.findtext(".//JournalIssue/PubDate/MedlineDate", "")
        abst = " ".join("".join(a.itertext()) for a in art.findall(".//Abstract/AbstractText"))
        recs.append({"pmid": pmid, "title": title.strip(), "journal": journal.strip(),
                     "year": year.strip(), "abstract": abst.strip()})
    return recs

def main():
    genes = [g.strip() for g in open(GENES) if g.strip()]
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    todo = [g for g in genes if g not in cache]
    print(f"total {len(genes)} | cached {len(cache)} | to fetch {len(todo)}", flush=True)
    for i, g in enumerate(todo, 1):
        count, pmids = esearch(g); time.sleep(DELAY)
        recs = efetch(pmids); time.sleep(DELAY)
        cache[g] = {"count": count, "records": recs}
        if i % 10 == 0 or i == len(todo):
            json.dump(cache, open(CACHE, "w"))
            print(f"[{i}/{len(todo)}] {g}: count={count} pmids={len(recs)}", flush=True)
    json.dump(cache, open(CACHE, "w"))
    n_hits = sum(1 for g in cache if cache[g]["records"])
    print(f"DONE. {len(cache)} genes cached; {n_hits} have >=1 abstract; "
          f"{len(cache)-n_hits} have 0 SCZ hits.", flush=True)

if __name__ == "__main__":
    main()
