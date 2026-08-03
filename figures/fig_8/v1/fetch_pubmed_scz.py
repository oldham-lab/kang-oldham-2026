#!/usr/bin/env python3
"""
Fetch PubMed abstracts for SCZ Agent 1 literature search.
Query: "GENE"[TIAB] AND "schizophrenia"[TIAB], 2000-present, top 5 by relevance.
Output: pubmed_scz_raw.json (gene -> list of {pmid, title, abstract, year})
"""

import os
import json
import time
import requests

GENES_FILE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/agent_artifacts/dfc_pipeline_track.txt")
OUT_FILE   = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/pubmed_scz_raw.json")
BASE_URL   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_gene(gene, retries=3):
    query = f'"{gene}"[TIAB] AND "schizophrenia"[TIAB]'
    # esearch
    for attempt in range(retries):
        try:
            r = requests.get(f"{BASE_URL}/esearch.fcgi", params={
                "db": "pubmed", "term": query,
                "retmax": 5, "retmode": "json",
                "sort": "relevance",
                "mindate": "2000", "maxdate": "2100", "datetype": "pdat"
            }, timeout=30)
            r.raise_for_status()
            break
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(10 if attempt == 0 else 30)
            else:
                return {"gene": gene, "pmids": [], "error": str(e)}

    esearch = r.json().get("esearchresult", {})
    pmids = esearch.get("idlist", [])
    total = int(esearch.get("count", 0))
    time.sleep(0.4)

    if not pmids:
        return {"gene": gene, "total": total, "papers": []}

    # efetch abstracts
    for attempt in range(retries):
        try:
            r2 = requests.get(f"{BASE_URL}/efetch.fcgi", params={
                "db": "pubmed", "id": ",".join(pmids),
                "rettype": "abstract", "retmode": "xml"
            }, timeout=30)
            r2.raise_for_status()
            break
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(10 if attempt == 0 else 30)
            else:
                return {"gene": gene, "total": total, "pmids": pmids, "error": str(e)}

    time.sleep(0.4)

    # Parse XML manually (avoid lxml dependency)
    import re
    xml = r2.text
    articles = []
    for art in re.findall(r'<PubmedArticle>(.*?)</PubmedArticle>', xml, re.DOTALL):
        pmid_m  = re.search(r'<PMID[^>]*>(\d+)</PMID>', art)
        title_m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', art, re.DOTALL)
        abs_m   = re.search(r'<AbstractText[^>]*>(.*?)</AbstractText>', art, re.DOTALL)
        year_m  = re.search(r'<PubDate>.*?<Year>(\d{4})</Year>', art, re.DOTALL)
        pmid    = pmid_m.group(1)  if pmid_m  else ""
        title   = re.sub(r'<[^>]+>', '', title_m.group(1)) if title_m else ""
        abstract= re.sub(r'<[^>]+>', '', abs_m.group(1))   if abs_m   else ""
        year    = year_m.group(1) if year_m else ""
        articles.append({"pmid": pmid, "title": title, "abstract": abstract, "year": year})

    return {"gene": gene, "total": total, "papers": articles}


def main():
    with open(GENES_FILE) as f:
        genes = [g.strip() for g in f if g.strip()]

    print(f"Fetching PubMed data for {len(genes)} genes...")

    # Load existing results if resuming
    try:
        with open(OUT_FILE) as f:
            results = json.load(f)
        done = {r["gene"] for r in results}
        print(f"Resuming: {len(done)} already fetched")
    except FileNotFoundError:
        results = []
        done = set()

    for i, gene in enumerate(genes):
        if gene in done:
            continue
        result = fetch_gene(gene)
        results.append(result)
        n_papers = result.get("total", 0)
        n_fetched = len(result.get("papers", []))
        print(f"  [{i+1}/{len(genes)}] {gene:20s} total={n_papers:4d}  fetched={n_fetched}")

        # Save after every gene
        with open(OUT_FILE, "w") as f:
            json.dump(results, f, indent=2)

    print(f"\nDone. Results written to {OUT_FILE}")
    # Summary
    has_results = sum(1 for r in results if r.get("total", 0) > 0)
    print(f"Genes with ≥1 PubMed hit: {has_results}/{len(genes)}")


if __name__ == "__main__":
    main()
