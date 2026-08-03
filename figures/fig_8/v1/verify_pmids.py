#!/usr/bin/env python3
"""Fetch PubMed metadata for all PMIDs in Agent 1 output for Agent 2 verification."""
import os
import csv, json, time, re, requests

IN_FILE  = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/scz_db_agent1_dfc.csv")
OUT_FILE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v1/pmid_verification.json")
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def fetch_pmids(pmids):
    for attempt in range(3):
        try:
            r = requests.get(f"{BASE_URL}/efetch.fcgi", params={
                "db":"pubmed","id":",".join(pmids),
                "rettype":"abstract","retmode":"xml"}, timeout=30)
            r.raise_for_status(); time.sleep(0.4); return r.text
        except Exception as e:
            print(f"  attempt {attempt+1} failed: {e}")
            time.sleep(10 if attempt==0 else 30)
    return None

def parse_articles(xml):
    out = {}
    for art in re.findall(r'<PubmedArticle>(.*?)</PubmedArticle>', xml, re.DOTALL):
        pmid_m  = re.search(r'<PMID[^>]*>(\d+)</PMID>', art)
        title_m = re.search(r'<ArticleTitle>(.*?)</ArticleTitle>', art, re.DOTALL)
        abs_m   = re.search(r'<AbstractText[^>]*>(.*?)</AbstractText>', art, re.DOTALL)
        year_m  = re.search(r'<PubDate>.*?<Year>(\d{4})</Year>', art, re.DOTALL)
        jour_m  = re.search(r'<Journal>.*?<ISOAbbreviation>(.*?)</ISOAbbreviation>', art, re.DOTALL)
        auth_m  = re.search(r'<AuthorList[^>]*>(.*?)</AuthorList>', art, re.DOTALL)
        first_author = ""
        if auth_m:
            ln = re.search(r'<LastName>(.*?)</LastName>', auth_m.group(1))
            first_author = ln.group(1) if ln else ""
        pmid  = pmid_m.group(1)  if pmid_m  else ""
        title = re.sub(r'<[^>]+>','',title_m.group(1)) if title_m else "NOT FOUND"
        abstr = re.sub(r'<[^>]+>','',abs_m.group(1))   if abs_m   else ""
        year  = year_m.group(1) if year_m else ""
        jour  = re.sub(r'<[^>]+>','',jour_m.group(1))  if jour_m  else ""
        out[pmid] = {"pmid":pmid,"title":title,"year":year,"journal":jour,
                     "first_author":first_author,"abstract":abstr[:500]}
    return out

# Collect all unique PMIDs
all_pmids = set()
with open(IN_FILE) as f:
    for row in csv.DictReader(f):
        for p in row["References"].split("|"):
            if p.strip(): all_pmids.add(p.strip())

print(f"Total unique PMIDs to verify: {len(all_pmids)}")
pmid_list = sorted(all_pmids)

results = {}
for i in range(0, len(pmid_list), 20):
    batch = pmid_list[i:i+20]
    xml = fetch_pmids(batch)
    if xml:
        parsed = parse_articles(xml)
        results.update(parsed)
        # Check which ones were not found
        for p in batch:
            if p not in parsed:
                results[p] = {"pmid":p,"title":"NOT FOUND","year":"","journal":"",
                              "first_author":"","abstract":""}
        print(f"  Verified {min(i+20,len(pmid_list))}/{len(pmid_list)} PMIDs")

with open(OUT_FILE,"w") as f:
    json.dump(results, f, indent=2)
print(f"Written to {OUT_FILE}")
print(f"Found: {sum(1 for v in results.values() if v['title'] != 'NOT FOUND')}/{len(all_pmids)}")
