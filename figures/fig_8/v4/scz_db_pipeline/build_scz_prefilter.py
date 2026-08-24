#!/usr/bin/env python3
"""
Build SCZ gene prefilter database.
Sources: OMIM (all schizophrenia phenotype MIM entries via elink),
         Open Targets Platform (MONDO:0005090, all associations),
         ClinVar (schizophrenia, ≤10 genes per variant to exclude large CNVs),
         SZDB (www.szdb.org — attempted, may be inaccessible)
Output: scz_gene_databases_prefilter.csv
"""

import os
import csv
import json
import time
import requests

OUT_FILE = os.path.join(os.environ.get("REPO_DIR", "/home/gugene/code/git/kang-oldham-2026"), "figures/fig_8/v4/scz_db_pipeline/scz_gene_databases_prefilter.csv")
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
MAX_GENES_PER_CLINVAR_VARIANT = 10  # exclude large structural variants


def ncbi_get(url, params, delay=0.4):
    for attempt in range(3):
        try:
            r = requests.get(url, params=params, timeout=30)
            r.raise_for_status()
            time.sleep(delay)
            return r
        except Exception as e:
            print(f"    attempt {attempt+1} failed: {e}")
            time.sleep(10 if attempt == 0 else 30)
    return None


# ---------------------------------------------------------------------------
# OMIM: all schizophrenia phenotype MIM entries → gene IDs via elink
# ---------------------------------------------------------------------------
def get_omim_genes():
    print("=== OMIM (schizophrenia phenotype entries) ===")
    genes = {}

    # Step 1: search OMIM db for all schizophrenia phenotype entries
    params = {
        "db": "omim",
        "term": "schizophrenia[Disease/Phenotype]",
        "retmax": 1000,
        "retmode": "json"
    }
    r = ncbi_get(f"{BASE_URL}/esearch.fcgi", params)
    if not r:
        print("  OMIM esearch failed — skipping")
        return genes

    mim_ids = r.json().get("esearchresult", {}).get("idlist", [])
    print(f"  Found {len(mim_ids)} schizophrenia MIM entries")

    # Step 2: elink from all MIM entries to gene db (batches of 100)
    gene_ids = set()
    for i in range(0, len(mim_ids), 100):
        batch = mim_ids[i:i+100]
        params = {
            "dbfrom": "omim", "db": "gene",
            "id": ",".join(batch),
            "cmd": "neighbor",
            "retmode": "json"
        }
        r = ncbi_get(f"{BASE_URL}/elink.fcgi", params)
        if r:
            for ls in r.json().get("linksets", []):
                for lsd in ls.get("linksetdbs", []):
                    gene_ids.update(lsd.get("links", []))
        print(f"  elink batch {i//100+1}: {len(gene_ids)} gene IDs cumulative")

    print(f"  Total gene IDs from elink: {len(gene_ids)}")
    if not gene_ids:
        return genes

    # Step 3: esummary to convert gene IDs → symbols (batches of 200)
    gene_ids = list(gene_ids)
    for i in range(0, len(gene_ids), 200):
        batch = gene_ids[i:i+200]
        params = {"db": "gene", "id": ",".join(map(str, batch)), "retmode": "json"}
        r = ncbi_get(f"{BASE_URL}/esummary.fcgi", params)
        if not r:
            continue
        result = r.json().get("result", {})
        for gid in batch:
            rec = result.get(str(gid), {})
            symbol = rec.get("name", "")
            name = rec.get("description", "")
            taxid = str(rec.get("organism", {}).get("taxid", ""))
            if symbol and symbol != "Error" and taxid == "9606":
                genes[symbol] = name

    print(f"  OMIM total: {len(genes)} human genes")
    return genes


# ---------------------------------------------------------------------------
# Open Targets Platform: MONDO:0005090 (schizophrenia), all associations
# ---------------------------------------------------------------------------
OT_QUERY = """
query SCZGenes($efoId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $efoId) {
    associatedTargets(
      page: { index: $index, size: $size }
    ) {
      count
      rows {
        target { approvedSymbol approvedName }
        score
        datasourceScores { id score }
      }
    }
  }
}
"""


def get_opentargets_genes():
    print("=== Open Targets Platform (MONDO:0005090) ===")
    genes = {}
    page_size = 500
    page_index = 0
    total = None

    while True:
        variables = {"efoId": "MONDO_0005090", "index": page_index, "size": page_size}
        data = None
        for attempt in range(3):
            try:
                r = requests.post(OT_API,
                                  json={"query": OT_QUERY, "variables": variables},
                                  timeout=60)
                r.raise_for_status()
                data = r.json()
                break
            except Exception as e:
                print(f"    OT attempt {attempt+1} failed: {e}")
                time.sleep(10)

        if not data:
            print("  OT pagination stopped — query failed")
            break

        at = (data.get("data") or {}).get("disease", {}).get("associatedTargets", {})
        if total is None:
            total = at.get("count", 0)
            print(f"  Total OT associations: {total}")

        rows = at.get("rows", [])
        if not rows:
            break

        for row in rows:
            target = row.get("target", {})
            symbol = target.get("approvedSymbol", "")
            name = target.get("approvedName", "")
            score = row.get("score", 0)
            sources = [ds["id"] for ds in row.get("datasourceScores", [])
                       if ds.get("score", 0) > 0]
            if symbol:
                genes[symbol] = {"name": name, "score": score,
                                 "sources": "|".join(sources)}

        page_index += 1
        fetched = page_index * page_size
        print(f"  Fetched {min(fetched, total)}/{total}", end="\r", flush=True)
        if fetched >= total:
            break
        time.sleep(0.4)

    print(f"\n  Open Targets total: {len(genes)} genes")
    return genes


# ---------------------------------------------------------------------------
# ClinVar: schizophrenia variants, excluding large structural variants
# ---------------------------------------------------------------------------
def get_clinvar_genes():
    print(f"=== ClinVar (schizophrenia; max {MAX_GENES_PER_CLINVAR_VARIANT} genes/variant) ===")
    genes = {}

    params = {
        "db": "clinvar",
        "term": "schizophrenia[Disease/Phenotype]",
        "retmax": 0,
        "retmode": "json",
        "usehistory": "y"
    }
    r = ncbi_get(f"{BASE_URL}/esearch.fcgi", params)
    if not r:
        print("  esearch failed — skipping ClinVar")
        return genes

    esearch = r.json().get("esearchresult", {})
    total = int(esearch.get("count", 0))
    webenv = esearch.get("webenv", "")
    qkey = esearch.get("querykey", "")
    print(f"  ClinVar variants: {total}")

    skipped_cnv = 0
    for start in range(0, min(total, 20000), 500):
        params = {
            "db": "clinvar",
            "query_key": qkey,
            "WebEnv": webenv,
            "retstart": start,
            "retmax": 500,
            "retmode": "json"
        }
        r = ncbi_get(f"{BASE_URL}/esummary.fcgi", params)
        if not r:
            continue
        result = r.json().get("result", {})
        for uid in result.get("uids", []):
            variant = result.get(str(uid), {})
            gene_list = variant.get("genes", [])
            if len(gene_list) > MAX_GENES_PER_CLINVAR_VARIANT:
                skipped_cnv += 1
                continue
            for g in gene_list:
                sym = g.get("symbol", "")
                if sym:
                    genes[sym] = g.get("geneid", "")
        print(f"  Processed {min(start+500, total)}/{total} variants, "
              f"{len(genes)} genes so far", end="\r", flush=True)

    print(f"\n  ClinVar total: {len(genes)} genes ({skipped_cnv} large-CNV variants skipped)")
    return genes


# ---------------------------------------------------------------------------
# SZDB: curated schizophrenia candidate genes (www.szdb.org)
# ---------------------------------------------------------------------------
def get_szdb_genes():
    print("=== SZDB (www.szdb.org) ===")
    genes = {}
    urls = [
        "http://www.szdb.org/download.html",
        "http://www.szdb.org/genelist.php",
        "http://www.szdb.org/",
    ]
    for url in urls:
        try:
            r = requests.get(url, timeout=20,
                             headers={"User-Agent": "Mozilla/5.0"})
            print(f"  {url} → {r.status_code} ({len(r.text)} bytes)")
            if r.status_code == 200 and len(r.text) > 100:
                break
        except Exception as e:
            print(f"  {url} → failed: {e}")

    if not genes:
        print("  SZDB: not accessible — skipping")
    return genes


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    omim  = get_omim_genes()
    ot    = get_opentargets_genes()
    clinv = get_clinvar_genes()
    szdb  = get_szdb_genes()

    all_genes = sorted(set(omim) | set(ot) | set(clinv) | set(szdb))
    print(f"\nTotal unique genes across all sources: {len(all_genes)}")

    fieldnames = ["Gene", "Gene_name", "In_OMIM", "In_OpenTargets_Platform",
                  "OpenTargets_score", "OpenTargets_sources",
                  "In_ClinVar", "In_SZDB", "Databases"]

    rows = []
    for gene in all_genes:
        in_omim  = "Y" if gene in omim  else "N"
        in_ot    = "Y" if gene in ot    else "N"
        in_clinv = "Y" if gene in clinv else "N"
        in_szdb  = "Y" if gene in szdb  else "N"

        gene_name = (omim.get(gene)
                     or ot.get(gene, {}).get("name", "")
                     or "")
        ot_score   = ot[gene]["score"]   if gene in ot else ""
        ot_sources = ot[gene]["sources"] if gene in ot else ""

        dbs = []
        if in_omim  == "Y": dbs.append("OMIM")
        if in_ot    == "Y": dbs.append("OpenTargets_Platform")
        if in_clinv == "Y": dbs.append("ClinVar")
        if in_szdb  == "Y": dbs.append("SZDB")

        rows.append({
            "Gene": gene,
            "Gene_name": gene_name,
            "In_OMIM": in_omim,
            "In_OpenTargets_Platform": in_ot,
            "OpenTargets_score": ot_score,
            "OpenTargets_sources": ot_sources,
            "In_ClinVar": in_clinv,
            "In_SZDB": in_szdb,
            "Databases": "|".join(dbs)
        })

    with open(OUT_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} genes to {OUT_FILE}")
    print(f"  OMIM:            {sum(1 for r in rows if r['In_OMIM']=='Y')}")
    print(f"  OpenTargets:     {sum(1 for r in rows if r['In_OpenTargets_Platform']=='Y')}")
    print(f"  ClinVar:         {sum(1 for r in rows if r['In_ClinVar']=='Y')}")
    print(f"  SZDB:            {sum(1 for r in rows if r['In_SZDB']=='Y')}")


if __name__ == "__main__":
    main()
