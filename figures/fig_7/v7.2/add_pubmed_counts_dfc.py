import csv
import time
import urllib.parse
import requests

INPUT  = "ad_db_summary_table_dfc.csv"
OUTPUT = "ad_db_summary_table_dfc_pubmed_counts.csv"
NEW_COL = "Pubmed_total_hits"

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
DELAY = 0.4  # seconds between queries


def get_hit_count(gene: str) -> int:
    query = f'"{gene}"[TIAB] AND "Alzheimer"[TIAB]'
    params = {
        "db": "pubmed",
        "term": query,
        "mindate": "2000",
        "maxdate": "2026",
        "retmax": 0,
        "rettype": "json",
        "retmode": "json",
    }
    for attempt in range(3):
        try:
            resp = requests.get(BASE_URL, params=params, timeout=15)
            resp.raise_for_status()
            data = resp.json()
            return int(data["esearchresult"]["count"])
        except Exception as e:
            wait = 10 if attempt == 0 else 30
            print(f"  [retry {attempt+1}] {gene}: {e} — waiting {wait}s")
            time.sleep(wait)
    print(f"  [FAILED] {gene} — recording -1")
    return -1


def main():
    with open(INPUT, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    out_fieldnames = list(fieldnames) + [NEW_COL]

    genes = [r["Gene"] for r in rows]
    print(f"Querying PubMed for {len(genes)} genes...\n")

    counts = {}
    for i, gene in enumerate(genes, 1):
        count = get_hit_count(gene)
        counts[gene] = count
        print(f"  [{i:>2}/{len(genes)}] {gene}: {count}")
        time.sleep(DELAY)

    with open(OUTPUT, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=out_fieldnames)
        writer.writeheader()
        for row in rows:
            row[NEW_COL] = counts[row["Gene"]]
            writer.writerow(row)

    print(f"\nDone. Written to {OUTPUT}")
    failed = [g for g, c in counts.items() if c == -1]
    if failed:
        print(f"Failed genes: {failed}")


if __name__ == "__main__":
    main()
