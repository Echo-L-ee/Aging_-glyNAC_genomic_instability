
import os
import time
import argparse
import pandas as pd
from Bio import Entrez
from dotenv import load_dotenv
from pathlib import Path

# --- Setup ---
load_dotenv()  # loads .env if present
Entrez.email = os.getenv("EMAIL", "apriloctober17@hotmail.com")
Entrez.api_key = os.getenv("NCBI_API_KEY","48035e8fafb81cfbe03a235c56c14ecc3909")

# Simple rate limit to be polite to NCBI
SLEEP_SEC = 0.34  # ~3 calls/sec if you have API key; 0.5-1s if no key

def search_pubmed(query, retmax=1000, sort="relevance"):
    # 1) search PubMed by query → get PMIDs
    handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=retmax,
        sort=sort
    )
    results = Entrez.read(handle)
    handle.close()
    pmids = results.get("IdList", [])
    total = int(results.get("Count", 0))
    return pmids, total


def _species_block(species: str) -> str:
    # **NEW** helper to build animal species filter (always exclude humans)
    species = (species or "any").lower()
    if species == "mice":
        return "(Mice[MeSH Terms] NOT Humans[MeSH Terms])"  # **NEW**
    if species == "rats":
        return "(Rats[MeSH Terms] NOT Humans[MeSH Terms])"  # **NEW**
    return "(Animals[MeSH Terms] NOT Humans[MeSH Terms])"   # **NEW**


def _classify_tier(article) -> str:
    """
    **NEW** Classify record as 'clinical' or 'animal' (or '' if neither)
    Logic:
      - clinical: PublicationType includes a clinical trial type AND MeSH has Humans
      - animal: MeSH has Animals and NOT Humans
    """
    try:
        med = article.get("MedlineCitation", {})
        art = med.get("Article", {})
        # Publication types (strings)
        clinical_types = {
            "Clinical Trial",
            "Randomized Controlled Trial",
            "Controlled Clinical Trial",
            "Clinical Study",
            "Pragmatic Clinical Trial",
            "Multicenter Study",
        }  # **NEW** (broad but reasonable)
        ptypes = [str(pt) for pt in art.get("PublicationTypeList", [])]
        mesh = med.get("MeshHeadingList", [])
        mesh_terms = {str(mh.get("DescriptorName", "")) for mh in mesh}
        has_humans = "Humans" in mesh_terms
        has_animals = "Animals" in mesh_terms

        if any(pt in clinical_types for pt in ptypes) and has_humans:
            return "clinical"
        if has_animals and not has_humans:
            return "animal"
        return ""
    except Exception:
        return ""


def fetch_details(pmids):
    # 2) fetch details (title, authors, abstract, etc.) for those PMIDs
    if not pmids:
        return []
    joined = ",".join(pmids) # -> ex."39200123,38011234,37890001"	
    handle = Entrez.efetch(db="pubmed", id=joined, rettype="xml", retmode="xml") # This calls NCBI’s efetch to download the full article records for those PMIDs from the pubmed database, asking for XML
    records = Entrez.read(handle) 
    handle.close() #Biopython parses the XML into a nested Python object (dict/list-ish). Then closes the connection to free resources.
    articles = []
    for article in records["PubmedArticle"]: #Each PubMed article is inside the PubmedArticle list. Loop through them
        # Inside the loop, we defensively dig fields out of the nested structure:
        med = article.get("MedlineCitation", {})
        art = med.get("Article", {})
        journal = art.get("Journal", {})
        jl = journal.get("JournalIssue", {})
        pub_date = jl.get("PubDate", {} ) #The publication date lives inside JournalIssue → PubDate
        year = (
            pub_date.get("Year")
            or pub_date.get("MedlineDate")
            or ""
        ) #Try to get a clean publication year.  If Year is missing (older articles), MedlineDate might have something like “1998 Jan–Feb”. If both missing, fallback to empty string

        title = art.get("ArticleTitle", "")
        abstract = ""
        if "Abstract" in art and "AbstractText" in art["Abstract"]:
            # AbstractText could be list of sections
            parts = art["Abstract"]["AbstractText"]
            if isinstance(parts, list):
                abstract = " ".join(str(p) for p in parts)
            else:
                abstract = str(parts)

        authors = []
        for a in art.get("AuthorList", []):
            last = a.get("LastName") or ""
            fore = a.get("ForeName") or ""
            if last or fore:
                authors.append(f"{last} {fore}".strip())
        authors_str = "; ".join(authors)

        journal_title = journal.get("Title", "")
        pmid = med.get("PMID", "")
        # DOI lives in ArticleIdList sometimes
        doi = ""
        art_ids = article.get("PubmedData", {}).get("ArticleIdList", [])
        for aid in art_ids:
            if aid.attributes.get("IdType") == "doi":
                doi = str(aid)

        # **NEW** classify tier
        tier = _classify_tier(article)  # **NEW**

        articles.append({
            "PMID": str(pmid),
            "Year": str(year),
            "Title": str(title),
            "Journal": str(journal_title),
            "Authors": authors_str,
            "DOI": doi,
            "Abstract": abstract,
            "Tier": tier,  # **NEW**
        })
        time.sleep(SLEEP_SEC)
    return articles

def main():
# 3)Main program: parse options, build query, run search, save CSV
    parser = argparse.ArgumentParser(description="Search PubMed for GlyNAC → Genomic Instability")
    parser.add_argument("--retmax", type=int, default=1000, help="Max number of results")
    parser.add_argument("--sort", type=str, default="relevance", choices=["relevance", "pub+date", "most+recent"], help="Sort order",)
    parser.add_argument("--out", type=str, default=r"C:\Users\apartment\OneDrive\Echo\Aging Project\glyNAC_genomic_instability.csv", help="Output CSV filename")
    parser.add_argument("--query", type=str, default=None, help="Custom PubMed query")

    # **NEW** prioritization flags
    parser.add_argument("--filter_priority", choices=["clinical", "animal", "none"], default="none",
                        help="Prioritize clinical trials (humans) or animal studies (non-human).")  # **NEW**
    parser.add_argument("--species", choices=["any", "mice", "rats"], default="any",
                        help="Species focus when --filter_priority animal; 'any' = all non-human animals.")  # **NEW**

    args = parser.parse_args()

     # Default query targeting GlyNAC and genomic instability concepts
    default_query = (
    '("GlyNAC"[tiab] '
    'OR (glycine[tiab] AND ("N-acetylcysteine"[tiab] OR "N acetylcysteine"[tiab] OR acetylcystein*[tiab] OR NAC[tiab])) '
    'OR "Acetylcysteine"[MeSH]) '
    'AND ('
    '"genomic instabil*"[tiab] OR "chromosomal instabil*"[tiab] OR "DNA damage"[tiab] '
    'OR "DNA Damage"[MeSH] OR "Genome Stability"[MeSH] OR "Chromosome Instability"[MeSH] '
    'OR "double-strand break*"[tiab] OR "Double-Stranded DNA Breaks"[MeSH] '
    'OR micronucle*[tiab] OR "Micronuclei, Chromosome-Defective"[MeSH] '
    'OR aneuploid*[tiab] OR "Aneuploidy"[MeSH] '
    'OR "comet assay"[tiab] OR "Single-Cell Gel Comet Assay"[MeSH] '
    'OR γH2AX[tiab] OR H2AX[tiab] '
    'OR "oxidative stress"[tiab] OR "Oxidative Stress"[MeSH] '
    'OR ROS[tiab] OR "Reactive Oxygen Species"[MeSH] '
    'OR glutathione[tiab] OR GSH[tiab] OR "Glutathione"[MeSH] '
    'OR "mitochondrial dysfunction"[tiab] OR mitochondr*[tiab] OR "Mitochondria"[MeSH]'
    ')'
    )
    query = args.query or default_query

    # **NEW** apply priority filters to the query
    if args.filter_priority == "clinical":
        query += (
            " AND (Clinical Trial[ptyp] OR Randomized Controlled Trial[ptyp] "
            "OR Controlled Clinical Trial[ptyp] OR clinicaltrial[Filter]) "
            "AND Humans[MeSH Terms]"
        )  # **NEW**
    elif args.filter_priority == "animal":
        query += " AND " + _species_block(args.species)  # **NEW**

    print(f"Using query:\n{query}\n")

    # 1) Search PubMed: get list of PMIDs (limited by retmax) and total match count
    pmids, total = search_pubmed(query, retmax=args.retmax, sort=args.sort)
    print(f"PubMed total matches for this query: {total}")
    print(f"IDs returned (up to retmax): {len(pmids)}")

    # 2) Download article details for these PMIDs
    articles = fetch_details(pmids)

    # **NEW**: sort clinical first, then animal, then others; then by Year desc
    tier_order = {"clinical": 0, "animal": 1, "": 2}  # **NEW**
    df = pd.DataFrame(articles, columns=["PMID", "Year", "Title", "Journal", "Authors", "DOI", "Abstract", "Tier"])
    # convert year to numeric safely
    df["Year_num"] = pd.to_numeric(df["Year"].str[:4], errors="coerce")  # **NEW**
    df["Tier_rank"] = df["Tier"].map(tier_order).fillna(3).astype(int)   # **NEW**
    df = df.sort_values(["Tier_rank", "Year_num"], ascending=[True, False]).drop(columns=["Year_num", "Tier_rank"])  # **NEW**

    # 3) Save to CSV (Tier column at end)
    df.to_csv(args.out, index=False, encoding="utf-8")
    print(f"Saved {len(df)} records to {args.out}")


if __name__ == "__main__":
    main()






