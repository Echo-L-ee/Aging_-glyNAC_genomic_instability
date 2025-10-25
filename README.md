# PubMed Search for Aging Hallmarks & Supplement Evidence


## Goal
Searches PubMed for papers supporting that GlyNAC rescues genomic instability, and saves a CSV (PMID, title, authors, journal, DOI, abstract, year).

## Install(VS Code + venv)
**Windows (PowerShell)**
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install biopython pandas python-dotenv


## Setup
Create `.env`:
EMAIL=XXX@hotmail.com
NCBI_API_KEY=XXXXXX

## Run
python pubmed_search.py

## Output
CSV with: PMID, Year, Title, Journal, Authors, DOI, Abstract.retmax=200

## Notes
- Default query targets GlyNAC + genomic/genome/chromosomal instability and DNA damage terms.
- Add an API key to avoid NCBI rate limits.
