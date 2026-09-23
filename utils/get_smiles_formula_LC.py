"""
Retrieve SMILES and Molecular Formula for chemicals via PubChem PUG-REST.
Tries InChIKey first (most reliable, unambiguous), falls back to CAS, then name.

Input:  Excel file with columns 'Compound name', 'CAS', 'INCHIKEY'
Output: CSV with SMILES, Molecular Formula, and which identifier/method succeeded

Usage:
    python get_smiles_formula.py input.xlsx output.csv
"""

import sys
import time
import requests
import pandas as pd

PUG_INCHIKEY_SMILES = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/CanonicalSMILES/TXT"
PUG_INCHIKEY_FORMULA = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/MolecularFormula/TXT"
PUG_NAME_SMILES = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/CanonicalSMILES/TXT"
PUG_NAME_FORMULA = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/MolecularFormula/TXT"

REQUEST_DELAY = 0.25  # ~4 requests/sec
TIMEOUT = 10
MAX_RETRIES = 2


def query_single(url_template: str, identifier: str):
    """Query PubChem PUG-REST for a single property. Returns the value or None."""
    if not identifier or str(identifier).strip() == "" or str(identifier).strip().lower() == "not provided":
        return None
    url = url_template.format(requests.utils.quote(str(identifier).strip()))
    for _ in range(MAX_RETRIES):
        try:
            resp = requests.get(url, timeout=TIMEOUT)
            if resp.status_code == 200:
                val = resp.text.strip()
                return val.splitlines()[0] if val else None
            elif resp.status_code == 404:
                return None
            else:
                time.sleep(1)
        except requests.exceptions.RequestException:
            time.sleep(1)
    return None


def query_pubchem(smiles_url: str, formula_url: str, identifier: str):
    """Fetch SMILES and Molecular Formula as two separate calls (TXT format doesn't
    reliably support multiple properties at once - this was the bug in the first version)."""
    smiles = query_single(smiles_url, identifier)
    time.sleep(REQUEST_DELAY)
    formula = query_single(formula_url, identifier)
    time.sleep(REQUEST_DELAY)
    return smiles, formula


def main(input_path: str, output_path: str):
    df = pd.read_excel(input_path)
    df.columns = [c.strip() for c in df.columns]

    results = []
    total = len(df)
    for i, row in df.iterrows():
        name = row.get("Compound name")
        cas = row.get("CAS")
        inchikey = row.get("INCHIKEY")

        smiles, formula = None, None
        method = None

        # Try InChIKey first - unambiguous, no name-matching risk
        if pd.notna(inchikey):
            smiles, formula = query_pubchem(PUG_INCHIKEY_SMILES, PUG_INCHIKEY_FORMULA, inchikey)
            if smiles:
                method = "InChIKey"

        # Fall back to CAS
        if not smiles and pd.notna(cas):
            smiles, formula = query_pubchem(PUG_NAME_SMILES, PUG_NAME_FORMULA, cas)
            if smiles:
                method = "CAS"

        # Fall back to name
        if not smiles and pd.notna(name):
            smiles, formula = query_pubchem(PUG_NAME_SMILES, PUG_NAME_FORMULA, name)
            if smiles:
                method = "Name"

        if not smiles:
            method = "NOT FOUND"

        results.append({
            "Compound name": name,
            "CAS": cas,
            "INCHIKEY": inchikey,
            "SMILES": smiles,
            "Molecular Formula": formula,
            "Retrieved via": method,
        })

        print(f"[{i+1}/{total}] {name} -> {formula or 'NOT FOUND'} ({method})")

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_path, index=False)

    n_found = (out_df["Retrieved via"] != "NOT FOUND").sum()
    print(f"\nDone. {n_found}/{total} retrieved.")
    print(f"Saved to {output_path}")

    not_found = out_df[out_df["Retrieved via"] == "NOT FOUND"]
    if len(not_found) > 0:
        print(f"\n{len(not_found)} compounds not found via InChIKey, CAS, or name:")
        for _, r in not_found.iterrows():
            print(f"  - {r['Compound name']} ({r['CAS']})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_smiles_formula.py <input.xlsx> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
