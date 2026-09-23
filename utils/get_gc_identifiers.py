"""
Retrieve missing identifiers (CAS or InChIKey) plus SMILES and Molecular Formula
for the GC standards list via PubChem PUG-REST.

Logic per row:
  - If InChIKey is available -> use it to look up CAS, SMILES, Molecular Formula
  - Else if CAS is available -> use it to look up InChIKey, SMILES, Molecular Formula
  - Rows with neither (shouldn't occur) are skipped with a warning

Input:  Excel file with columns 'Compound name', 'INCHIKEY', 'CAS Number'
Output: CSV with SMILES, Molecular Formula, and which identifier/method succeeded

Usage:
    python get_gc_identifiers.py input.xlsx output.csv
"""

import sys
import time
import requests
import pandas as pd

PUG_INCHIKEY = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{}/property/{}/TXT"
PUG_NAME = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{}/property/{}/TXT"

REQUEST_DELAY = 0.25  # ~4 requests/sec
TIMEOUT = 10
MAX_RETRIES = 2


def query_single(url_template: str, identifier: str, prop: str):
    """Query PubChem PUG-REST for a single property. Returns the value or None."""
    if not identifier or str(identifier).strip() == "" or str(identifier).strip().upper() in ("N/A", "NONE"):
        return None
    url = url_template.format(requests.utils.quote(str(identifier).strip()), prop)
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


def get_all_properties(base_url_template: str, identifier: str):
    """Fetch SMILES and Molecular Formula, plus whichever cross-identifier is
    needed. Each property is a separate call - PubChem TXT format doesn't
    reliably support multiple properties in one request (confirmed bug from
    an earlier version of this script)."""
    smiles = query_single(base_url_template, identifier, "CanonicalSMILES")
    time.sleep(REQUEST_DELAY)
    formula = query_single(base_url_template, identifier, "MolecularFormula")
    time.sleep(REQUEST_DELAY)
    inchikey = query_single(base_url_template, identifier, "InChIKey")
    time.sleep(REQUEST_DELAY)
    return smiles, formula, inchikey


def main(input_path: str, output_path: str):
    df = pd.read_excel(input_path)
    df.columns = [c.strip() for c in df.columns]

    results = []
    total = len(df)
    for i, row in df.iterrows():
        name = row.get("Compound name")
        inchikey_in = row.get("INCHIKEY")
        cas_in = row.get("CAS Number")

        has_inchikey = pd.notna(inchikey_in) and str(inchikey_in).strip().upper() not in ("N/A", "NONE", "")
        has_cas = pd.notna(cas_in) and str(cas_in).strip().upper() not in ("N/A", "NONE", "")

        smiles, formula, cas_out, inchikey_out = None, None, None, None
        method = None

        if has_inchikey:
            # Look up via InChIKey: get SMILES, Formula, and confirm/cross-check InChIKey itself
            smiles, formula, _ = get_all_properties(PUG_INCHIKEY, inchikey_in)
            inchikey_out = inchikey_in
            # CAS isn't retrievable as a "property" via PUG-REST (it's a synonym, not a property field)
            # so CAS stays None here - flagged in output for manual lookup if needed
            if smiles:
                method = "InChIKey (input)"
        elif has_cas:
            # Look up via CAS: get SMILES, Formula, InChIKey
            smiles, formula, inchikey_out = get_all_properties(PUG_NAME, cas_in)
            cas_out = cas_in
            if smiles:
                method = "CAS (input)"
        else:
            method = "NO IDENTIFIER PROVIDED"

        if method not in ("InChIKey (input)", "CAS (input)"):
            if not smiles:
                method = "NOT FOUND"

        results.append({
            "Compound name": name,
            "CAS Number": cas_out if cas_out else cas_in,
            "INCHIKEY": inchikey_out if inchikey_out else inchikey_in,
            "SMILES": smiles,
            "Molecular Formula": formula,
            "Retrieved via": method,
        })

        print(f"[{i+1}/{total}] {name} -> {formula or 'NOT FOUND'} ({method})")

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_path, index=False)

    n_found = out_df["SMILES"].notna().sum()
    print(f"\nDone. {n_found}/{total} compounds resolved with SMILES/Formula.")
    print(f"Saved to {output_path}")

    not_found = out_df[out_df["Retrieved via"].isin(["NOT FOUND", "NO IDENTIFIER PROVIDED"])]
    if len(not_found) > 0:
        print(f"\n{len(not_found)} compounds not resolved:")
        for _, r in not_found.iterrows():
            print(f"  - {r['Compound name']} (CAS: {r['CAS Number']}, InChIKey: {r['INCHIKEY']})")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python get_gc_identifiers.py <input.xlsx> <output.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
