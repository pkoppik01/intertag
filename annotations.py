#!/usr/bin/env python3
"""
annotations.py

Queries InterProScan (EBI) to retrieve domain, TM, signal peptide,
and PTM-like regions from a protein sequence.

Requires 'requests' to be installed: pip install requests

Example usage in your main pipeline:
  domain_ranges, tm_ranges, signal_pep_ranges, ptm_positions = get_domain_tm_signal_ptm(sequence)
"""

import time
import requests
import json

IPRSCAN_RUN_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run"
IPRSCAN_STATUS_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status"
IPRSCAN_RESULT_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result"

def run_interproscan(sequence, email="someone@somewhere.com", app_name="MyApp"):
    """
    Submit sequence to InterProScan (EBI) via REST API.
    Returns a job_id if submission is successful, else None.
    
    Modification:
      - Added 'iprlookup': 'true' to the payload to instruct the service to perform
        the integrated lookup step. This mapping returns the corresponding InterPro entry
        (e.g., IPR009017) so that the 'entry' field is populated in the results.
    """
    payload = {
        'email': email,
        'title': 'TagFinderJob',
        'sequence': sequence,
        'app': app_name,
        'iprlookup': 'true'  # Added parameter to request integrated lookup
    }
    try:
        r = requests.post(IPRSCAN_RUN_URL, data=payload)
        r.raise_for_status()
        job_id = r.text.strip()
        return job_id
    except requests.exceptions.RequestException as e:
        print("Error submitting to InterProScan:", e)
        return None

def poll_job(job_id, max_wait=300):
    wait_time = 0
    poll_interval = 5
    while wait_time < max_wait:
        status_url = f"{IPRSCAN_STATUS_URL}/{job_id}"
        try:
            r = requests.get(status_url)
            r.raise_for_status()
            status = r.text  # e.g. "QUEUED", "RUNNING", "PENDING", "FINISHED"

            if status == "FINISHED":
                return True
            elif status in ["RUNNING", "PENDING", "QUEUED"]:
                time.sleep(poll_interval)
                wait_time += poll_interval
            else:
                print("Job finished with status:", status)
                return False

        except requests.exceptions.RequestException as e:
            print("Error polling InterProScan:", e)
            return False
    
    print(f"Job did not finish within {max_wait} seconds.")
    return False

def fetch_iprscan_result(job_id, fmt="json"):
    """
    Retrieve the InterProScan result in JSON format (by default).
    Possible 'fmt' values include: json, xml, tsv, etc.
    """
    result_url = f"{IPRSCAN_RESULT_URL}/{job_id}/{fmt}"
    try:
        r = requests.get(result_url)
        r.raise_for_status()
        return r.json()
    except requests.exceptions.RequestException as e:
        print("Error retrieving InterProScan result:", e)
        return None

def merge_overlaps(ranges):
    """
    Merge overlapping intervals in a list of (start, end) (1-based).
    Returns a sorted, merged list of (start, end).
    For example: [(1,5),(3,7)] -> [(1,7)]
    """
    if not ranges:
        return []
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged = [sorted_ranges[0]]
    for cur in sorted_ranges[1:]:
        last = merged[-1]
        # If current start is <= last.end + 1, overlap or adjacent
        if cur[0] <= last[1] + 1:
            merged[-1] = (last[0], max(last[1], cur[1]))
        else:
            merged.append(cur)
    return merged

def get_domain_tm_signal(sequence):
    """
    End-to-end function:
      1) Submits 'sequence' to InterProScan via run_interproscan()
      2) Polls until completion (poll_job())
      3) Fetches JSON (fetch_iprscan_result())
      4) Parses domain, TM, signal peptide, PTM info from the JSON

    Returns:
      domain_ranges          : list of (start, end) for domain regions
      tm_ranges              : list of (start, end) for transmembrane segments
      signal_peptide_ranges  : list of (start, end)
      ptm_positions          : list of single residue positions (e.g., phospho)
    """
    job_id = run_interproscan(sequence)
    if not job_id:
        print("Failed to submit job to InterProScan.")
        return [], [], [], []

    if not poll_job(job_id):
        print("InterProScan job did not complete successfully.")
        return [], [], [], []

    results_json = fetch_iprscan_result(job_id, fmt="json")
    if not results_json:
        print("No JSON results returned.")
        return [], [], [], []

    # Debug: see the actual structure
    # print("Raw InterProScan JSON:", results_json)

    # The top-level is a dict. We expect a key "results" which is a list.
    if "results" not in results_json:
        print("JSON does not contain 'results' at top level.")
        return [], [], [], []

    results_list = results_json["results"]
    if not results_list:
        print("No items in 'results' list.")
        return [], [], [], []

    # Take the first item in the results array
    top_result = results_list[0]
    if "matches" not in top_result:
        print("No 'matches' key in the first item of 'results'.")
        return [], [], [], []

    matches = top_result["matches"]

    domain_ranges = []
    tm_ranges = []
    signal_peptide_ranges = []
    ptm_positions = []

    # Parse each match
    for match in matches:
        # We'll try to detect domain vs. TM vs. signal vs. PTM 
        # by checking 'entry.type' and 'entry.description'.
        # Some older logic might check 'signature.type' or 'signatureDesc'.
        #print("DEBUG: match:", match)

        signature = match.get("signature", {})
        #print("DEBUG: signature:", signature)

        entry_data = signature.get("entry")
        #print("DEBUG: entry_data:", entry_data)  # Should now be populated if iprlookup worked
        
        signature_type = ""
        signature_desc = ""
        if entry_data:
            signature_type = entry_data.get("type", "")  # e.g., "FAMILY", "DOMAIN", etc.
            signature_desc = entry_data.get("description", "") or ""
            signature_name = entry_data.get("name", "") or ""
        else:
            # Fallback in case no integrated entry is provided
            signature_desc = signature.get("description", "") or ""
            signature_name = signature.get("name", "") or ""

        for loc in match.get("locations", []):
            start = loc["start"]
            end = loc["end"]

            # Simple logic to classify:
            if (any(kw in signature_desc.lower() for kw in ["transmembrane region", "transmembrane domain"]) or
                any(kw in signature_name.lower() for kw in ["transmembrane region", "transmembrane domain"]) or
                any(kw in signature_type.lower() for kw in ["transmembrane region", "transmembrane domain"])):
                tm_ranges.append((start, end))
            if (any(kw in signature_desc.lower() for kw in ["signal peptide", "signal region"]) or
                any(kw in signature_name.lower() for kw in ["signal peptide", "signal region"]) or
                any(kw in signature_type.lower() for kw in ["signal peptide", "signal region"])):
                signal_peptide_ranges.append((start, end))
            if signature_type.lower() in ["domain", "family", "homologous_superfamily"]:
                # We'll treat "homologous_superfamily" also as domain
                domain_ranges.append((start, end))
            # PTM check: naive approach looks for known keywords
            # You can refine this logic as needed
            if any(kw in signature_desc.lower() for kw in ["phospho", "glyco", "acetyl", "ubiquitin"]):
                for pos in range(start, end + 1):
                    ptm_positions.append(pos)

    # Remove duplicates & merge overlapping intervals
    print(domain_ranges)
    domain_ranges = merge_overlaps(domain_ranges)
    tm_ranges = merge_overlaps(tm_ranges)
    signal_peptide_ranges = merge_overlaps(signal_peptide_ranges)
    ptm_positions = sorted(set(ptm_positions))

    return domain_ranges, tm_ranges, signal_peptide_ranges

def get_ptm_positions(sequence):
    cutoff = 0.5
    modeloptions = ["Phosphoserine_Phosphothreonine",
                "Phosphotyrosine",
                "N-linked_glycosylation",
                "O-linked_glycosylation",
                "Ubiquitination",
                "SUMOylation",
                "N6-acetyllysine",
                "Methylarginine",
                "Methyllysine",
                "Pyrrolidone_carboxylic_acid",
                "S-palmitoyl_cysteine",
                "Hydroxyproline",
                "Hydroxylysine"]

    model=modeloptions[0]+";"+modeloptions[4] #for multiple models
    url = "http://api.musite.net/musitedeep/"+model+"/"+sequence
    myResponse = requests.get(url)
    if(myResponse.ok):
        # In this Example, jData are prediction results from MusiteDeep predictor
        data = json.loads(myResponse.content.decode('utf-8'))
        if "Error" in data.keys(): 
            print(data["Error in retrieving PTM data"])
    else:
        myResponse.raise_for_status()
    """
    Parses the JSON-like dictionary and returns a list of residue positions
    where at least one modification score exceeds the specified cutoff.
    """
    positions = []
    for result in data.get("Results", []):
        # Retrieve the modification score string
        mod_scores = result.get("Cutoff=0.5", "")
        # Split the string into individual modifications using semicolon as separator.
        # Each modification entry is expected to be in the format "ModificationName:score"
        for mod in mod_scores.split(';'):
            if not mod:
                continue  # Skip empty parts if any.
            try:
                # Split the modification entry into name and score.
                mod_name, score_str = mod.split(':')
                score = float(score_str)
                # Check if the score exceeds the cutoff.
                if score > cutoff:
                    positions.append(int(result.get("Residue")))
                    # Once we have a score above the cutoff for this residue,
                    # we can stop checking further modifications for the current result.
                    break
            except ValueError:
                print(f"Error processing modification: {mod}")

    ptm_positions = sorted(set(positions))

    return ptm_positions
