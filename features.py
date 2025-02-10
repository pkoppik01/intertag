#!/usr/bin/env python3
"""
features.py

Compute Shannon entropy, predict disorder (aiupred), and now parse
secondary structure (STRIDE) + SASA (freesasa) for the user-provided PDB.
"""

import math
from collections import Counter
import subprocess
import freesasa

def calculate_shannon_entropy(msa_sequences, ignore_gaps=True):
    if not msa_sequences:
        return []
    alignment_length = len(msa_sequences[0])
    entropies = []
    for pos in range(alignment_length):
        column = [seq[pos] for seq in msa_sequences]
        if ignore_gaps:
            column = [res for res in column if res != '-']
        if not column:
            entropies.append(0.0)
            continue
        counts = Counter(column)
        entropy = 0.0
        for aa, count in counts.items():
            p = count / len(column)
            entropy -= p * math.log2(p)
        entropies.append(entropy)
    return entropies

def predict_disordered_binding_regions(sequence):
    """
    Predict per-residue disorder using the aiupred library.
    Returns {pos: (aa, disorder_score, None)} with anchor_score=None.
    """
    from aiupred import aiupred_lib

    embedding_model, regression_model, device = aiupred_lib.init_models()
    disorder_scores = aiupred_lib.predict_disorder(sequence, embedding_model, regression_model, device)

    predictions = {}
    for i, score in enumerate(disorder_scores, start=1):
        aa = sequence[i-1]
        predictions[i] = (aa, float(score), None)
    return predictions

def get_secondary_structure_and_sasa(pdb_file, stride_executable="stride", chain_id="A"):
    """
    Uses STRIDE for secondary structure, freesasa for SASA.

    Returns a dict:  { (chain, residue_number) : (ss_code, sasa_value) }

    Example SS codes from STRIDE:
      H = alpha-helix, E = beta-strand, C = coil/other, ...
    """
    # --- 1. Parse secondary structure with STRIDE ---
    # Run stride on the PDB
    cmd = [stride_executable, pdb_file]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print("Error running STRIDE:", e)
        return {}

    ss_sasa = {}

    # STRIDE lines for residue annotations typically start with "ASG"
    # Format example:
    # ASG  MET A   1    H ...
    # ^^^  ^^^ ^   ^   ^ 
    #  1   2   3   4   5
    # We'll parse the residue number, chain, and SS code.
    # Example fix inside get_secondary_structure_and_sasa():
    for line in result.stdout.splitlines():
        if not line.startswith("ASG "):
            continue
        tokens = line.split()

        # tokens[3] = '1' (res num)
        # tokens[4] = '1' (some extra field, e.g. sequence index)
        # tokens[5] = 'C' (the single-letter SS code)
        # tokens[6] = 'Coil' (the descriptive name)

        if len(tokens) < 6:
            continue  # skip malformed lines
        
        chain = tokens[2]               # e.g. 'A'
        res_num = int(tokens[3])        # e.g. 1
        ss_code = tokens[5]            # e.g. 'C' or 'E' or 'H'
        
        # Store in your dict:
        ss_sasa[(chain, res_num)] = (ss_code, None)


    # --- 2. Parse SASA with freesasa ---
    structure = freesasa.Structure(pdb_file)
    result_freesasa = freesasa.calc(structure)
    residue_areas = result_freesasa.residueAreas()

    # residue_areas is a dict: chain -> {resnum: ResidueArea}, etc.
    for chain, residues_dict in residue_areas.items():
        for res_num, res_area_obj in residues_dict.items():
            sasa_val = res_area_obj.total

            res_num_int = int(res_num)
            # If (chain, res_num) already in ss_sasa => update. 
            # If not found, add it anyway.
            if (chain, res_num) in ss_sasa:
                old_ss, _ = ss_sasa[(chain, res_num)]
                ss_sasa[(chain, res_num)] = (old_ss, sasa_val)
            else:
                ss_sasa[(chain, res_num)] = (None, sasa_val)

    return ss_sasa

def parse_ss_sasa_for_chain(ss_sasa_dict, chain_id='A'):
    # First gather all entries (rnum, ss, sasa)
    entries = []
    for (c, rnum), (ss, sasa) in ss_sasa_dict.items():
        if c == chain_id:
            # Convert rnum to int if it might be string
            rnum_int = int(rnum)
            entries.append((rnum_int, ss, sasa))

    # Now merge any duplicates on rnum
    merged_map = {}  # rnum -> (ss, sasa)
    for (rnum_int, ss, sasa) in entries:
        if rnum_int not in merged_map:
            merged_map[rnum_int] = [None, None]  # [ss, sasa]
        # If we see a non-None ss, store it
        if ss is not None:
            merged_map[rnum_int][0] = ss
        # If we see a non-None sasa, store it
        if sasa is not None:
            merged_map[rnum_int][1] = sasa

    # Build final sorted list
    final_list = []
    for rnum_int in sorted(merged_map.keys()):
        ss_val, sasa_val = merged_map[rnum_int]
        final_list.append((ss_val, sasa_val))  # e.g. ('C', 250.45)

    return final_list

def combine_features(query_alignment, entropy_list, disorder_dict, ss_sasa_list=None):
    """
    Build a list of feature dicts for each residue in the aligned query sequence.
    We'll skip positions that are '-' in the alignment.
    """
    results = []
    query_pos = 0
    for i, char in enumerate(query_alignment):
        if char == '-':
            continue
        query_pos += 1

        # Basic features
        feat = {
            "position": query_pos,
            "residue": char,
            "entropy": entropy_list[i] if i < len(entropy_list) else 0.0,
        }

        # Disorder from aiupred
        if query_pos in disorder_dict:
            _, dscore, ascore = disorder_dict[query_pos]
            feat["disorder"] = dscore
            feat["anchor"] = ascore
        else:
            feat["disorder"] = 0.0
            feat["anchor"] = None

        # Secondary structure + SASA if provided
        if ss_sasa_list and (query_pos - 1) < len(ss_sasa_list):
            ss, sasa = ss_sasa_list[query_pos - 1]
            feat["secondary_structure"] = ss  # e.g. 'H', 'E', 'C', None
            feat["sasa"] = sasa              # numeric
        else:
            feat["secondary_structure"] = None
            feat["sasa"] = None

        results.append(feat)
    return results
