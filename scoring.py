#!/usr/bin/env python3
"""
scoring.py

Now also considers secondary_structure & sasa in the final penalty.
"""

def mark_excluded_positions(length, ranges):
    excluded = set()
    for (s, e) in ranges:
        for i in range(s, e+1):
            if 1 <= i <= length:
                excluded.add(i)
    return excluded

def assign_domain_penalties(length, domain_ranges, penalty=2.0):
    d = {i: 0.0 for i in range(1, length+1)}
    for (s, e) in domain_ranges:
        for i in range(s, e+1):
            if i in d:
                d[i] = penalty
    return d

def assign_ptm_penalties(length, ptm_positions, base_penalty=1.0, distance_decay=True):
    d = {i: 0.0 for i in range(1, length+1)}
    if not ptm_positions:
        return d

    if not distance_decay:
        for p in ptm_positions:
            if 1 <= p <= length:
                d[p] = base_penalty
    else:
        cutoff = 3
        for i in range(1, length+1):
            dist = min(abs(i - p) for p in ptm_positions)
            if dist <= cutoff:
                d[i] = base_penalty / (1 + dist)
    return d

def combine_all_scores(query_seq, combined_features, domain_penalty_dict, ptm_penalty_dict, exclude_positions):
    results = []
    for feat in combined_features:
        i = feat["position"]
        if i in exclude_positions:
            continue

        # Existing penalty contributions
        domain_p = domain_penalty_dict.get(i, 0.0)
        ptm_p = ptm_penalty_dict.get(i, 0.0)
        entropy = feat.get("entropy", 0.0) or 0.0
        disorder = feat.get("disorder", 0.0) or 0.0

        # SS & SASA
        ss = feat.get("secondary_structure", None)
        sasa = feat.get("sasa", None)

        # Example scoring:
        conservation_pen = max(0.0, 5.0 - entropy)  # penalize strong conservation
        disorder_pen = 1.0 - disorder  # prefer disordered => lower penalty if disorder is high

        # Let's penalize buried residues. We'll define "normalized_sasa" ~ 0..1
        # E.g. if we guess max SASA for a residue is around 250-300. 
        # This is simplistic. Adjust as you see fit.
        sasa_pen = 0.0
        if sasa is not None:
            max_possible_sasa = 250.0
            normalized_sasa = min(sasa / max_possible_sasa, 1.0)
            # If you want to penalize BURIED residues, you can do penalty = (1 - normalized_sasa).
            sasa_pen = 1.0 - normalized_sasa

        # Suppose you also want to penalize alpha-helices more than coil.
        # We'll do a small example:
        # H => +0.5 penalty, E => +0.3, everything else => 0
        ss_pen = 0.0
        if ss == "H":
            ss_pen = 0.5
        elif ss == "E":
            ss_pen = 0.3

        total_pen = domain_p + ptm_p + conservation_pen + disorder_pen + sasa_pen + ss_pen

        results.append({
            "position": i,
            "residue": feat["residue"],
            "total_penalty": total_pen,
            "domain_penalty": domain_p,
            "ptm_penalty": ptm_p,
            "entropy": entropy,
            "disorder": disorder,
            "secondary_structure": ss,
            "sasa": sasa
        })

    # Sort ascending => "least-worst"
    results.sort(key=lambda x: x["total_penalty"])
    return results
