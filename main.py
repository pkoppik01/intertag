#!/usr/bin/env python3
"""
main.py
"""

import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import msa
import features
import annotations
import scoring

def get_query_sequence():
    print("You will need either your protein sequence or a UniProt ID, as well as the corresponding PDB file for that sequence. \nPlease use AlphaFold Server to generate the PDB if you do not already have one.\n")

    print("Choose input type:")
    print("1. Provide an amino acid sequence")
    print("2. Provide a UniProt ID")
    choice = input("Enter 1 or 2: ").strip()
    if choice == "1":
        seq_input = input("Enter your amino acid sequence: ").strip()
        query_record = SeqRecord(Seq(seq_input), id="Query", description="User provided sequence")
    elif choice == "2":
        uniprot_id = input("Enter the UniProt ID: ").strip()
        query_record = msa.fetch_sequence_from_uniprot(uniprot_id)
    else:
        print("Invalid choice. Exiting.")
        sys.exit(1)
    return query_record

def main():
    # 1. Get query sequence
    query_record = get_query_sequence()
    full_seq = str(query_record.seq)
    SeqIO.write(query_record, "query.fasta", "fasta")

    # 2. Perform MSA
    msa_result, hits = msa.perform_msa(query_record)
    msa_seqs = [str(seq) for seq in msa_result]
    if not msa_seqs:
        print("No MSA sequences found, aborting.")
        sys.exit(1)
    query_alignment = msa_seqs[0]

    # 3. Compute features: Shannon entropy + AIUPred disorder
    entropy_list = features.calculate_shannon_entropy(msa_seqs)
    disorder_dict = features.predict_disordered_binding_regions(full_seq)

    # 4. Ask about PDB for STRIDE + freesasa
    ss_sasa_list = None
    
    pdb_file = input("Enter the path to your PDB file: ").strip()
    stride_path = "stride/stride"  # or the path to your stride executable
    chain_id = "A"          # assume we want chain A
    ss_sasa_dict = features.get_secondary_structure_and_sasa(pdb_file, stride_executable=stride_path, chain_id=chain_id)
    ss_sasa_list = features.parse_ss_sasa_for_chain(ss_sasa_dict, chain_id=chain_id)

    # 5. Combine all features into a single list
    combined_feats = features.combine_features(query_alignment, entropy_list, disorder_dict, ss_sasa_list)

    # 6. Annotate with InterProScan
    print("\nSubmitting to InterProScan for domain/tm/signal. Please wait...")
    domain_ranges, tm_ranges, signal_pep_ranges = annotations.get_domain_tm_signal(full_seq)

    print("\nSubmitting to MusiteDeep for ptm. Please hold...")
    ptm_positions = annotations.get_ptm_positions(full_seq)
    print("Annotation complete.")
    print(f"Domains: {domain_ranges}")
    print(f"Transmembrane: {tm_ranges}")
    print(f"Signal peptides: {signal_pep_ranges}")
    print(f"PTM positions: {ptm_positions}")

    # 7. Build domain + PTM penalties
    length = len(full_seq)
    domain_penalties = scoring.assign_domain_penalties(length, domain_ranges, penalty=2.0)
    ptm_penalties = scoring.assign_ptm_penalties(length, ptm_positions, base_penalty=1.0, distance_decay=True)

    # 8. Exclude signal + TM regions
    excluded_signal = scoring.mark_excluded_positions(length, signal_pep_ranges)
    excluded_tm = scoring.mark_excluded_positions(length, tm_ranges)
    exclude_positions = excluded_signal.union(excluded_tm)

    # 9. Combine everything => final ranking
    final_scores = scoring.combine_all_scores(
        full_seq,
        combined_feats,
        domain_penalties,
        ptm_penalties,
        exclude_positions
    )

    # 10. Print top hits
    print("\n=== Top 10 'least-worst' insertion sites ===")
    print("Pos\tAA\tTotPen\tDomPen\tPTMPen\tEntropy\tDisorder\tSS\tSASA")
    for row in final_scores[:10]:
        print(f"{row['position']}\t{row['residue']}\t"
              f"{row['total_penalty']:.2f}\t"
              f"{row['domain_penalty']:.2f}\t{row['ptm_penalty']:.2f}\t"
              f"{row['entropy']:.2f}\t{row['disorder']:.2f}\t"
              f"{row.get('secondary_structure','NA')}\t"
              f"{row.get('sasa','NA')}")

    # 11. Write all positions to CSV
    csv_filename = "final_scores.csv"
    with open(csv_filename, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Position", "Residue", "TotalPenalty",
                         "DomainPenalty", "PTMPenalty",
                         "Entropy", "Disorder", "SecondaryStructure", "SASA"])
        for row in final_scores:
            writer.writerow([
                row["position"],
                row["residue"],
                f"{row['total_penalty']:.2f}",
                f"{row['domain_penalty']:.2f}",
                f"{row['ptm_penalty']:.2f}",
                f"{row['entropy']:.2f}",
                f"{row['disorder']:.2f}",
                row.get("secondary_structure","NA"),
                row.get("sasa","NA")
            ])

    print(f"\nAll scored positions written to {csv_filename}. \nPlease check for annotation accuracy!")

if __name__ == "__main__":
    main()
