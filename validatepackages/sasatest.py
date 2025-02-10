import freesasa

# Load the PDB file and calculate SASA
structure = freesasa.Structure("examples/mng3kfp11pred1.pdb")
result = freesasa.calc(structure)

# Get per-residue SASA information (nested dict: chain -> residue number -> ResidueArea)
residue_areas = result.residueAreas()

# Iterate over each chain and residue to print the total SASA
for chain_id, residues in residue_areas.items():
    for res_num, res_area_obj in residues.items():
        # Access the 'total' attribute from the ResidueArea object
        total_area = res_area_obj.total
        print(f"Chain {chain_id} Residue {res_num} has total SASA: {total_area:.2f}")

# Get the total solvent accessible surface area
total_area = result.totalArea()
print("Total SASA:", total_area)