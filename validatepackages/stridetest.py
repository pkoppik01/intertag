import subprocess

# Run STRIDE (assuming 'stride' is in your PATH) on a PDB file:
result = subprocess.run(['stride/stride', 'examples/mng3kfp11pred1.pdb'], capture_output=True, text=True)
print(result.stdout)
