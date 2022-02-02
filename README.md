# CellSys
An open-source tool for building initial structures for bio-membranes and drug-delivery systems. 🐍🧬  
📝 Paper: https://onlinelibrary.wiley.com/doi/10.1002/jcc.26793

## Introduction
This Python-based open-source package is designed to make the process of initial structure generation for drug-delivery systems easier.  
Focusing on:
1. Phospholipid bi- and mono- layers
2. Micelles
3. Liposomes
## Example
```python
import cellsys
"""
Example:
Build an 8x8 mixed phospholipid bilayer, containing DPPC and DMPC lipids.
"""
# specify monomers
monomers = ["DPPC", "DMPC"]

# upper and lower composition
comp_upper = [40, 24]
comp_lower = [32, 32]

# make bilayer object, based on "CHARMM36" forcefield data
membrane = cellsys.bilayer(8, 8, "charmm36")

for monomer in monomers:
    membrane.load_monomer(monomer)

membrane.make(2, 4, comp_upper, comp_lower) # spaing, thickness, upper and lower composition
# write output in a *.gro file
membrane.write_gro()
# generate topology file (*.top)
cellsys.utills.gmx.make_topology(membrane)
```
