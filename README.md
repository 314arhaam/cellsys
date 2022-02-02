# CellSys
An open-source tool for building initial structures for bio-membranes and drug-delivery systems. 🐍🧬  
📝 Paper: https://onlinelibrary.wiley.com/doi/10.1002/jcc.26793

## Introduction
This Python-based open-source package is designed to make the process of initial structure generation for drug-delivery systems easier.  
Focusing on:
1. Phospholipid bi- and mono- layers
2. Micelles
3. Liposomes

Built-in utility tools and file formats are currently based on [GROMACS](https://www.gromacs.org/) software.
## Example
```python
import cellsys
"""
Example:
Build an 8x8 mixed phospholipid bilayer, containing DPPC and DMPC lipids.
"""
# specify monomers
monomers = ["DPPC", "DPPS"]

# upper and lower composition
comp_upper = [32, 32]
comp_lower = [64, 0]

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
## Result
Rendered visualization of the membrane, made using [VMD](https://www.ks.uiuc.edu/Research/vmd/) software.
<p align="center">
  <img src="https://github.com/314arhaam/cellsys/blob/main/graphics/dppc-dpps.png" width="400" title="DPPC/DPPS membrane">
</p>
