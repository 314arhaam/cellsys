# CellSys
An open-source tool for building initial structures for bio-membranes and drug-delivery systems. 🐍🧬  
📝 Paper: https://onlinelibrary.wiley.com/doi/10.1002/jcc.26793  

**This repository is the initial release of the software**  

## Introduction
### Purpose
This Python-based open-source package is designed to make the process of initial structure generation for drug-delivery systems easier.  
Focusing on:
1. Phospholipid bi- and mono- layers
2. Micelles
3. Liposomes

Built-in utility tools and file formats are currently based on [GROMACS](https://www.gromacs.org/) software.

### Data
CellSys needs the following data to build initial structures:
1. `[ moleculetype ]` and `[ atom ]` names based on an `*.itp` file.
2. xyz-Coordinate of the residue, as a `numpy` array. (`*.npy` file)

To add a new forcefield, put above data in a folder and add to `cellsys_data` directory. If correct, it will be appeared after importing CellSys in `python`.

### How to use
1. Download `main` folder and copy in a directory.
2. 💻 **Terminal Mode**: Open `python` in the directory and `import cellsys`  
Or  
2. 📝 **Script Mode**: Run your desired commands as a `python` script.


## Example
```python
import cellsys
"""
Example:
Build an 8x8 mixed phospholipid bilayer, containing DPPC and DMPC lipids.
"""
# specify monomers
monomers = ["DPP", "DPS"]

# upper and lower composition
comp_upper = [32, 32]
comp_lower = [64, 0]

# make bilayer object, based on "CHARMM36" forcefield data
membrane = cellsys.bilayer(8, 8, "charmm36")

for monomer in monomers:
    membrane.load_monomer(monomer)

membrane.make(2., 4., comp_upper, comp_lower) # spaing, thickness, upper and lower composition
# write output in a *.gro file
membrane.write_gro()
# generate topology file (*.top)
cellsys.utills.gmx.make_topology(membrane)
```

Also, systems.py could be run inside Linux/Unix Terminal or Windows Command Prompt.
```shell
$ python systems.py
```
## Result
Rendered visualization of the membrane, made using [VMD](https://www.ks.uiuc.edu/Research/vmd/) software.
<p align="center">
  <img src="https://github.com/314arhaam/cellsys/blob/main/graphics/dppc-dpps.png" width="400" title="DPPC/DPPS membrane">
</p>

## Citation
```
@article{https://doi.org/10.1002/jcc.26793,
author = {Abbasi, Ali and Amjad-Iranagh, Sepideh and Dabir, Bahram},
title = {CellSys: An open-source tool for building initial structures for bio-membranes and drug-delivery systems},
journal = {Journal of Computational Chemistry},
volume = {43},
number = {5},
pages = {331-339},
keywords = {cell membranes, molecular dynamics simulations, phospholipid membranes, python programming language},
doi = {https://doi.org/10.1002/jcc.26793},
url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.26793},
eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/jcc.26793},
year = {2022}
}
```
