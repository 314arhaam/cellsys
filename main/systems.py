import cellsys

membrane = cellsys.bilayer(8, 9, "charmm36")

membrane.load_monomer("DMP") # add DMPC
membrane.load_monomer("DMS") # add DMPS

membrane.make(vector = 1.,
              z_dist = 2.5,
              comp_upper = [71,1],
              comp_lower = [36,36])

membrane.write_gro()
