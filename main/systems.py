import cellsys

membrane=cellsys.bilayer(8,9,"charmm36")
membrane.load_monomer("DMP")
membrane.load_monomer("DMS")
membrane.make(1.,2.5,[71,1],[36,36])
membrane.write_gro()
