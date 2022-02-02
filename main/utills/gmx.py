import numpy as np
import pandas as pd
import os

def gro_to_numpy(filename,suffix=''):
    gro_file=open(filename,"r")
    gro_text=gro_file.read()
    gro_file.close()
    data=dict()
    data["residue_number"]=[]
    data["residue_name"]=[]
    data["atom_name"]=[]
    data["atom_number"]=[]
    data["x"]=[]
    data["y"]=[]
    data["z"]=[]
    for line in gro_text.split("\n")[2:-2]:
        data["residue_number"].append(line[:5])
        data["residue_name"].append(line[5:10])
        data["atom_name"].append(line[10:15])
        data["atom_number"].append(line[15:20])
        data["x"].append(float(line[20:28]))
        data["y"].append(float(line[28:36]))
        data["z"].append(float(line[36:44]))
    data=pd.DataFrame(data)
    coords=np.array(data[["x","y","z"]])
    np.save(filename.split(".")[0],coords)
    return coords

def editconf(mol):
    """
    Put system into cg of its box
    Args:
        mol: system
    """
    box="  %8.3f %8.3f %8.3f"%tuple(mol.box)
    address="CellSys-%s/%s"%(mol.filename[:-4],mol.filename)
    os.system('gmx editconf -f %s -o %s -c -box %s'%(address,address,box))

def solvate(mol):
    pass

def compress_membrane(mol,tarea=80.,th=3.7):
    editconf(mol)
    address="CellSys-%s"%(mol.filename[:-4])
    os.system('cp -r utills/em.mdp %s'%(address))
    address=address+"/"+mol.filename
    h=[mol.thickness()]
    area=[mol.apl(h[-1])]
    while area[-1]>tarea or h[-1]>th:
        area_not_ok=int(area[-1]>tarea)
        h_not_ok=int(h[-1]>th)
        vec=np.array([[0.1*area_not_ok,0.1*area_not_ok,0.1*h_not_ok]])
        mol.compress(vector=vec)
        mol.write_gro()
        editconf(mol)
        grompp_success=not os.system("gmx grompp -f CellSys-%s/em.mdp -c %s -o %s[EMC].tpr -p CellSys-%s/topology.top -maxwarn 2"%(mol.filename[:-4],address,address[:-4],mol.filename[:-4]))
        mdrun_success=not os.system("gmx mdrun -deffnm %s[EMC] -v 1"%(address[:-4]))
        if not mdrun_success and not grompp_success:
            print("Error caused while running.")
            break
        area.append(mol.apl(h[-1]))
        h.append(mol.thickness())
        coords=gro_to_numpy("%s[EMC]"%(address[:-4]))
        mol.set_coords(coords)
        print("******AREA %6.3f %6.3f*******"%(area[-1],h[-1]))
    editconf(mol)
    return area

def make_topology(mol):
    """
    Make topology files. *.top

    Args:
        mol: An object of class system.
    """
    text='\n'
    # calculate the number of monomers. NOTE: it's just for bilayer
    for i,resname in enumerate(mol.ffresname):
        num=0
        for _,stack_name in enumerate(mol.coords):
            num+=mol.composition[stack_name][i]
        text+="%s %d\n"%(resname,num)
    # copy forcefield files
    address="CellSys-%s"%(mol.filename[:-4])
    os.system("cp -r data/forcefield/gromos43a1_lipid.ff %s"%(address))
    body_file=open('data/forcefield/gromos43a1_lipid.top','r')
    # generate text
    body=body_file.read()
    body_file.close()
    body+=text
    topology=open("%s/topology.top"%(address),'w')
    topology.write(body)
    topology.close()
    return 0

def numpy_to_gro(coords,filename):
    new=open(filename+"_new.gro","w")
    old=open(filename+".gro","r")
    old_text=old.read()
    old.close()
    new.write("numpy to gro\n")
    old_text_split=old_text.split("\n")
    new.write(old_text_split[1])
    text_tmp=""
    for i,line in enumerate(old_text_split[2:-2]):
        line_split=line.split()
        new_line="%s%s%s%8.3f%8.3f%8.3f\n"%(line_split[0],line_split[1],line_split[2],
                                            coords[i,0],
                                            coords[i,1],
                                            coords[i,2])
        text_tmp+=new_line
    text_tmp+=" 10 10 10"
    new.write(text_tmp)
    return 0

def mix_gro(filename1,filename2,box):
    gro1=open(filename1,"r")
    gro2=open(filename2,"r")
    gro_text1,gro_text2=gro1.read(),gro2.read()
    gro1.close()
    gro2.close()
    n1,n2=int(gro_text1.split("\n")[1]),int(gro_text2.split("\n")[1])
    new="\n".join(gro_text1.split("\n")[2:-2]+gro_text2.split("\n")[2:-2])
    new="Title\n %d \n %s\n %8.3f %8.3f %8.3f\n"%(n1+n2,new,box[0],box[1],box[2])
    new_file=open("mixer.gro","w")
    new_file.write(new)
    new_file.close()
    return 0

def read_itp(filename="lipids_43A1-S3.itp"):
    itp_file=open(filename)
    itp_txt=itp_file.read()
    itp_file.close()
    data=itp_txt#.replace("#",";")#[itp_txt.find("atoms"):itp_txt.find("bonds")]
    text=""
    for line in data.split("\n"):
        if ";" in line:
            semi=line.find(";")
            if semi:
                text+=line[:semi]
                text+="\n"
        elif line!="":
            text+=line
            text+="\n"
    out=open("raw.itp","w")
    out.write(text)
    out.close()
    res_dict={}
    bond_dict={}
    lines=text.split("\n")
    for i,line in enumerate(lines):
        if "moleculetype" in line:
            residue=lines[i+1].split()[0]
            res_dict[residue]=[]
            bond_dict[residue]=[]
        elif "atoms" in line:
            j=i+1
            if "ifndef" in lines[j]:
                print("Error: C++ preprocessor is detected in the [ atoms ] section in the input files. This means different types of residues are existed in the file which currently is not supported.")
                return None
            while lines[j].split()[0].isnumeric():
                res_dict[residue].append(lines[j].split()[4])
                j+=1
        elif "bonds" in line:
            j=i+1
            if "ifndef" in lines[j]:
                print("Error: C++ preprocessor is detected in the [ atoms ] section in the input files. This means different types of residues are existed in the file which currently is not supported.")
                return None
            while lines[j].split()[0].isnumeric():
                bonds_index=lines[j].split()
                bond_dict[residue].append((bonds_index[0],bonds_index[1]))
                j+=1
    n_residues=len(res_dict.keys())
    print("\n"+"\033[;7m`%s`\033[0;0m\nAnalysis has been finished\n%-3d residue%s found.\n"%(filename,n_residues,[" was", "s were"][n_residues!=1]))
    for key in res_dict:
        print("\033[1;36mresidue\033[0;0m %-6s %-3d atoms"%(key,len(res_dict[key])))
        #print(res_dict[key],"\n")
    return res_dict,bond_dict
