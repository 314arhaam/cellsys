import os
import sys
import time
import warnings
import utills.gmx
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore",module="matplotlib")
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from scipy.spatial import voronoi_plot_2d
from scipy.spatial.transform import Rotation as rot

force_data = {}
bonds_data = {}
package_banner=r"""

         ██████╗███████╗██╗     ██╗     ███████╗██╗   ██╗███████╗
        ██╔════╝██╔════╝██║     ██║     ██╔════╝╚██╗ ██╔╝██╔════╝
        ██║     █████╗  ██║     ██║     ███████╗ ╚████╔╝ ███████╗
        ██║     ██╔══╝  ██║     ██║     ╚════██║  ╚██╔╝  ╚════██║
        ╚██████╗███████╗███████╗███████╗███████║   ██║   ███████║
         ╚═════╝╚══════╝╚══════╝╚══════╝╚══════╝   ╚═╝   ╚══════╝
                                                                                                                                                                                                                                                                         
Version:     0.1.0 - Jan. 2021
Github:      https://github.com/314arhaam
Paper:       T.B.A!

"""
"""

[Structure]
            Section
                    Stack
                            Atom

Structure  ->  Section  ->  Stack  ->  xyz-coordinates
M.coords      ["upper"]    ["DPPC"]   (n_molecules,n_atoms,3)

"""

class component:
    """An easy tool to create phospholipid membranes

    Some additional functions could be performed automatically on the
    membranes by GROMACS software

        Attributes
        ----------
        nX : int
            Number of monomers in the x-direction

        nY : int
            Number of monomers in the y-direction

        name : str
            Name of the membrane

        coords : dict
            xyz-coordinates of monomers.

        resname : list
            name of the monomers/residues

        ffresname : list
            name of the monomers/residues according to the forcefield's
            atom name.

        monomer : dict
            raw xyz-coordinates for each monomer in the membrane

        data : dict
            forcefield data according to membrane composition.
            actually it is _atom_types for the current component.

        _atom_types : DataFrame
            membrane forcefield database

        composition : dict
            composition of the membrane

        box : list
            bounding box of the unit cell

        filename : str

        n_atoms : list
            number of the atoms
    """
    def __init__(self, forcefield: str = "charmm36") -> None:
        self.forcefield = forcefield
        self.coords, self.monomer, self.composition = {}, {}, {}
        self.residue_name = []
        self.box = [10., 10., 10.]
        self.filename = 'cellsys-%s.gro'
        self.name = "generic-component"
        #self.database=[]#
        #self.data=dict()#
        #self.ffresname=[]#
        #self.number_of_atoms=[]#
        #self.voronoi_plots=[]#
        '''
        if forcefield in [i.split(".")[0] for i in os.listdir('data/forcefield/') if '.csv' in i]:
            ffaddress='data/forcefield/%s.csv'%(forcefield.lower())#gromos43a1_lipid.csv'
            #self._atom_types=pd.read_csv(ffaddress,sep=',')
        else:
            raise KeyError("forcefield %s doesn't exist, available files: %s"%(forcefield, [i for i in os.listdir('data/forcefield/') if '.csv' in i]))
        '''
    
    def __repr__(self) -> str:
        repr = self.name
        for section in self.coords:
            repr += "\n" + " "*4 + section
            for stack in self.residue_name:
                try:
                    n = len(self.coords[section][stack])
                except KeyError:
                    n = 0
                repr += "\n" + " "*4*2 + stack+"%3d"%(n)
        return repr
    
    def load_monomer(self, resname: str) -> None:
        """Load monomer xyz-coords and forcefield atom names to create membrane
        Args:
            rensame (str): residue name or molecule name.
        """
        # avoid duplicate monomer loading
        if resname in self.residue_name:
            print("monomer %s is already existed"%(resname))
            return None
        try:
            # PACK-DATA: read coordinates of monomer/residue from a numpy array
            # default residue structure is coordinates of its atoms.
            # mols=[i for i in os.listdir('data/structure/') if '.npy' in i]
            # res=[mol for mol in mols if resname in mol]
            try:
                numpy_coords=np.load('cellsys_data/'+self.forcefield+'/'+resname+'.npy')
            except FileNotFoundError:
                numpy_coords=np.random.rand(len(force_data[self.forcefield][resname]),3)
                print("randomized")
            # add name of the residue into a list
            self.residue_name.append(resname)
            # save default structure of residue into a dictionary
            self.monomer[resname]=numpy_coords#-numpy_coords.mean(0)
            # put geometrical center of the residue on the origins (0,0,0)
            # not doing this will cause a messy translation in assemble process.
            self.monomer[resname]-=self.monomer[resname].mean(0)
            # ffresname is name of the residues based on the forcefield.
            # same residues may have different names in different forcefields
            # the 2nd row in data.csv is ffresname = _atom_types[resname][0]
            # self.ffresname.append(self._atom_types[resname][0])
            # data file is actually the [atomtypes]
            # self.data[resname]=self._atom_types[resname][1:].dropna()
            # number of atoms in each residue; an important parameter
            # self.number_of_atoms=[self.monomer[res].shape[0] for res in self.residue_name]
            ### finished ###
            print("monomer added: %s\n"%(resname))
            self.__repr__()
        except (KeyError, FileNotFoundError):
            print("404 Error!\nAvailable monomers:")
            print([i[:-4] for i in os.listdir('data/structure/') if '.npy' in i])
        return None
    
    # not finished ; experimental
    def read_gro(self, filename: str) -> int:
        gro_file = open(filename, "r")
        gro_text = gro_file.read()
        gro_file.close()
        data = {}
        data["residue_number"] = []
        data["residue_name"] = []
        data["atom_name"] = []
        data["atom_number"] = []
        data["x"], data["y"], data["z"] = [], [], []
        for line in gro_text.split("\n")[2:-2]:
            data["residue_number"].append(line[:5])
            data["residue_name"].append(line[5:10])
            data["atom_name"].append(line[10:15])
            data["atom_number"].append(line[15:20])
            data["x"].append(float(line[20:28]))
            data["y"].append(float(line[28:36]))
            data["z"].append(float(line[36:44]))
        self.ffresname=data["residue_name"]
        self.coords["generic"] = {}
        self.residue_name.append("whole")
        data = pd.DataFrame(data)
        self.coords["generic"]["whole"] = np.expand_dims(np.array(data[["x", "y", "z"]]), 0)
        return 0
    
    def trim(self, func) -> None:
        new_coords = self.coords.copy()
        new_composition = self.composition.copy()
        for section in self.coords:
            for stack in self.coords[section]:
                mask = func(self.coords[section][stack][:, :, 0],
                            self.coords[section][stack][:, :, 1],
                            self.coords[section][stack][:, :, 2]).all(1)
                new_coords[section][stack] = self.coords[section][stack][mask, :, :]
                new_composition[section] = mask.shape[0]
        self.coords = new_coords
        self.composition = new_composition
        # `self.update` was here
        return None
    
    def toDataFrame(self) -> pd.DataFrame:
        """
        Make a dataframe containing" x,y,z,section,stack,atomname
        """
        data_dict = {}
        data_dict["x"], data_dict["y"], data_dict["z"] = [], [], []
        data_dict["atom"], data_dict["stack"], data_dict["section"] = [], [], []
        for section in self.coords:
            for stack in self.coords[section]:
                for molecule in self.coords[section][stack]:
                    for i,atom in enumerate(molecule):
                        data_dict["x"].append(atom[0])
                        data_dict["y"].append(atom[1])
                        data_dict["z"].append(atom[2])
                        data_dict["atom"].append(force_data[self.forcefield][stack][i+1])#self.data[stack][i+1])
                        data_dict["stack"].append(stack)
                        data_dict["section"].append(section)
        return pd.DataFrame(data = data_dict)
    
    def make_grid(self, min_dist: float, nX: int, nY: int) -> np.array:
        """Mathematically a membrane is an n*m matrix of monomers.
        this method makes a nX*nY grid with equal manhattan-distance `min_dist`
        between 2 adjacent elements.

        Args:
            min_dist (float): distance between neighbouring nodes

        Returns:
            grid_matrix (np.array): a matrix

        Note:
            In future versions, z-axis will be considered.
        """
        grid_matrix = [[i*min_dist,j*min_dist,0]  for i in range(nX)
                       for j in range(nY)]
        grid_matrix = np.array(grid_matrix)
        return grid_matrix
    
    def all_coords(self) -> np.array:
        """Show xyz-coords of all atoms as an M*N matrix

        Returns:
            np.array
        """
        coords = []
        for section in self.coords:
            for stack in self.coords[section]:
                x, y, _ = self.coords[section][stack].shape
                coords.append(self.coords[section][stack].reshape(x*y,3))
        return np.concatenate(coords, axis = 0)
    
    def set_coords(self, new_coords) -> None:
        """Change the values of `coords` attribute to new_coords

        Args:
            new_coords (np.array)
                new_coords.shape = total number of atoms , 3
                unlike the self.coords, this is a 2D array. to have monomer-wise
                coordinates, new_coords must be divided into different chunks.
                each chunk has the shape of:
                (# of atoms in monomer[i] * # of monomer[i] molecules in
                membrane,3) composition of upper and lower section are the same
                so first element of each chunk.shape is an even number,
                indicating that exact number of the monomers in the upper-
                section exist in lower section.
        """
        coords_temp = {}
        coords_temp.fromkeys(self.coords)
        before = 0
        for k,resname in enumerate(self.residue_name):
            for j in self.coords: #self.coords
                #chunk_len=self.composition[j][k]*self.number_of_atoms[k] 
                chunk_len = self.composition[j][k]*len(force_data[self.forcefield][self.residue_name[k]])
                x, y, _ = self.coords[j][resname].shape
                self.coords[j][resname] = (new_coords[before:before+chunk_len,:]).reshape(x,y,3)
                before += chunk_len
        # `self.update` was here
        return None
    
    def fit_box(self, delta_h: float = 5.) -> None:
        """Fit a simulation box
        """
        c = []
        section = list(self.coords.keys())[0]
        for stack in self.coords[section]:
            temp = self.coords[section][stack]
            nX, nY, _ = temp.shape
            c.append(temp.reshape(nX*nY, 3))
        m = np.concatenate(c, axis = 0)
        self.box = m.max(0) - m.min(0)
        self.box[2] = delta_h * 1.5
        '''
        mean_box=self.all_coords().mean(0)
        vec=self.box/2-mean_box
        self.translate(vec)
        '''
        return None
    
    def write_gro(self) -> None:
        # `self.update` was here
        filename=self.filename#%(self.name)
        print("writing data to %s"%(filename))
        box=" %8.3f %8.3f %8.3f\n"%tuple(self.box)
        print("box = ", box, self.box)
        count,text_=1,""
        #self.set_coords(self.all_coords()-self.all_coords().mean(0)+self.box/2)# test
        self.translate(-self.all_coords().mean(0)+self.box/2)
        for stack_index,stack in enumerate(self.residue_name):
            for section in self.coords:
                for residue_index,residue in enumerate(self.coords[section][stack]):
                    for atom_index,atom in enumerate(residue):
                        text_+="%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n"%(1+residue_index,
                                                                 self.residue_name[stack_index],
                                                                 force_data[self.forcefield][stack][atom_index],#self.data[stack][atom_index+1],
                                                                 count,
                                                                 atom[0],atom[1],atom[2])
                        count+=1
        gro_file=self.name+" Built by CellSys\n"+str(count-1)+'\n'+text_+box
        address="CellSys-%s/%s"%(filename[:-4],filename)
        final_membrane=open(address, 'w')
        final_membrane.write(gro_file)
        final_membrane.close()
        return None
    ###
    def translate(self,vector):
        """Translate all monomers-based cartesian translation

        Args:
            vector: (np.array) a 1x3 numpy array which is translation vector.
        """
        for section in self.coords:
            for stack in self.coords[section]:
                self.coords[section][stack]+=vector#test
        return 0
    def vmd(self):
        address="CellSys-%s/%s"%(self.filename[:-4],self.filename)
        os.system("vmd %s"%(address))
        return 0

class membrane(component):
    def __init__(self,nX,nY,forcefield):
        super().__init__(forcefield)
        self.nX,self.nY=nX,nY
        self.voronoi_plots=[]
    def apl(self,h):
        """
        """
        self.fit_box(delta_h=h)
        return np.prod(self.box[:2])/(self.nX*self.nY)*100.0
    def assemble(self,vector,comp=None,z=0):
        """Assemble a monolayer according to input data
        Args:
            vector: positions
            comp: composition
            z: membrane thickness

        Returns:
            coords: coordinates of atoms in a section assembly
        """
        if comp == None:
            comp = [1*int(self.nX*self.nY)]
        coords, xyz = {}, 0
        index = [*range(self.nX*self.nY)]
        np.random.shuffle(index)
        vector=vector[index]
        for i, n in enumerate(comp):
            if n:
                res = self.residue_name[i]
                #I,J=np.ones([self.number_of_atoms[i],1]), np.ones([n,1]) # len(force_data[self.forcefield][self.residue_name[k]])
                I,J = np.ones([len(force_data[self.forcefield][self.residue_name[i]]),1]), np.ones([n,1])
                xyz = np.kron(J,self.monomer[res])+np.kron(vector[:n],I)
                #xyz=xyz.reshape(n,self.number_of_atoms[i],3) # len(force_data[self.forcefield][self.residue_name[k]])
                xyz = xyz.reshape(n,len(force_data[self.forcefield][self.residue_name[i]]),3)
                coords[self.residue_name[i]]=xyz
                vector = vector[n:, :]
            else:
                coords[self.residue_name[i]] = np.zeros((0,0,3))
        self.name = "[%s]-[PY]"%('-'.join(self.residue_name))
        if z:
            temp = {}
            for resname in coords:
                temp[resname] = (coords[resname] - np.array([[[0, 0, z]]]))*np.array([[[1, 1, -1]]])
            coords = temp
        return coords
    def voronoi(self):
        for section in self.coords:
            cg = np.concatenate(tuple(self.coords[section][resname].mean(1)[:, :2]
                                      for resname in self.residue_name))
            self.voronoi_plots.append(Voronoi(cg))
        return self.voronoi_plots
    def membrane_map(self):
        number_of_sections=len(self.voronoi_plots)
        _,ax = plt.subplots(number_of_sections)
        if number_of_sections == 1:
            ax = [ax]
        for i, vplot in enumerate(self.voronoi_plots):
            voronoi_plot_2d(vplot, ax = ax[i])
        plt.show()
        return 0
    
    def compress(self,vector = np.array([[0.1, 0.1, 0]])):
        for section in self.coords.keys():
            for i, stack in enumerate(self.coords[section]):
                if self.composition[section][i]:
                    coords = self.coords[section][stack]
                    coords -= np.expand_dims(coords.mean(1), 1) * vector
                    self.coords[section][stack] = coords
        # `self.update` was here
        return 0

class monolayer(membrane):
    def __init__(self,nX,nY,forcefield) -> None:
        super().__init__(nX, nY, forcefield)
    
    def make(self, vector, comp, name = "main_layer") -> None:
        if type(vector) in {float, int}:
            vector = self.make_grid(vector, self.nX, self.nY)
        self.composition[name] = comp
        self.coords[name] = self.assemble(vector,comp, z = 1)
        self.fit_box(delta_h=5.0)
        # `self.update` was here
        n = ''
        for i in range(len(self.residue_name)):
            N = 0
            for j in self.coords:
                N += self.composition[j][i]
            n = n + str(N) + '-'
        self.filename = f"{['[P]-', '[M]-'][len(self.residue_name)>1]+self.name+n[: -1]:s}.gro"
        os.system(f"mkdir CellSys-{self.filename[:-4]:s}")
        return None

class bilayer(membrane):
    def __init__(self, nX, nY, forcefield):
        super().__init__(nX, nY, forcefield)
    
    def apl(self, h):
        """
        """
        thickness=self.thickness()
        self.fit_box(delta_h=thickness)
        return np.prod(self.box[:2])/(self.nX*self.nY)*100.0
    
    def make(self,vector, z_dist, comp_upper=None, comp_lower=None):
        if type(vector) in {float,int}:
            vector=self.make_grid(vector,self.nX,self.nY)
        if sum(comp_lower)!=self.nX*self.nY or sum(comp_upper)!=self.nX*self.nY:
            warnings.warn("There is a hole in the membrane.")
        self.composition["upper"]=comp_upper
        self.composition["lower"]=comp_lower
        self.coords["upper"]=self.assemble(vector, comp_upper)
        self.coords["lower"]=self.assemble(vector, comp_lower, z_dist)
        h=self.thickness()
        self.fit_box(delta_h=h)
        # `self.update` was here
        n = ""
        for i, _ in enumerate(self.residue_name):
            N = 0
            for j in self.coords:
                N += self.composition[j][i]
            n = n + str(N) + '-'
        self.filename=['[P]-','[M]-'][len(self.residue_name)>1]+self.name+n[:-1]+'.gro'
        os.system("mkdir CellSys-%s"%(self.filename[:-4]))
        return None
    
    def thickness(self,ref=5):
        z_value=dict()
        z_value["upper"],z_value["lower"]=[],[]
        lipids=[resname for resname in self.residue_name if resname not in ["CHOL"]]
        for resname in lipids:
            for monomer in self.coords["upper"][resname]:
                for sweep in lipids:
                    if self.coords["lower"][sweep].shape[0]:
                        distance=monomer[ref,:]-self.coords["lower"][sweep][:,ref,:]
                        nearest=(distance[:,:2]**2).sum(1).argmin(0)
                        z_value["upper"].append(np.array(distance[nearest,2]))
            for monomer in self.coords["lower"][resname]:
                for sweep in lipids:
                    if self.coords["upper"][sweep].shape[0]:
                        distance=monomer[ref,:]-self.coords["upper"][sweep][:,ref,:]
                        nearest=(distance[:,:2]**2).sum(1).argmin(0)
                        z_value["lower"].append(np.array(distance[nearest,2]))
        return abs(np.append(np.array(z_value["upper"]),np.array(z_value["lower"]))).mean()

class micelle(component):
    def __init__(self,N,M):
        super().__init__()
        self.name='micelle'
        self.N=N
        self.M=M
    def assemble(self):
        """Assemble a monolayer according to input data
        Args:
            vector: positions
            comp: composition
            z: membrane thickness

        Returns:
            coords: coordinates of atoms in a section assembly
        """
        N,M=self.N,self.M
        coords,xyz=dict(),0
        res=self.residue_name[0]
        mono=self.monomer[res]*np.array([[1,1,-1]])+np.array([[0,0,5]])
        xyz=np.array([mono,]*N*M)
        m=rot.from_euler('zx',np.array([[i/N*2*np.pi, (1+j)/(M+1)*np.pi] for i in range(N) for j in range(M)]))
        xyz=np.matmul(xyz,m.as_matrix())
        coords[self.residue_name[0]]=xyz
        self.name=[self.name,"[%s]-[PY]"%('-'.join(self.residue_name))][self.name=='Membrane']
        self.coords["all"]=coords
        self.composition["all"]=[N*M]
        return coords
    def compress(self,vector=np.array([[0.1,0.1,0.1]])):
        """
        """
        for section in self.coords.keys():
            for i,stack in enumerate(self.coords[section]):
                if self.composition[section][i]:
                    coords=self.coords[section][stack]
                    coords-=np.expand_dims(coords.mean(1),1)*vector
                    self.coords[section][stack]=coords
        return 0


if __name__=="cellsys":
    os.system("clear")
    print(package_banner)#isometric2
    print("\nAvailable Forcefields: ")
    forcefield_folder = os.listdir("cellsys_data")
    for i, forcefield in enumerate(forcefield_folder):
        print("%2d) \033[1;32m%s\033[0;0m"%(i+1, forcefield))
    for forcefield in forcefield_folder:
        force_data[forcefield] = {}
        bonds_data[forcefield] = {}
        for f in os.listdir(f"cellsys_data/{forcefield:s}"):
            if f[-4:] == ".itp":
                res_analysis, bond_analysis = utills.gmx.read_itp(f"cellsys_data/{forcefield:s}/{f:s}")
                force_data[forcefield].update(res_analysis)
                bonds_data[forcefield].update(bond_analysis)
    print("\nForcefield data")
    for forcefield in force_data:
        print("\033[1;32m%-15s\033[0;0m"%(forcefield) ,end="")
        print("\033[1;34m", end = "")
        print(*force_data[forcefield].keys(), sep = " ")
        print("\033[0;0m", end = "")

