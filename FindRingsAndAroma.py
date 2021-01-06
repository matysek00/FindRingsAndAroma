
import numpy as np
import ase
import ase.io
import ase.neighborlist
import matplotlib.pyplot as plt


class AtomRing():
    """
    A class used to represent rings of atoms
    
    Atributes
    ---------
    Optimal_rs : dict
       size of a ring and a corespoding optimal bond length data taken from On the Harmonic Oscillator Model of Electron Delocalization (HOMED) Index and its Application to Heteroatomic π-Electron Systems
    self.atom_ind : array
       indicies of atoms in the ase.atoms.Atoms system
    self.size : int
       size of the ring
    self.distances : array
       bond lengths [at1-at2, at2-at3 ....., atn-at1]
    self.center : float
       gemetrical central of the ring Not a center of mass
    
    Methods
    -------
    Set_AromHOMA(self, alpha, r_opt=None):
        Calculate HOMA aromaticity
    """
    
    
    Optimal_rs = {3: 1.508,
                  4: 1.554,
                  6: 1.540,
                  7: 1.536,
                  }
    

    def __init__(self, ind, Atoms):
        """
        Parameters
        ----------
        ind : array like
           indicies of atoms in the Atoms system
        Atoms : ase.atoms.Atoms
           system of atoms where the ring is
        """
        self.atom_ind = np.array(ind) #indexes of atom in the AtomRing
        self.size = len(self.atom_ind)
    
        pos = Atoms.get_positions()[self.atom_ind] 
        #distance = [at1-at2, at2-at3 ....., atn-at1]
        self.distances = np.linalg.norm(pos - np.roll(pos, -1, axis=0), axis=1)
        self.center = pos.sum(axis=0)/self.size


    def Set_aromHOMA(self, alpha, r_opt=None):
        """
        AromHOMA(self, alpha, r_opt=None):
        Calculate HOMA aromaticity
        Sets self.aromHOMA = 1 - (alpha/self.size)*((r_opt - self.distances)**2).sum()
        r_opt defaulted to value in AtomRing.Optimal_rs coresponding to the size of AtomRing                                       """

        if r_opt == None:
            r_opt = AtomRing.Optimal_rs[self.size]

        self.aromHOMA = 1 - (alpha/self.size)*((r_opt -self.distances)**2).sum()




def Get_rings(Atoms):
    """
    Get_rings(Atoms):
    Finds all rings in the atoms system and returns them in a list of AtomRing objects.
    
    Parameters
    ----------
    Atoms : ase.atoms.Atoms
      System of atoms where the rings will found
           
    Returns
    -------
    Rings : list
      List of AtomRing object representing all rings in the system
    """
    n= len(Atoms.get_atomic_numbers())
    Neighb =  ase.neighborlist.build_neighbor_list(Atoms, bothways=True)
    
    matrix = Neighb.get_connectivity_matrix(sparse=False) - np.identity(n)
    
    indexes = [np.where(x == 1)[0].tolist() for x in matrix]
    dist_mat = Atoms.get_all_distances()
    
    #now remove all dead ends
    dead_ends = [1]
    
    while len(dead_ends) >= 1:
        dead_ends = [i for i, x in enumerate(indexes) if len(x) == 1]
        #keep empyt lists as place holders
        for j in dead_ends:
            indexes[j] = []
        #remove all refences to dead ends
        indexes = [list(filter(lambda x: x not in dead_ends, ind)) for ind in indexes]
                    
    Rings = []
    
    for i1, neighs1 in enumerate(indexes):
        if len(neighs1) != 0:
            dist1 = dist_mat[i1].copy()#distance of all atoms to the first atom
            
            for i2 in neighs1:
                neighs2 = indexes[i2].copy()
                neighs2.remove(i1)
                        
                for i3 in neighs2:
                    new_Ring = [i1, i2, i3]
                    
                    cont = True       
                    for ring in Rings:
                        #check if the ring does not already exists
                        if i1 in ring and i2 in ring and i3 in ring:
                            cont = False
                            break
                            
                    while cont:
                        
                        paths = indexes[new_Ring[-1]].copy()
                        paths.remove(new_Ring[-2])
                        #posible next atoms excluding the previous one
                        
                        dist = [dist1[i1] for i1 in paths]
                        #their distance to original atom
                        
                        Atnew = paths[dist.index(min(dist))] # newly added atom
                        new_Ring.append(Atnew)
                        cont = new_Ring[-1] not in new_Ring[:-1] # check whether ring closed            
                        
                        if new_Ring[0] == new_Ring[-1] and len(new_Ring) > 4:
                            #check if the ring is okay(at least size 4 and a close loop)
                            Rings.append(new_Ring[:-1])  
    
    Rings = [AtomRing(x,Atoms) for x in Rings]
    return Rings


def PlotRings(Rings, Atoms, inplot=''):
    """
    A simple plot of AtomRings.
    
    Parameters
    ----------
    Rings : list
      List of AtomRing object representing all rings in the system
    Atoms : ase.atoms.Atoms
      System of atoms where the rings will found           
    inplot : str
      name of a variable of AtomRing object that will be writen in the middle of the ring
    """

    for i, ring in enumerate(Rings):
        positions = Atoms[ring.atom_ind].get_positions()
        
        for j in range(ring.size):
            p1 = positions[j][:-1]
            p2 = positions[(j+1) % ring.size][:-1]
            
            plt.arrow(*p1, *(p2-p1),  head_width=0.0, head_length=0.0,  fc='none', color='black')

        try:
            tex = getattr(ring, inplot)
            tex = round(tex, 3) if type(tex) == np.float64 else tex
        except:
            tex = ''

        plt.text(*ring.center[:-1], str(tex), horizontalalignment='center', verticalalignment='center')



