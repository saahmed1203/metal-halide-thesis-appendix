
# coding: utf-8

# In[1]:


print("Running file")
import sys, os
import numpy as np
from datetime import datetime
import itertools

from tqdm import tqdm_notebook, tnrange, tqdm

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D

import ase
from ase.visualize import view
from ase import Atoms
import ase.neighborlist

from quippy import Atoms as qpAtoms
import quippy as qp

from sitator.util import PBCCalculator
from sitator.visualization import plot_atoms

from samos.trajectory import Trajectory
from samos.analysis.rdf import RDF
from samos.plotting.plot_rdf import plot_rdf


# In[ ]:


#can make traj files with gui
#can only view traj files by loading it onto an array


# In[ ]:


#traj file is a list version of xyz file
#i will use traj file as practice input to understand this code   2/20/2020


# In[10]:


T = 600
#contains all the frames, coordinates, atom types, positions from a simulation
traj_fname='./total_AgI600beta.xyz'
#contains ??
traj_samos='./AgIinput.samos' 

#opens traj_samos file
if os.path.exists(traj_samos):
    t = Trajectory.load_file(traj_samos)
#else loop changes with diff traj file                             
else:
    #finding the atomic species and the number of atoms automatically
    xyz_file = open(traj_fname,'r')

    cations = ['Ag', 'Cu']
    anions = ['I', 'Cl', 'Br']

    num_ion_1 = 0
    num_ion_2 = 0
    cation_is_first_atom = False

    #look at just the first frame to count the number of species
    for ind, line in enumerate(xyz_file):
        if ind >= 2: #3rd row
            line = line.split()
            #start counting the number of ion 1
            if line[0] in cations:
                if ind == 2 or cation_is_first_atom:
                    cation_is_first_atom = True
                    num_ion_1 += 1
                    sym_ion_1 = line[0]
                else:
                    num_ion_2 += 1
                    sym_ion_2 = line[0]


            #start counting the number of ion 2    
            elif line[0] in anions:
                if cation_is_first_atom:
                    num_ion_2 += 1
                    sym_ion_2 = line[0]
                else:
                    num_ion_1 += 1
                    sym_ion_1 = line[0]

            if "ATOMIC_POSITIONS" in line:
                break

   # print('Ion 1:', sym_ion_1, '; number of that ion:',num_ion_1)
   # print('Ion 2:', sym_ion_2, '; number of that ion:',num_ion_2)
   
    #rest of old code
    traj_dict = {
            #what are these bottom numbers? A: lattice vectors for LLZO(specific to every structure)
            "cell":[[13.91120388,0.0,0.0],[0.0,16.06327461,0.0],[0.0,0.0,15.39067040]],  #making a matrix to represent lattice vectors
            #
            "trajectory_format":"axsf", 'timestep':2.9, #time step changes with diff files too
            #each atom is atom object, w/ position, index, ect...to build array struc
            "species":[sym_ion_1]*num_ion_1+[sym_ion_2]*num_ion_2 #(i opened xyz file and counted each I and Ag. is there another way to find out #of species?)
        }

    
    #nat stands for: number of atoms
    #print nat
    
    nat = len(traj_dict['species'])  #species represents all the atoms in the simulation
    
    print("nat:", nat)
    
    header_per_frame = 2 #(input file has 2 info lines before)
    nlines_per_frame = nat + header_per_frame #lines per frame? lines of what? 
    positions = list() #creating empty list to append calculated positions
    trajectory_format = traj_dict['trajectory_format']
    
    #OPENS and READS traj file~~~~~~~~~~~~~~~~~~~~~~~~(might have to fix!!!)
    with open(traj_fname) as ftraj: 
        #what is istep? A: iterating through frames????????
        for istep in range(0, traj_dict.get('skip_steps',0)): 
            [ftraj.readline() for l in range(nlines_per_frame)]
        while True: #while WHAT is true?
            sth_ = [ftraj.readline() for _ in range(header_per_frame)][-1] #defining every frame
            #sth stands for:
#             If line is not header frame:
            if not sth_:
                break
            positions_t = np.empty((nat, 3)) 
            
            for iat in range(nat):
                line = ftraj.readline().split()
                print(line)#PRINT LINE
                positions_t[iat, :] = line[1:4]  
            positions.append(positions_t) #appending positions into previously empty positions list
    atoms = Atoms(traj_dict['species']) #load species onto Atoms class
    atoms.set_positions(positions[0])
    atoms.cell = np.array(traj_dict['cell']) #idk what this line is doing A: creating unit cell?????
    t = Trajectory() #i don't understand what Trajectory() is #
    t.set_timestep(traj_dict['timestep'])
    t.set_atoms(atoms)
    t.set_positions(np.array(positions))
    t.atoms.set_pbc(True)
    t.recenter([sym_ion_1,sym_ion_2]) #change this for diff traj file
    t.save(traj_samos)
    t.write(traj_samos)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

 #info from bottom line comes from traj file   
print "Trajectory length: %i steps (%f ps at dt = %f fs)" % (t.nstep, t.nstep * t.get_timestep() / 1000.0, t.get_timestep())


#input for this cell: trajectory file(output from md)
#output for this cell: breaks file into lists(to initialize list)
#Output: steps, ps(length of simulation), dt(time step)


# In[1]:


# t.atoms.get_pbc()
# t.atoms.set_pbc(True)
# t.atoms.get_cell()


# In[2]:


#t.atoms.get_positions()

