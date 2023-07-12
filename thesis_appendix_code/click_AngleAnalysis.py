# Usage:  $python3 AngleAnalysis.py path_to_.xyz path_to_.traj path_to_realjumpsorder.npy


import os
import sys
import shutil
import tkinter as tk
from tkinter import filedialog

# from ase.neighborlist import NeighborList as nl#.old
import matplotlib.pyplot as plt
# import pandas
# import seaborn as sns;
# sns.set_theme()
import numpy as np
from ase.io import *
from ase.neighborlist import neighbor_list as nl

# import data files
xyzPresent = False
npyPresent = False

if len(sys.argv) > 1: #if the user inputs the files on the command-line
    print('system arguments detected.')
    #print(len(sys.argv))
    for i in sys.argv:
        #	print(i)
        if '.traj' in i:
            traj = read(str(i), index=":")  # opens the simulation data
        if '.npy' in i:
            jumps = np.load(str(i))
            npyPresent = True
        if '.xyz' in i:
            xyz = open(str(i), "r")
            xyzPresent = True
            
else: #if the user doesn't input the file names (open file manager instead)
    #file manager
    root = tk.Tk()
    root.withdraw()
    file_list = []
    print('no system arguments detected.')
    print('opening file manager window...')
    traj_filename = filedialog.askopenfilename(title = 'Choose the wannier trajectory file',
                                               filetypes = [('trajectory file','.traj')])
    print('file chosen:',traj_filename)
    file_list.append(traj_filename)

    numpy_filename = filedialog.askopenfilename(title = 'Choose the jump list numpy file',
                                               filetypes = [('numpy file','.npy')])
    print('file chosen:',numpy_filename)
    file_list.append(numpy_filename)

    xyz_filename = filedialog.askopenfilename(title = 'Choose the MD xyz file (if necessary)',
                                               filetypes = [('xyz file','.xyz')])
    print('file chosen:',xyz_filename)
    file_list.append(xyz_filename)
    root.destroy()

    #print('\n',file_list, '\n')

    #going through each file in the list to make sure it's the correct type
    for i in file_list:
        #	print(i)
        if '.traj' in i:
            traj = read(str(i), index=":")  # opens the simulation data
            #print(traj[0], '\n')
        if '.npy' in i:
            jumps = np.load(str(i))
            npyPresent = True
        if '.xyz' in i:
            xyz = open(str(i), "r")
            xyzPresent = True


# Functions
# function for calculating analysis progress
def progresspercent(frame, length):
    if (frame % (length / 10) == 0):
        print(str((frame / length) * 100) + "%")
        
if not xyzPresent and npyPresent:
    JumpFrameOffset = (input("Input JumpFrameOffset: "))
    if JumpFrameOffset != 'n':
        JumpFrameOffset = int(JumpFrameOffset)
#JumpFrameOffset = 'n'#(input(Input JumpFrameOffset, to find it from xyz enter NA\n"))
if xyzPresent and npyPresent:
    JumpFrameOffset = 'n'#(input("Input JumpFrameOffset, to find it from xyz enter NA\n"))
    print("Finding JumpFrameOffset")
    xyzFrames = -1
    frameNum = 0
    found = False
    for i in xyz:
        if i == 'ATOMIC_POSITIONS 0 0 0\n':
            xyzFrames = xyzFrames + 1
        for n in traj[0]:
            if str(n.position[0]) in i and str(n.position[1]) in i and str(n.position[2]) in i:
                frameNum = xyzFrames
                print(str(i), str(n), str(frameNum))
                found = True
                break
        if found == True:
            break
    JumpFrameOffset = frameNum
    print(JumpFrameOffset, type(JumpFrameOffset), frameNum)
    JumpFrameOffset = int(JumpFrameOffset)

symPass = True
flipCatOffset = False

print(traj[0].symbols)
cations = ['Ag', 'Cu']
anions = ['I', 'Cl', 'Br']
metadata = str(traj[0].symbols)
(metadata[0])
count = 1
sym = ''
for i in range(count, len(metadata)):
    if metadata[i].isdigit() == False:
        numWC = int(metadata[1:(i)])
        count = i
        break
for i in range(count, len(metadata)):
    if metadata[i].isdigit() == True:
        sym = metadata[count:i]
        count = i
        break
if sym in cations:
    Cation = sym
    if symPass == True:
        symPass = False
elif sym in anions:
    Anion = sym
    if symPass == True:
        flipCatOffset = True
for i in range(count, len(metadata)):
    if metadata[i].isdigit() == False:
        numE1 = int(metadata[count:i])
        count = i
        break
for i in range(count, len(metadata)):
    if metadata[i].isdigit() == True:
        sym = metadata[count:i]
        count = i
        break
if sym in cations:
    Cation = sym
elif sym in anions:
    Anion = sym

numE2 = int(metadata[count:len(metadata)])

print(numWC, numE1, numE2)

if flipCatOffset == True:
    CationIDOffset = numWC + numE1
else:
    CationIDOffset = numWC
CationIDOffset = int(CationIDOffset)
print(CationIDOffset, type(CationIDOffset))

print(Cation,Anion)

str(traj[0][int(numWC + (numE1 + numE2) / 2)].symbol)
Xstart = 0
Xend = numWC - 1
if str(traj[0][int(numWC + (numE1 + numE2) / 2)].symbol) in anions:
    CATstart = numWC
    CATend = numWC + numE1
    ANstart = numWC + numE1 + 1
    ANend = ANstart + numE2 - 1
else:
    ANstart = numWC
    ANend = numWC + numE1
    CATstart = numWC + numE1 + 1
    CATend = CATstart + numE2 - 1

print(Xstart, Xend, ANstart, ANend, CATstart, CATend)

# organize traj file into frames
frames = []  # array for organizing traj file into frames
for atoms in traj:
    frames.append(atoms)

# Set distances for checking neighbors
# ID: X 0-485, I 486-539, Cu 540-593
# can be found in the pair dist functions and can be automated!!!!!!
if 'Ag' in Cation and 'Cl' in Anion:
    radiusCutoff = {('Ag', 'Cl'): 4.0, ('Ag', 'X'): .8, ('Cl', 'X'): 1.0}  # Radius cutoffs changed by NAL
if 'Ag' in Cation and 'Br' in Anion:
    radiusCutoff = {('Ag', 'Br'): 4.0, ('Ag', 'X'): .8, ('Br', 'X'): 1.0}
if 'Cu' in Cation and 'Cl' in Anion:
    radiusCutoff = {('Cu', 'Cl'): 4.0, ('Cu', 'X'): .8, ('Cl', 'X'): 1.0}
if 'Cu' in Cation and 'Br' in Anion:
    radiusCutoff = {('Cu', 'Br'): 4.0, ('Cu', 'X'): .8, ('Br', 'X'): 1.0}  # Radius cutoffs changed by NAL
if 'Ag' in Cation and 'I' in Anion:
    radiusCutoff = {('Ag', 'I'): 4.0, ('Ag', 'X'): .8, ('I', 'X'): 1.0}
if 'Cu' in Cation and 'I' in Anion:
    radiusCutoff = {('Cu', 'I'): 4.0, ('Cu', 'X'): .8, ('I', 'X'): 1.0}  # Radius cutoffs changed by NAL

print(radiusCutoff)


# Angle Analysis by Nicole
print('Initiating...')
# for a in range(len(frames)):#loop through all the frames
allAngleData = []
# for a in range(1):#Test 1 frame
for a in range(len(frames)):
    # progresspercent(a,len(frames))
    centerAtoms, nLIST, dist = nl('ijd', frames[a], radiusCutoff, self_interaction=False)
    
    # Get all angles associated with the anion, I
    if ANend < CATstart: #anion is listed first in the MD
        for anion in range(ANstart,ANend):     #for silver systems (anion first)
            neighs = nLIST[centerAtoms == anion]  # get ALl anion neighbors
            neighsX = neighs[neighs < Xend + 1]  # all WC - Anion neighbors
            neighsCAT = neighs[np.logical_and(neighs < CATend + 1, neighs >= CATstart - 1)]  #for silver systems (anion first)
            
            ##Find closest WC to cation, OR smallest angle. How many angles do we want? All? Start by going by cation:
            for cation in neighsCAT:
                cationarray = np.empty([4, 6])  # make array for angle with each WC
                for WCn, WC in enumerate(neighsX):
                    distbyC = frames[a].get_distance(cation, WC, mic=True)  # mic = true includes periodic boundaries
                    thetabyC = frames[a].get_angle(cation, anion, WC, mic=True)
                    if (WCn < 4):
                        cationarray[WCn] = [distbyC, thetabyC, cation, anion, WC, a]
                cat_sortDist = cationarray[np.argsort(cationarray[:, 0])]
                cat_sortAng = cationarray[np.argsort(cationarray[:, 1])]
                noPref = abs(cat_sortDist[0, 0] - cat_sortDist[1, 0]) < 0.2
                if noPref:
                    allAngleData.append(cat_sortDist[0])
                    allAngleData.append(cat_sortDist[1])
                else:
                    allAngleData.append(cat_sortDist[0])
                    
    else:       #cation is listed first in the MD
        for anion in range(ANstart -1 ,ANend):     #for silver systems (anion first)
            neighs = nLIST[centerAtoms == anion]  # get ALl anion neighbors
            neighsX = neighs[neighs < Xend + 1]  # all WC - Anion neighbors
            neighsCAT = neighs[np.logical_and(neighs < CATend, neighs >= CATstart)]  # for copper systems (cation first)

            ##Find closest WC to cation, OR smallest angle. How many angles do we want? All? Start by going by cation:
            for cation in neighsCAT:
                cationarray = np.empty([4, 6])  # make array for angle with each WC
                for WCn, WC in enumerate(neighsX):
                    distbyC = frames[a].get_distance(cation, WC, mic=True)  # mic = true includes periodic boundaries
                    thetabyC = frames[a].get_angle(cation, anion, WC, mic=True)
                    if (WCn < 4):
                        cationarray[WCn] = [distbyC, thetabyC, cation, anion, WC, a]
                cat_sortDist = cationarray[np.argsort(cationarray[:, 0])]
                cat_sortAng = cationarray[np.argsort(cationarray[:, 1])]
                noPref = abs(cat_sortDist[0, 0] - cat_sortDist[1, 0]) < 0.2
                if noPref:
                    allAngleData.append(cat_sortDist[0])
                    allAngleData.append(cat_sortDist[1])
                else:
                    allAngleData.append(cat_sortDist[0])

allAngleArray = np.asarray(allAngleData)
np.save("./allAngleArrayNLA.npy", allAngleArray)
allAnglePlt = []
for i in range(len(allAngleArray)):
    allAnglePlt.append(allAngleData[i][1])

print(len(allAngleArray))

if os.path.exists('Outputs'):
    shutil.rmtree('Outputs')
    os.makedirs('Outputs')
else:
    os.makedirs('Outputs')

#Data

#All Angles vs count of atoms
plt.xlabel("degrees", fontsize=12)
plt.ylabel("# of atoms", fontsize=12)
plt.hist(allAnglePlt, 100,
         density=False,
         histtype='bar',
         facecolor='b',
         alpha=1)
plt.legend(['Cation - WC - Anion angles'])
plt.savefig('Outputs/' + Cation + Anion + 'AllAngles.png')
# plt.show()
plt.clf()

#All Angles vs Probability density
plt.xlabel("degrees", fontsize=12)
plt.ylabel("Probability density", fontsize=12)
plt.hist(allAnglePlt, 100,
         density=True,
         histtype='bar',
         facecolor='b',
         alpha=1)
plt.legend(['Cation - WC - Anion angles'])
plt.savefig('Outputs/' + Cation + Anion + 'AllAnglesPD.png')
# plt.show()
plt.clf()

#All Distances vs count of atoms
plt.xlabel("angstroms", fontsize=12)
plt.ylabel("# of atoms", fontsize=12)
plt.hist(allAngleArray[:, 0], 100,
         density=False,
         histtype='bar',
         facecolor='b',
         alpha=1)
plt.legend(['Cation - WC - Anion distances'])
plt.savefig('Outputs/' + Cation + Anion + 'AllDistances.png')
# plt.show()
plt.clf()


#All Distances vs probability Density
plt.xlabel("angstroms", fontsize=12)
plt.ylabel("Probability density", fontsize=12)
plt.hist(allAngleArray[:, 0], 100,
         density=True,
         histtype='bar',
         facecolor='b',
         alpha=1)
plt.legend(['Cation - WC - Anion distances'])
plt.savefig('Outputs/' + Cation + Anion + 'AllDistancesPD.png')
# plt.show()
plt.clf()




if npyPresent:
    ##Check jumps
    print('jump analysis')
    # for IDdata, data in enumerate(allAngleArray):
    jumpAngles = []
    jumpframe = []
    jumpDistances = []
    for j in jumps:
        frames = []
        # print(j, (j[3] + int(JumpFrameOffset/10)))
        jumpCAT = j[0] + CationIDOffset
        # jumpframe = [*range((j[3] + int(JumpFrameOffset/10) - 25), (j[3] + int(JumpFrameOffset/10) + 25))]
        # print(jumpframe[0:30])
        jumpmatch = allAngleArray[allAngleArray[:, 2] == jumpCAT]
        # any(item in allAngleArray[:,5] for item in jumpframe))
        for angleROW in allAngleArray:
            if jumpCAT == angleROW[2]:
                if j[3] + (JumpFrameOffset / 100) + 5 > angleROW[5] > j[3] + ((JumpFrameOffset) / 100) - 5:
                    # print(angleROW)
                    # angleRowSort = np.argsort(angleROW, axis= -1, kind=None, order=None)
                    # print(angleROW)
                    frames.append(angleROW)
        # print(len(frames))
        angleRowSort = np.sort(frames, axis=0, kind='quicksort', order=None)
        # print(angleRowSort)
        if len(angleRowSort) > 1 and abs(angleRowSort[0][0] - angleRowSort[1][0]) < .5:
            jumpAngles.append(angleRowSort[0][1])
            jumpAngles.append(angleRowSort[1][1])
            jumpDistances.append(angleRowSort[0][0])
            jumpDistances.append(angleRowSort[1][0])
        elif len(angleRowSort) > 0:
            jumpAngles.append(angleRowSort[0][1])
            jumpDistances.append(angleRowSort[0][0])
    # Jumping cation bond distances vs count
    plt.xlabel("angstroms", fontsize=12)
    plt.ylabel("# of atoms", fontsize=12)
    plt.hist(jumpDistances, 100,
             density=False,
             histtype='bar',
             facecolor='b',
             alpha=1)
    plt.legend(['Cation - WC - Anion distances'])
    plt.savefig('Outputs/' + Cation + Anion + 'JumpDistances.png')
    # plt.show()
    plt.clf()

    # Jumping cation bond distances vs probability density
    plt.xlabel("angstroms", fontsize=12)
    plt.ylabel("Probability density", fontsize=12)
    plt.hist(jumpDistances, 100,
             density=True,
             histtype='bar',
             facecolor='b',
             alpha=1)
    plt.legend(['Cation - WC - Anion distances'])
    plt.savefig('Outputs/' + Cation + Anion + 'JumpDistancesPD.png')
    # plt.show()
    plt.clf()

    # Jumping cation bond angles vs count

    plt.xlabel("angles", fontsize=12)
    plt.ylabel("# of atoms", fontsize=12)
    plt.hist(jumpAngles, 100,
             density=False,
             histtype='bar',
             facecolor='b',
             alpha=1)
    plt.legend(['Cation - WC - Anion angles'])
    plt.savefig('Outputs/' + Cation + Anion + 'JumpAngles.png')
    # plt.show()
    plt.clf()

    # Jumping cation bond angles vs Probability density

    plt.xlabel("angles", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.hist(jumpAngles, 100,
             density=False,
             histtype='bar',
             facecolor='b',
             alpha=1)
    plt.legend(['Cation - WC - Anion angles'])
    plt.savefig('Outputs/' + Cation + Anion + 'JumpAnglesPD.png')
    # plt.show()
    plt.clf()

    #
    # plt.hist([AllAngles, pltAngles], 100,
    labels = ["Angles of diffusing cation interacting with nearest anions",
              "Angles of only jumping cations interacting with nearest anions"]
    plt.xlabel("Degrees", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.hist([allAngleArray[:, 0], jumpAngles], 100,
             density=True,
             histtype='step',
             facecolor='b',
             alpha=1)
    # plt.legend(labels)
    plt.savefig('Outputs/' + Cation + Anion + 'AllvsJump_DistancesPD.png')
    # plt.show()
    plt.clf()

    # plt.hist([AllAngles, pltAngles], 100,
    labels = ["Distances of diffusing cation interacting with nearest anions",
              "Distances of only jumping cations interacting with nearest anions"]
    plt.xlabel("Angstroms", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    plt.hist([allAnglePlt, jumpDistances], 100,
             density=True,
             histtype='step',
             facecolor='b',
             alpha=1)
    # plt.legend(labels)
    plt.savefig('Outputs/' + Cation + Anion + 'AllvsJump_AnglesPD.png')
    # plt.show()
    plt.clf()

    # All angles and All distances
    plt.hist2d(allAngleArray[:, 1], allAngleArray[:, 0], bins=(50, 50), cmap=plt.cm.jet, range=[[0, 90], [0, 3.5]])
    cb = plt.colorbar()
    cb.set_label('Number of entries')
    plt.savefig('Outputs/' + Cation + Anion + 'AllAnglesAllDistances_HeatMap.png')
    # plt.show()
    
    plt.clf()

    # Jump angles & jump distances
    plt.hist2d(jumpAngles, jumpDistances, bins=(25, 25), cmap=plt.cm.jet, range=[[0, 90], [0, 3.5]])
    cb = plt.colorbar()
    cb.set_label('Number of entries')

    plt.savefig('Outputs/' + Cation + Anion + 'JumpAnglesJumpDistances_HeatMap.png')
    # plt.show()
    
    plt.clf()

    # plt.hist([AllAngles, pltAngles], 100,
    labels = ["Angles of diffusing cation interacting with nearest anions",
              "Angles of only jumping cations interacting with nearest anions"]
    plt.xlabel("Degrees", fontsize=12)
    plt.ylabel("Probability Density", fontsize=12)
    n, bins, patsh = plt.hist([allAnglePlt, jumpAngles], 100,
                              density=True,
                              histtype='step',
                              facecolor='b',
                              alpha=1)
    # plt.legend(labels)
    plt.savefig('Outputs/' + Cation + Anion + 'AllvsJumps.png')
    # plt.show()
    plt.clf()




print('Done')
