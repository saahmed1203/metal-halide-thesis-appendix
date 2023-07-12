import sys, os
#os.chdir("/home/shared/Li3OBr_Cl/xdats/1vacSM/1100/")
sys.version_info

import sys, os
import numpy as np
import numpy.ma as ma
from math import *
import ase
from ase import *
from ase.io import *
import sitator
from sitator import SiteNetwork
from sitator import SiteTrajectory
from samos.trajectory import Trajectory
from samos.analysis import dynamics
reload(dynamics)
from samos.plotting.plot_dynamics import plot_msd_isotropic
from ase.data.colors import jmol_colors
from sitator import voronoi
from sitator.voronoi import VoronoiSiteGenerator as VSG
from sitator import site_descriptors as sd
from sitator.site_descriptors import SiteTypeAnalysis as STA
from sitator.dynamics import JumpAnalysis as JA
from sitator.dynamics import MergeSitesByDynamics as MSBD
from sitator.misc import GenerateAroundSites as GAS
from sitator.misc import NAvgsPerSite as NAPS
from sitator.landmark import LandmarkAnalysis
from sitator.util import PBCCalculator
from sitator.visualization import plot_atoms, layers, grid, plot_points
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
sys.path.append("/home/shared/Python_Scripts/lethal_lithium")
sys.path.append("/home/shared/Python_Scripts/lethal_lithium/carlo/codes/voronoi")
#sys.path.append("/home/salma/angle_analysis/jump_analysis/AgCl800rs")
from carlo.codes.voronoi import decompositionNN
from md_analysis import jumpanalysis, rotations


# Make real jump array
def makesitearrays(site_traj):
    Ns = len(site_traj)
    alljumps = []
    backjumps = []
    realjumps = []
    sitearray = site_traj
    for i in range(len(sitearray[0,:])):
#    for i in range(14,15):
	# for each Li, find the first frame where a site is assigned 
	# find frame where site does not equal -1
        if sitearray[0,i]==-1:
            Lisite = sitearray[np.argmax(sitearray[0:,i]>-1),i]
        else:
            Lisite = sitearray[0,i]
        prevSITE = Lisite
	# now the code looks for backjumps
        alljumpsi = []
        backjumpsi = []
        realjumpsi = []
#        for n, s in enumerate(sitearray[:100,i]):
        for n, s in enumerate(sitearray[:,i]):
            if s != Lisite and s!= -1:
                jumpi = [i,Lisite,s,n-1]
                Lisite = s
		backjumpBOOL = False
                if np.any(sitearray[0:n-1,i]==s):
                    repeatfr = np.max(np.argwhere(sitearray[0:n-1,i]==s))
                    if n - repeatfr < 50 :
                        firstbackj = [x for x in alljumpsi if x[0]==i and x[1]==s][-1]
                        if firstbackj not in backjumpsi:
                            backjumpBOOL = True
                            backjumpsi.append(firstbackj)
                            backjumpsi.append(jumpi)
                if backjumpBOOL == True:
                    realjumpsi.remove(firstbackj)
                else:
                    realjumpsi.append(jumpi)
                alljumpsi.append(jumpi)
        alljumps = alljumps + alljumpsi
        backjumps = backjumps + backjumpsi
        realjumps = realjumps + realjumpsi
    realjumpsarray = np.array(realjumps)
    backjumpsarray = np.array(backjumps)
    alljumparray = np.array(alljumps)
    backjumpsorder = backjumpsarray[backjumpsarray[:,3].argsort()]
    flatjumpsorder = alljumparray[alljumparray[:,3].argsort()]
    realjumpsorder = realjumpsarray[realjumpsarray[:,3].argsort()]
    np.save("flatjumpsorder.npy",flatjumpsorder)
    np.save("backjumpsorder.npy",backjumpsorder)
    np.save("realjumpsorder.npy",realjumpsorder)
    return flatjumpsorder, backjumpsorder, realjumpsorder

if __name__=="__main__":
# Load data from previous sitator run
    site_traj = np.load("AgI_traj.npy")
    Ns = len(site_traj)
    print(Ns)
#  Here is where you make your jumps lists
    alljumpsorder, backjumpsorder, realjumpsorder = makesitearrays(site_traj) 
 #   realjumpsorder = np.load("realjumpsorder.npy")
 #   backjumpsorder = np.load("backjumpsorder.npy")
 #   alljumpsorder = np.load("flatjumpsorder.npy")
    print len(backjumpsorder), len(realjumpsorder), len(alljumpsorder)
