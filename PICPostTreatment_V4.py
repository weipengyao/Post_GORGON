# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 09:50:30 2018

@author: kapichu
"""

import time
import PICclass as PIC
import numpy as np
from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
import random



rc('xtick', labelsize=22)
rc('ytick', labelsize=22)
rcParams['lines.linewidth'] = 4
rcParams['font.size'] = 35
clist = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#7f7f7f', '#e377c2', '#8c564b', '#17becf', '#bcbd22']


###############################################################################
### Data ###
"""Here are all the stuff to declare before running the script. It consists in the name ogf the data files to be analyzed but also
in the details of the simulation box as well as how you want to analyze them """

#C> - Name of the simulation where to find the MHD data
#C>---------------------------------------------------- 
MHD_prefix = "PA77"
MHD_workdir = "/home/kapichu/mesopsl_work/Results/" + MHD_prefix + "/"

#C> - Name of the simulation
#C>-------------------------
prefix = "PA77"
nickname = "PA77"
print "Studying simulation", prefix



#C> - Path to get the data
#C>-----------------------
workdir = "/media/kapichu/data/Work_D2/Weibul/" + prefix + "/PIC_reduced/"


#C> - Path to save the figures
#C>---------------------------
savepath = "/media/kapichu/data/Work_D2/Weibul/" + prefix + "/Treatment/"

###############################################################################

#C> - Required arguments to pass to the PICclass script
#C>----------------------------------------------------
Z = 1.
atomic_weight = 1.


#C> - Informations about times
#C>---------------------------
MHD_Time_start = 0.0
MHD_Time_end = 18.802
MHD_time_out_frequency = 0.5

PIC_Time_start = 2.0                                                          #Time when the PIC module was activated in the simulation; nanosecond unit
PIC_Time_end = 18.802 #17.2                                                    #Time when the PIC module wrote its last output files; nanosecond unit.
PIC_time_out_frequency = 0.1

PIC_iter_start = 0
PIC_iter_stop = int(PIC_Time_end/PIC_time_out_frequency)

###############################################################################

pic = PIC.PIC(workdir, prefix, 0, savepath, Z, atomic_weight)
elementary_charge, gamma, K_Boltzmann, proton_mass, electron_mass, c_light, mu_zero = pic.load_cst()


###############################################################################
####                                Functions                              ####
###############################################################################
#
#
###############################################################################
####                        End of the functions                           ####
###############################################################################




###############################################################################
###                                   MAIN                                  ###
###############################################################################
timing_start = time.time()

pic = PIC.PIC(workdir, prefix, 0, savepath, Z, atomic_weight)

pic.SecurityCheckPaths("")

ipx,ipy,ipz,iPx,iPy,iPz,iBx,iBy,iBz,iEx,iEy,iEz,imass,icharge,itag,iMHD_rho,iMHD_Ti,iMHD_Te,iMHD_Vx,iMHD_Vy,iMHD_Vz,idx_nb_part,idx_PIC_dt,idx_PIC_time = pic.Header()
idx_wrong_dt = idx_nb_part
idx_nb_wrong_dt = idx_PIC_dt


plt.rc('figure.subplot', top=0.94)
plt.rc('figure.subplot', bottom=0.16)
plt.rc('figure.subplot', left=0.14)
plt.rc('figure.subplot', right=0.92)
plt.figure(figsize=(14.,14.))

PIC_selected = np.asarray([])

selection_timing_start = time.time()
for i in range(PIC_iter_start, PIC_iter_stop):
    pic = PIC.PIC(workdir, prefix, i, savepath, Z, atomic_weight)
#
    nb_particles = pic.CountParticles_V2(idx_nb_part)
#
    """Reading all the data for all the particles"""
    data_to_load = np.asarray([ipx, ipy, ipz, iPx, iPy, iPz, iBx, iBy, iBz, iEx, iEy, iEz, imass, icharge, idx_wrong_dt, idx_nb_wrong_dt, itag])
    PICData = pic.ExtractAllData_V2(data_to_load, PIC_time_out_frequency, nb_particles, 15, 16, info=True)  #C> - indices 15 and 16 are for the nb of wrong dt of the particle AND for the tags
#
    """Computing and plotting the energies of all the particles"""
    Energy = pic.ComputeEnergies(PICData[:,3:6],PICData[:,12], PICData[:,-1])
#
    """Creating a sub-category of particles based on a criterion (here: 10.*Eref --> energetic particles) and saving their tags
    WARNING: The function removes the particles selected from the initial PICData array"""
    #PICData, PICData_up = pic.DistinguishEnergy(PICData, Energy, 8.0e6) # 0.*Eref)         #C> - Here is a function using numpy.where function: it is very slow if the number of particles is large
    PICData_up = Energy[(Energy[:,0] > 8.0e4),:]                                            #C> - This way creates a mask and is way faster than the other.
#
    try :
        PIC_selected = np.append(PIC_selected, PICData_up[:,-1])
        print "Number of particles distinguished with the criterion: " + np.size(PICData_up[:,-1])
    except :
        pass

plt.close()
print '  *  '

    
"""Removes multiple appearance of the same particles"""
PIC_selected = np.unique(PIC_selected)
PIC_selected = PIC_selected[1:]
print "PIC first selection: ", np.size(PIC_selected), " particles selected"


selection_timing_end = time.time()
print "All times swept and data read - selection of reduced sets of particles done | Time taken: " + str( "%.2f" % ((selection_timing_end-selection_timing_start)/60.) ) + " minutes"


"""Takes 300 particles randomly in the array of the selected particles"""
random_list = random.sample(xrange(np.size(PIC_selected)), 300)
PIC_selection = np.sort(PIC_selected[random_list])


"""This is an array to save the quantities related to the particles for all the times"""
PICData_sub1 = np.zeros([np.size(PIC_selection),16,PIC_iter_stop+1])         #C> - The '16' is to be changed accordingly to the number of quantities taken
PICData_sub1[:,0,0] = PIC_selection


"""This is an array to save mainly the positions of the selected particles for all the times"""
BasicQ_sub1 = np.zeros([np.size(PIC_selection),14,PIC_iter_stop+1])
BasicQ_sub1[:,0,0] = PIC_selection


PIC_selection = PIC_selection.astype(int)

selection_timing_start = time.time()
for i in range(PIC_iter_start, PIC_iter_stop):
    pic = PIC.PIC(workdir, prefix, i, savepath, Z, atomic_weight)
#C> ----
    nb_particles = pic.CountParticles_V2(idx_nb_part)
    data_to_load = np.asarray([ipx,ipy,ipz,iPx,iPy,iPz,iBx,iBy,iBz,iEx,iEy,iEz,imass,icharge,iMHD_rho,iMHD_Ti,iMHD_Te,iMHD_Vx,iMHD_Vy,iMHD_Vz,idx_wrong_dt,idx_nb_wrong_dt,itag]) #  [ipx, ipy, ipz, iPx, iPy, iPz, iBx, iBy, iBz, iEx, iEy, iEz, imass, icharge, idx_wrong_dt, idx_nb_wrong_dt, itag])
    PICData = pic.ExtractAllData_V2(data_to_load, PIC_time_out_frequency, nb_particles, 21, 22)  #C> - indices 15 and 16 are for the nb of wrong dt of the particle AND for the tags
#C> ----
    max_tag = np.max(PICData[:,-1])
    up_size = np.size(PIC_selection[(PIC_selection < max_tag)])
    if (up_size > 0):
        PICData_first_set = PICData[PIC_selection[(PIC_selection < max_tag)],:]
        print "Time: " + str("%.2f" % (i*PIC_time_out_frequency)) + ", number of particles up found: " + str(np.shape(PICData_first_set))

        PICData_sub1[:up_size, 15, i+1] = PICData_first_set[:, -1]
        PICData_sub1[:up_size, 14, i+1] = PICData_first_set[:, 13]
        PICData_sub1[:up_size, 13, i+1] = PICData_first_set[:, 12]
        """Get and save the spatial coordinates, velocity and B-field components of the particles - for a later use to then study their trajectories and others"""
        BasicQ_sub1[:up_size, :, i+1] = PICData_first_set[:, :14]
        """Get the local (to the particle) MHD quantities: density (kg/m^3), ionic temperature(eV), electronic temperature (eV) and flow velocity (m/s)"""
####    PICData_sub1[:,7,1:] = 1.0e-4                                              #Vacuum of GORGON from 1e-5kg/m3 to 1e-4kg/m3
        PICData_sub1[:up_size, 7:10, i+1] = PICData_first_set[:, 14:17]
        PICData_sub1[:up_size, 10,   i+1] = np.sqrt(PICData_first_set[:, 17]**2. + PICData_first_set[:, 18]**2. + PICData_first_set[:, 19]**2.)
        """Get the local B- and E-fields (in this order!) moduli"""
        PICData_sub1[:up_size, 3, i+1] = np.linalg.norm(PICData_first_set[:, 6: 9], axis=1)
        PICData_sub1[:up_size, 4, i+1] = np.linalg.norm(PICData_first_set[:, 9:12], axis=1)
        """Get the perpendicular and parallel velocities (in this order!) a vector; the array for parallel velocities has four dimensions to save the sign with respect to B_local"""
        Vperp_first_set, Vpara_first_set = pic.V1_wrt_V2(PICData_first_set[:, 3:6], PICData_first_set[:, 6:9], PICData_sub1[:up_size, 0, 0])
        """Get the particle ion temperature (in eV) and then (in this order!) the square of the perpendicular velocity component modulus (in (m/s)^2)"""
        momenta_prime = np.linalg.norm(PICData_first_set[:,3:6], axis=1)*np.linalg.norm(PICData_first_set[:,3:6], axis=1)
        gamma_prime = np.sqrt(1. + (momenta_prime/(PICData_sub1[:up_size, 13, i+1]*c_light)**2))
        PICData_sub1[:up_size, 11, i+1] = (momenta_prime)/(elementary_charge*PICData_sub1[:up_size, 13, i+1]*gamma_prime**2)
        PICData_sub1[:up_size, 12, i+1] = (np.linalg.norm(Vperp_first_set[:, :3], axis=1)*np.linalg.norm(Vperp_first_set[:,:3], axis=1))/(PICData_sub1[:up_size, 13, i+1]*gamma_prime)**2
        """Get the perpendicular and parallel components (in this order!) of the electric field; the array for parallel electric field has four dimensions to save the sign with respect to B_local"""
        Efield_perp_first_set, Efield_para_first_set = pic.V1_wrt_V2(PICData_first_set[:, 9:12], PICData_first_set[:, 6:9], PICData_sub1[:up_size, 0, 0])
        """Save the perpendicular and parallel (in this order!) components of the electric field"""
        PICData_sub1[:up_size, 5, i+1] = np.linalg.norm(Efield_perp_first_set [:, :3],axis=1)
        PICData_sub1[:up_size, 6, i+1] = np.linalg.norm(Efield_para_first_set [:, :3],axis=1)
        """Get the energies: total, perpendicular and parallel (in this order!)"""
        PICData_sub1[:up_size, 0, i+1] = pic.ComputeEnergiesSelection(PICData_first_set[:, 3:6], PICData_sub1[:up_size, 13,i+1], PICData_sub1[:up_size, 15, i+1])[:,0]
        PICData_sub1[:up_size, 1, i+1], PICData_sub1[:up_size, 2, i+1] = pic.EpeEpa(Vperp_first_set, Vpara_first_set, PICData_sub1[:up_size, 13, i+1])
        print " ------------------ "
#C> - -----------------------------------------------------
    print "Time: " + str( "%.2f" % (i*PIC_time_out_frequency) ) + ", data read and quantities derived..."

selection_timing_end = time.time()
print "Reduced set of selected particles: data read and quantities derived... | Time taken: " + str( "%.2f" % ((selection_timing_end-selection_timing_start)/60.) ) + " minutes"
print "*"

selection_timing_start = time.time()


"""Conversion from momenta components to velocity components"""
BasicQ_sub1[:,3:6,1:] = pic.ConversionPj_to_Vj(BasicQ_sub1[:,3:6,1:])


"""Here the script is calling the PIC class deidcated to the plots (see PICclass.py script)"""
"""Plot of the spatial coordinates as function of time and as function of one another"""
picplot = PIC.PICPlot(workdir, nickname, savepath, 14., 14., .9, .1, .09, .96)
print "Calling for the function to plot the trajectories of the particles from the first subcategory"
pic.SecurityCheckPaths("Trajectories/Energetics/")
picplot.PlotTrajectories(PIC_iter_start, PIC_iter_stop, PIC_time_out_frequency, 8., 8., BasicQ_sub1[:,:3,:]*1e3, 50, "Trajectories/Energetics/", 'r-')   #(PIC_Time_start, PIC_Time_end, PIC_time_out_frequency, 8., 8., BasicQ_sub1[:,:3,:]*1e3, 20, "Trajectories/Energetics/", 'r-')


""" Compute and plot some plasma parameters/quantities: cyclotron time, Larmor radius, first adiabatic invariant """
PlasmaQ_1 = np.zeros([np.size(PIC_selected),3,PIC_iter_stop])
PlasmaQ_1[:,0,:] = pic.ComputeTau_ci(PICData_sub1[:,3,1:], PICData_sub1[:,13,1:], PICData_sub1[:,14,1:])
PlasmaQ_1[:,1,:] = pic.ComputeLarmorRadius(PICData_sub1[:,11,1:], PICData_sub1[:,3,1:], PICData_sub1[:,13,1:], PICData_sub1[:,14,1:])
PlasmaQ_1[:,2,:] = pic.ComputeFirstAdiabaticInvariant(PICData_sub1[:,12,1:], PICData_sub1[:,3,1:], PICData_sub1[:,13,1:])


"""Here the script is calling the PIC class deidcated to the plots (see PICclass.py script)"""
picplot = PIC.PICPlot(workdir, nickname, savepath, 14., 14., .9, .1, .14, .96)
print "I'm going to plot the plasma quantities of the energetic particles"
pic.SecurityCheckPaths("Mu_Rho_and_Tau/Energetics/")
picplot.PlotPlasmaQuantities(PIC_iter_start, PIC_Time_end, PIC_time_out_frequency, PICData_sub1[:,0,:], 1.0e3, PlasmaQ_1, 20, "Mu_Rho_and_Tau/Energetics/", 'b-')


selection_timing_end = time.time()
print "All plots relative to the selected particles done | Time taken: " + str( "%.2f" % ((selection_timing_end-selection_timing_start)/60.) ) + " minutes"




#C>========================================================
timing_end = time.time()
print "Performance study with reorganized files - execution duration: " + str((timing_end-timing_start)/60.), " minutes"
print "Simulation studied: " + prefix
print "End"
#C>========================================================




