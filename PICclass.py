#import PVTIclass as PVTI
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 16 16:18:44 2018

@author: kapichu
"""

#import PVTIclass as PVTI
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import time
import sys
import os
import csv


#from pyexp import normalize, dot, cross

###         CONSTANTS       ###
elementary_charge = 1.6022e-19
gamma = 5./3.
K_Boltzmann = 1.3807e-23
proton_mass = 1.6726e-27
electron_mass = 9.1094e-31
c_light = 2.9979e8/20.
mu_zero = 4.0*np.pi*1e-7
###         CONSTANTS       ###

class PIC:
    
    def __init__(self, wrkdir, pfx, tm, svpath, chrge, atm_mass, fname=False):
        self.workdir = wrkdir
        self.prefix = pfx
        self.time = tm
        self.savepath = svpath
        self.charge = chrge
        self.atomic_weight = atm_mass
        self.filename=fname
    
    def load_cst(self):
        
        qe = 1.6022e-19
        Gamma = 5./3.
        K_B = 1.3807e-23
        mp = 1.6726e-27
        me = 9.1094e-31
        c_light = 2.9979e8/20.
        mu_zero = 4.0*np.pi*1e-7
        
        return qe, Gamma, K_B, mp, me, c_light, mu_zero


###############################################################################
#
#
###############################################################################
    def WriteParaviewFiles(self, t_start, t_stop, array_to_write, file_to_save_to):
        """Function that writes files Paraview can read - it takes as arguments an array 3-columns array ([particle,1:3,1:]) containing the spatial coordinates of the particles whose
           tags are saved on the first line of the array ([particle,0,0]).
           WARNING: That may be too complicated, maybe I should change it for a simpler solution.
           -------------------------------------------------------------------------------------------
           That version writes as many files as PIC times and saves all the particles at the corresponding time - the name structure is 'file_name.csv.time'
           That way, when browsing the files - to make an animation for example - Paraview opens the file when it needs it and closes it when it is done with it
           - it just sweeps over them in time - which allows not to overload the memory; this is thanks to the structure of the file names.
           -------------------------------------------------------------------------------------------
           Another name structure like: 'file_name_time.csv' cannot work because Paraview does not sweep over the files as it advances in time, it loads them
           in row which implies a risk of memory overloading at some point."""
#
        for l in range(0,np.size(array_to_write[:,0,0])):
            try:
                tmp_idx = np.where(array_to_write[l,0,1:] == 0.0)[0][0]

                array_to_write[l,0,tmp_idx:] = array_to_write[l,0,tmp_idx]
                array_to_write[l,1,tmp_idx:] = array_to_write[l,1,tmp_idx]
                array_to_write[l,2,tmp_idx:] = array_to_write[l,2,tmp_idx]

            except:
                pass
#
        for i in range(t_start, t_stop):
            file_name = self.savepath + file_to_save_to + ".csv." + str(i)
#
            with open(file_name, 'wb') as csvfile:
                spamwriter = csv.writer(csvfile, delimiter=' ', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
                spamwriter.writerow(['#X-coordinate'] + ['Y-coordinate'] + ['Z-coordinate'] + ['Particle-tag'])
                for j in range(0,np.size(array_to_write[:,0,0])):
                    #if(array_to_write[j,0,i+1]!=0.0 and array_to_write[j,1,i+1]!=0.0 and array_to_write[j,2,i+1]!=0.0) :
                    #    pass
                    spamwriter.writerow([array_to_write[j,0,i+1], array_to_write[j,1,i+1], array_to_write[j,2,i+1], array_to_write[j,0,0]])
#
        return



###############################################################################
    def SecurityCheckPaths(self, dir_to_check):
        """Function checking the existence of the directories where the data are said to be stored and where the figures
        are said to be saved, that is the main directory and all the subdirectories of the tree according to the flags.
        -
        WARNING: For now it works without any arguments because it is all in the same script, but when going to a
                 modular way of working, it will be needed to add the arguments. """
        #
        if not os.path.exists(self.workdir):
            sys.exit("Directory to get the data does not exist... Stopping execution.")
        #
        if not os.path.exists(self.savepath + dir_to_check):
            os.makedirs(self.savepath + dir_to_check)
            print("SecurityCheckPaths function: creation of the directory: " + self.savepath + dir_to_check + ".")
        
        return
###############################################################################
#
#
###############################################################################
    def Header(self):
        """Gets the indices of the quantities by reading the header of a PIC data file.
        --> 3 for particle position, 3 for particle momentum components, 3 for E-field components, 3 for B-field components, 
            1 for the tag, 1 for the number of particles, 1 for PIC time step and 1 for the PIC time."""
#
        src = open(self.workdir+"PIC_output0.dat", "r")
        tmp_head = src.readline().rstrip('\n\r').split()
#
        idx_pos_x = tmp_head.index('#position_x')
        idx_pos_y = tmp_head.index('position_y')
        idx_pos_z = tmp_head.index('position_z')
#
        idx_imp_x = tmp_head.index('impulsion_x')
        idx_imp_y = tmp_head.index('impulsion_y')
        idx_imp_z = tmp_head.index('impulsion_z')
#
        idx_B_x = tmp_head.index('B_x')
        idx_B_y = tmp_head.index('B_y')
        idx_B_z = tmp_head.index('B_z')
#
        idx_E_x = tmp_head.index('E_x')
        idx_E_y = tmp_head.index('E_y')
        idx_E_z = tmp_head.index('E_z')
#
        idx_mass = tmp_head.index('mass')
        idx_charge = tmp_head.index('charge')
        idx_tag = tmp_head.index('tag')
#
        idx_MHD_rho = tmp_head.index('MHD_density')
        idx_MHD_Ti = tmp_head.index('MHD_Ti')
        idx_MHD_Te = tmp_head.index('MHD_Te')
        idx_MHD_Vx = tmp_head.index('MHD_Vx')
        idx_MHD_Vy = tmp_head.index('MHD_Vy')
        idx_MHD_Vz = tmp_head.index('MHD_Vz')
#
        idx_nb_part = idx_MHD_Vz + 1
        idx_PIC_dt = idx_MHD_Vz + 2
        idx_PIC_time = idx_MHD_Vz + 3

        return idx_pos_x,idx_pos_y,idx_pos_z,idx_imp_x,idx_imp_y,idx_imp_z,idx_B_x,idx_B_y,idx_B_z,idx_E_x,idx_E_y,idx_E_z,idx_mass,idx_charge,idx_tag, idx_MHD_rho,idx_MHD_Ti,idx_MHD_Te,idx_MHD_Vx,idx_MHD_Vy,idx_MHD_Vz, idx_nb_part,idx_PIC_dt,idx_PIC_time
        #return idx_pos_x,idx_pos_y,idx_pos_z,idx_imp_x,idx_imp_y,idx_imp_z,idx_B_x,idx_B_y,idx_B_z,idx_E_x,idx_E_y,idx_E_z,idx_tag, idx_MHD_rho,idx_MHD_Ti,idx_MHD_Te,idx_MHD_Vx,idx_MHD_Vy,idx_MHD_Vz, idx_nb_part,idx_PIC_dt,idx_PIC_time
###############################################################################
#
#
###############################################################################
    def CountParticles(self, ind_npart):
        """Function that counts the total number of particles in the whole domain at a given time
        --> Is used to pre-allocate the arrays to save computation time (optimization purpose)."""
#
        n_part = 0
        self.filename = self.workdir + "PIC_output" + str(self.time) + ".dat"
        dfile = open(self.filename, "r")
        tmp_head = dfile.readline().rstrip('\n\r').split()
        n_part = n_part + int(tmp_head[ind_npart])
#
        return n_part
###############################################################################
#
#
###############################################################################
    def ExtractAllData(self, datarray):
        """Function that extracts all the data about all the particles (positions, momenta, fields, tags) without any selection criterion/criteria.
        Is supposed to be used alongside CountParticles although not mandatory. Takes as argument an array of size the total number of particles
        and reads all the files from all the processors at a given time to fill the array with all the features of the array without any distinction.
        --> Not optimized at all."""
#
        tmp_ctr = 0
#
        self.filename = self.workdir + "PIC_output" + str(self.time) + ".dat"
        dfile = open(self.filename, "r")
        tmp_head = dfile.readline().rstrip('\n\r').split()
#
        for line in dfile:
            line = line.rstrip('\n\r').split()
#
            datarray[tmp_ctr,0] = float(line[0])        #ipx])
            datarray[tmp_ctr,1] = float(line[1])        #ipy])
            datarray[tmp_ctr,2] = float(line[2])        #ipz])
            datarray[tmp_ctr,3] = float(line[3])        #iPx])
            datarray[tmp_ctr,4] = float(line[4])        #iPy])
            datarray[tmp_ctr,5] = float(line[5])        #iPz])
            datarray[tmp_ctr,6] = float(line[9])        #iBx])
            datarray[tmp_ctr,7] = float(line[10])       #iBy])
            datarray[tmp_ctr,8] = float(line[11])       #iBz])
            datarray[tmp_ctr,9] = float(line[6])        #iEx])
            datarray[tmp_ctr,10] = float(line[7])       #iEy])
            datarray[tmp_ctr,11] = float(line[8])       #iEz])
            datarray[tmp_ctr,12] = float(line[12])      #itag])
#
            tmp_ctr = tmp_ctr + 1
        
        return datarray,tmp_ctr
###############################################################################
#
#
###############################################################################
    def CountParticles_V2(self, ind_npart):
        """Function that counts the total number of particles in the whole domain at a given time
        --> 1st argument = total number of particles: is used to pre-allocate the arrays to save computation time (optimization purpose).
        --> 2nd argument = 2D array of the processors (1st column) containing a given number of particles (2nd column) (optimization purpose)."""
#
        n_part = 0
#
        self.filename = self.workdir + "PIC_output" + str(self.time) + ".dat"
        dfile = open(self.filename, "r")
        tmp_head = dfile.readline().rstrip('\n\r').split()
#
        n_part = n_part + int(tmp_head[ind_npart])
#            
        return n_part
###############################################################################
#
#
###############################################################################
    def ExtractAllData_V2(self, columns, loc_PICdt, loc_nb_part, loc_nb_wrong_dt, loc_tags, info=False):
    #def ExtractAllData_V2(self, columns, loc_PICdt, loc_nb_part, loc_tags, info=False):
        """Version 2 of the function that extracts all the data about all the particles without any selection criterion/criteria.
        The last argument is the list of columns to extract.
        --> WARNING 1: you must be really careful with the order of the indices you are giving because the output array will have
            the values in the same order!!!
        --> WARNING 2: must be used alongside CountParticles_V2 in order to have the list of processors containing particles.
        --> WARNING 3: the last two arguments must always be given and must always be consistent with what was given in column. """
#
        if(loc_nb_part > 0):
            Trsf_array = np.zeros([1,np.size(columns)])
#
            self.filename = self.workdir + "PIC_output"+ str(self.time) + ".dat"
            Test = np.loadtxt(self.filename, skiprows=1, usecols=columns)
#
            try :
                Trsf_array = np.append(Trsf_array, Test, axis=0)
            except :
                Trsf_array = np.append(Trsf_array, np.array([Test]), axis=0)
#
            Trsf_array = np.delete(Trsf_array, 0, 0)
#
            nb_raws = int(np.max(Trsf_array[:,loc_tags])+1)
            nb_columns = int(np.size(columns))
#             print nb_raws
            Data_array = np.zeros([nb_raws, nb_columns])
            array_indices = (Trsf_array[:,loc_tags]).astype(int)
            Data_array[array_indices,:] = Trsf_array
#
            mask_wrong_dt = np.where(Data_array[:,loc_nb_wrong_dt] != 0.)
            Data_array[mask_wrong_dt,:] = 0.
#
            nb_particles_check = np.size(np.where(Data_array[:,loc_tags] != 0.))
            nb_particles_fucked_up = np.size(mask_wrong_dt)
#
            if(info == True):
                print "Time: " + str("%.2f" % (self.time*loc_PICdt)) + " ns, nb_particles = " + str(nb_particles_check) + "/" + str(loc_nb_part-nb_particles_fucked_up) + " with " + str(nb_particles_fucked_up) + " fucked up particles removed."
#
        else :
            Data_array = np.zeros([2, np.size(columns)])
#
        return Data_array #Trsf_array
###############################################################################
#
#
###############################################################################
    def DistinguishEnergy(self, datarray, dat_eng, erg_criterion):
        """Function reading an array of particles and selecting/distinguishing some of them with a criterion on their energy. It creates one array for each
           group of particles.
           --> 1st argument: array of particles.
           --> 2nd argument: array of particles' energy ==> Lines must match between both arrays!
           --> 3rd argument: value of the energy used for the selection, energies lower than the criterion are returned as first argument, energies larger 
               are returned as second argument. """
#
        tag_mask = np.where(dat_eng[:,0] - erg_criterion > 0.)[0]
        dataselected = np.array([datarray[i,:] for i in tag_mask])
#
        ante_tag_mask = np.where(dat_eng[:,0] - erg_criterion < 0.)[0]
        dataexcluded = np.array([datarray[j,:] for j in ante_tag_mask])
#
        return dataexcluded, dataselected
###############################################################################
#
#
###############################################################################
    def DistinguishTags(self, data, i_tag):
        """Function """
    #    
        mask_2Ei = np.where(np.mod(data[:,i_tag],2) == 0)[0]
        mask_1Ei = np.where(np.mod(data[:,i_tag],2) == 1)[0]
        
        Trsf_2Ei = np.array([data[l,:] for l in mask_2Ei])
        Trsf_1Ei = np.array([data[l,:] for l in mask_1Ei])
        
        return Trsf_1Ei, Trsf_2Ei
###############################################################################
#
#
###############################################################################
    def ReorganizeDistinguishedData(self, data, i_PICdt, i_PICtime, path):
        """Function that rewrites all of the PIC output files created by GORGON in order to make 1file/time for all the processors -
        the idea the reduction of the work that comes after by having only one file that can be loaded globaly using numpy.loadtxt -
        the header must be read by other means at the moment and there is still need to work on fast ways to distinguish the particles."""
    #
        tmp_timing_start = time.time()
    #
        """Once all of this is done you can save the global files """
        self.filename = self.workdir + "PIC_output" + str(self.time) + ".dat"
        dfile = open(self.filename, "r")
        f_head = dfile.readline().rstrip('\n\r').split()
    
        file_header = "#position_x position_y position_z impulsion_x impulsion_y impulsion_z E_x E_y E_z B_x B_y B_z tag MHD_density MHD_Ti MHD_Te MHD_Vx MHD_Vy MHD_Vz " + str(np.size(data[:,0])) + " " + str(f_head[i_PICdt]) + " " + str(f_head[i_PICtime])
    #
        """Saving the particles in a new global file"""
        new_filename = self.savepath + path + "PIC_output" + str(self.time) + ".dat"
        np.savetxt(new_filename, data, delimiter=' ', header=file_header)
    #
        tmp_timing_end = time.time()
        print("File:", new_filename, "written, nb of particles:", np.size(data[:,0]), "time taken:", tmp_timing_end-tmp_timing_start)
        print(" ")
    #
        return
###############################################################################
#
#
###############################################################################
    def PICSelection_V2(self, names_to_search, columns):
        """Function that extracts a subset of particles from their tag given in an array (2nd argument).
            --> Version 2 is a little bit slower than V1 that is in PICPostTreatment_V2, but it is more flexible considering the user can choose which
                quantities they want to create: the whole line from the file or just the positions for example, or else....
                I considered this as more powerful for the PICclass in the end."""
#
        data_to_get = np.zeros([np.size(names_to_search),np.size(columns)])
        data_to_get[:,-1] = names_to_search[:]
#
        self.filename = self.workdir + "PIC_output" + str(self.time) + ".dat"
        Test = np.loadtxt(self.filename, skiprows=1, usecols=columns)
#        
        for l in range(0,np.size(names_to_search)) :
            try :
                index = np.where(Test[:,-1] == names_to_search[l])[0]
            except :
                index = np.where(Test[-1] == names_to_search[l])[0]            #Let's say there is only 1part remaining in the whole (extreme case - for some reason)
#
            if(np.size(index) != 0) :                                          #If one of the particle the script is looking for has actually left the domain
                try:
                    data_to_get[l,:] = Test[index,:]
                except:
                    data_to_get[l,:] = Test[:]
#
        return data_to_get
###############################################################################
#

#
###############################################################################
    def ComputeEnergiesSelection(self, datmomenta, datmass, datags):
        """Function that computes the energies of the particles at a given time.
        --> WARNING: works with the momenta!"""
#
        ergs = np.zeros([np.size(datmomenta[:,0]),2])
        momenta = datmomenta[:,0]**2 + datmomenta[:,1]**2 + datmomenta[:,2]**2

        ergs[:,0] = 0.5*(datmass[:]*(c_light**2.0)*momenta[:])/(elementary_charge*((datmass[:]*c_light)**2 + momenta[:]))
        ergs[(datmass == 0.),0] = 0.
        ergs[:,1] = datags
        
      ###  print ergs[(ergs[:,0] > 8.0e6),:]
        
#                    mask_nan = np.where(datags != 0.)       #C> - handles rows where there are no particles.
#                    ergs = np.zeros([np.size(mask_nan),2])
#                    momenta = datmomenta[mask_nan,0]**2 + datmomenta[mask_nan,1]**2 + datmomenta[mask_nan,2]**2
#            #
#                #gamma = (np.sqrt(momenta[:] + (atomic_weight*proton_mass*c_light)**2))/(atomic_weight*proton_mass*c_light)
#                #ergs = (gamma[:] - 1.)*(proton_mass*atomic_weight)*c_light**2
#                #ergs = ergs/elementary_charge
#            #
#            ####Formula: E = (mc^2P^2)/{2e[(mc)^2 + P^2]} [eV]   (another formula was: E=(Gamma - 1.)*m*c^2.... relativistic case)
#                    ergs[:,0] = 0.5*(datmass[mask_nan]*(c_light**2.0)*momenta[:])/(elementary_charge*((datmass[mask_nan]*c_light)**2 + momenta[:]))
#                    ergs[:,1] = datags[mask_nan]
#
        return ergs
###############################################################################
#


#
###############################################################################
    def ComputeEnergies(self, datmomenta, datmass, datags):
        """Function that computes the energies of the particles at a given time.
        --> WARNING: works with the momenta!"""
#
        ergs = np.zeros([np.size(datags),2])
        momenta = datmomenta[:,0]**2 + datmomenta[:,1]**2 + datmomenta[:,2]**2
        
        ergs[:, 0] = 0.5*(datmass[:]*(c_light**2.0)*momenta[:])/(elementary_charge*((datmass[:]*c_light)**2 + momenta[:]))
        ergs[(datmass == 0.), 0] = 0.
        ergs[:, 1] = datags[:]

#                            mask_nan = np.where(datags != 0.)       #C> - handles rows where there are no particles.
#                            ergs = np.zeros([np.size(mask_nan),2])
#                            momenta = datmomenta[mask_nan,0]**2 + datmomenta[mask_nan,1]**2 + datmomenta[mask_nan,2]**2
#                    #
#                        #gamma = (np.sqrt(momenta[:] + (atomic_weight*proton_mass*c_light)**2))/(atomic_weight*proton_mass*c_light)
#                        #ergs = (gamma[:] - 1.)*(proton_mass*atomic_weight)*c_light**2
#                        #ergs = ergs/elementary_charge
#                    #
#                    ####Formula: E = (mc^2P^2)/{2e[(mc)^2 + P^2]} [eV]   (another formula was: E=(Gamma - 1.)*m*c^2.... relativistic case)
#                            ergs[:,0] = 0.5*(datmass[mask_nan]*(c_light**2.0)*momenta[:])/(elementary_charge*((datmass[mask_nan]*c_light)**2 + momenta[:]))
#                            ergs[:,1] = datags[mask_nan]
#
        return ergs
###############################################################################
#
#
###############################################################################
    def ConversionPj_to_Vj(self, dat_P):
        """Function that converts the momenta components into velocity components.
        --> Details: that version works globally over a larger than one (>1) number of particles and a larger than one (>1) number of times.
        --> for the version working in sequential mode over the times see git commit date: 6/6/18"""
#
        gamma_prime = (dat_P[:,0,:]**2 + dat_P[:,1,:]**2 + dat_P[:,2,:]**2)/(self.atomic_weight*proton_mass*c_light)**2
        #gamma_prime = (dat_P[:,0,:]**2 + dat_P[:,1,:]**2 + dat_P[:,2,:]**2)/(atomic_weight*proton_mass*c_light)**2
        gamma_prime = np.sqrt(1. + gamma_prime)
#
        dat_P[:,0,:] = dat_P[:,0,:]/(self.atomic_weight*proton_mass*gamma_prime)
        dat_P[:,1,:] = dat_P[:,1,:]/(self.atomic_weight*proton_mass*gamma_prime)
        dat_P[:,2,:] = dat_P[:,2,:]/(self.atomic_weight*proton_mass*gamma_prime)

        return dat_P
###############################################################################
#
#
###############################################################################
    def V1_wrt_V2(self, Vector1, Vector2, Vtags):
        """Function that computes the parallel and perpendicular components of 'Vector1' with respect to 'Vector2' """
#
        tmp_Vector2 = np.zeros([np.size(Vector2[:,0]),np.size(Vector2[0,:])])
        tmp_Vector2[:,0] = Vector2[:,0]
        tmp_Vector2[:,1] = Vector2[:,1]
        tmp_Vector2[:,2] = Vector2[:,2]
#
        V1_para = np.zeros([np.size(Vector1[:,0]),5])
        V1_perp = np.zeros([np.size(Vector1[:,0]),4])
#
        V2_mag = np.sqrt(Vector2[:,0]**2 + Vector2[:,1]**2 + Vector2[:,2]**2)
        tmp_Vector2[:,0] = Vector2[:,0]/V2_mag
        tmp_Vector2[:,1] = Vector2[:,1]/V2_mag
        tmp_Vector2[:,2] = Vector2[:,2]/V2_mag
#
        V1_para[:,0] = (Vector1[:,0]*tmp_Vector2[:,0] + Vector1[:,1]*tmp_Vector2[:,1] + Vector1[:,2]*tmp_Vector2[:,2])*tmp_Vector2[:,0]
        V1_para[:,1] = (Vector1[:,0]*tmp_Vector2[:,0] + Vector1[:,1]*tmp_Vector2[:,1] + Vector1[:,2]*tmp_Vector2[:,2])*tmp_Vector2[:,1]
        V1_para[:,2] = (Vector1[:,0]*tmp_Vector2[:,0] + Vector1[:,1]*tmp_Vector2[:,1] + Vector1[:,2]*tmp_Vector2[:,2])*tmp_Vector2[:,2]
        V1_para[:,3] = (V1_para[:,0]*tmp_Vector2[:,0] + V1_para[:,1]*tmp_Vector2[:,1] + V1_para[:,2]*tmp_Vector2[:,2])/np.sqrt(V1_para[:,0]**2 + V1_para[:,1]**2 + V1_para[:,2]**2)      #Sign with respect to the local B-field
#
        V1_perp[:,0] = Vector1[:,0] - V1_para[:,0]
        V1_perp[:,1] = Vector1[:,1] - V1_para[:,1]
        V1_perp[:,2] = Vector1[:,2] - V1_para[:,2]
#
        V1_para[:,4] = Vtags
        V1_perp[:,3] = Vtags
#
        return V1_perp, V1_para
###############################################################################
#
#
###############################################################################
    def EpeEpa(self, dat_Vperp, dat_Vpara, dat_mass):
        """Function that computes the perpendicular and parallel repartition of the energy.
        It should take as argument the perpendicular and parallel (in that order) components given by function V1_wrt_V2, which the energy
        repartition computed is perp. and para. to Vector2. In practice you should generally take the magnetic field as V2..."""
        erg_para = np.zeros([np.size(dat_Vperp[:,0])])
        erg_perp = np.zeros([np.size(dat_Vperp[:,0])])
#
        perp_mod = dat_Vperp[:,0]**2 + dat_Vperp[:,1]**2 + dat_Vperp[:,2]**2
        para_mod = dat_Vpara[:,0]**2 + dat_Vpara[:,1]**2 + dat_Vpara[:,2]**2
#
        momenta_mod = perp_mod + para_mod
#
        erg_perp[:] = 0.5*perp_mod[:]/(elementary_charge*dat_mass[:]*(1. + momenta_mod[:]/(dat_mass[:]*c_light)**2))
        erg_para[:] = 0.5*para_mod[:]/(elementary_charge*dat_mass[:]*(1. + momenta_mod[:]/(dat_mass[:]*c_light)**2))
#
        return erg_perp, erg_para
###############################################################################
#
#
###############################################################################
    def ComputeTau_ci(self, dat_B, dat_mass, dat_charge):
        """Function that computes the ion gyroperiod and converts it into nanosecond."""
#
        tmp_gyroperiod = (dat_charge[:,:]*dat_B[:,:])/dat_mass[:,:]                 #cgs fromula: ion gyrofrequency [rad/s]
        tmp_gyroperiod = 1.0e9*2.0*np.pi/tmp_gyroperiod                             #Derivation of the gyroperiod and conversion into nanosecond
#
        return tmp_gyroperiod
###############################################################################
#
#
###############################################################################
    def ComputeLarmorRadius(self, dat_Ti, dat_B, dat_mass, dat_charge):
        """Function that computes the ion gyrodaius and converts it into µ-meter."""
#
        tmp_Larmor_radii = (1.02e2*np.sqrt((dat_mass[:,:]/proton_mass)*dat_Ti[:,:])/((dat_charge[:,:])*dat_B[:,:]*1e4))               #cgs formula: ion gyroradius
        tmp_Larmor_radii = tmp_Larmor_radii*1e4                                                                     #Conversion into µ-meter
#
        return tmp_Larmor_radii
###############################################################################
#
#
###############################################################################
    def ComputeFirstAdiabaticInvariant(self, dat_Vperp, dat_B, dat_mass):
        """Function that the value of the first adiabatic invariant."""
#
        tmp_adiabatic_invariant = (dat_mass[:,:]*dat_Vperp[:,:])/(2.*dat_B[:,:])       #First adiabatic invariant in ... units
#
        return tmp_adiabatic_invariant
###############################################################################
#
#
###############################################################################





class PICPlot:
    
    def __init__(self, wrkdir, nkname, svpath, wdth, hght, tp, btm, lft, rgt):
        self.workdir = wrkdir
        self.nickname = nkname
        self.savepath = svpath
        self.width = wdth
        self.height = hght
        self.top = tp
        self.bottom = btm
        self.left = lft
        self.right = rgt
        self.filename=None
    
    
###############################################################################
    def PlotTrajectories(self, t_start, t_end, t_freq, Xlim, Ylim, dat_pos, freq_b_points, file_to_save_to, plt_specification):
        """Function """
#
        tmp_radius = np.zeros([np.size(dat_pos[:,0,0]),np.size(dat_pos[0,0,:])])
        tmp_radius[:,1:] = np.sqrt((dat_pos[:,0,1:]-Xlim)**2 + (dat_pos[:,1,1:]-Ylim)**2)
#
        X_axis = np.arange(t_start,t_end)*t_freq
#
        plt.figure(figsize=(self.width, self.height))
        plt.subplots_adjust(top=self.top, bottom=self.bottom, left=self.left, right=self.right)
#
        for i in range(0,np.size(dat_pos[:,0,0])):
            gs = gridspec.GridSpec(nrows=3, ncols=4, wspace=0.4, hspace=0.4)
            plt.subplot(gs[0, 2:])
            #ax1 = plt.xlim([-Xlim,Xlim])
            #ax1 = plt.ylim([0.,Zlim])
            plt.xlabel('x-coordinate (mm)', fontsize=22)
            plt.ylabel('z-coordinate (mm)', fontsize=22)
            plt.plot(dat_pos[i,0,np.where(dat_pos[i,0,1:] != 0.)[0]+1]-Xlim, dat_pos[i,2,np.where(dat_pos[i,2,1:] != 0.)[0]+1], plt_specification)
            plt.plot(dat_pos[i,0,1]-Xlim, dat_pos[i,2,1], 'k*')
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(dat_pos[i,0,j+freq_b_points+1]-Xlim, dat_pos[i,2,j+freq_b_points+1], 'go')
                except:
                    pass
#
            plt.subplot(gs[1, 2:])
            #ax2 = plt.xlim([-8.,8.])
            #ax2 = plt.ylim([0.,Zlim])
            plt.xlabel('radius (mm)', fontsize=22)
            plt.ylabel('z-coordinate (mm)', fontsize=22)
            plt.plot(tmp_radius[i,np.where(dat_pos[i,2,1:] != 0.)[0]+1], dat_pos[i,2,np.where(dat_pos[i,2,1:] != 0.)[0]+1], plt_specification)
            plt.plot(tmp_radius[i,1], dat_pos[i,2,1], 'k*')
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(tmp_radius[i,j+freq_b_points+1], dat_pos[i,2,j+freq_b_points+1], 'go')
                except:
                    pass
#
            plt.subplot(gs[2, 2:])
            #ax3 = plt.xlim([-Xlim,Xlim])
            #ax3 = plt.ylim([-Ylim,Ylim])
            plt.xlabel('x-coordinate (mm)', fontsize=22)
            plt.ylabel('y-coordinate (mm)', fontsize=22)
            plt.plot(dat_pos[i,0,np.where(dat_pos[i,0,1:] != 0.)[0]+1]-Xlim, dat_pos[i,1,np.where(dat_pos[i,1,1:] != 0.)[0]+1]-Ylim, plt_specification)
            plt.plot(dat_pos[i,0,1]-Xlim, dat_pos[i,1,1]-Ylim, 'k*')
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(dat_pos[i,0,j+freq_b_points+1]-Xlim, dat_pos[i,1,j+freq_b_points+1]-Ylim, 'go')
                except:
                    pass
#
            plt.subplot(gs[0, :2])
            plt.xlim(X_axis[0], X_axis[-1])
            plt.ylabel('x-coordinate (mm)', fontsize=22)
            plt.plot(X_axis, dat_pos[i,0,1:]-Xlim, plt_specification)
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], dat_pos[i,0,j+freq_b_points+1]-Xlim, 'go')
                except:
                    pass
#
            plt.subplot(gs[1, :2])
            plt.xlim(X_axis[0], X_axis[-1])
            plt.ylabel('y-coordinate (mm)', fontsize=22)
            plt.plot(X_axis, dat_pos[i,1,1:]-Ylim, plt_specification)
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], dat_pos[i,1,j+freq_b_points+1]-Ylim, 'go')
                except:
                    pass
#
            plt.subplot(gs[2, :2])
            plt.xlim(X_axis[0], X_axis[-1])
            plt.ylabel('z-coordinate (mm)', fontsize=22)
            plt.xlabel('Time (ns)', fontsize=22)
            plt.plot(X_axis, dat_pos[i,2,1:], plt_specification)
            for j in range(0,np.size(dat_pos[i,0,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], dat_pos[i,2,j+freq_b_points+1], 'go')
                except:
                    pass
#   
            title = 'Particle: ' + str(int(dat_pos[i,0,0]*1e-3))
            plt.suptitle(title, fontsize=20)
            #plt.tight_layout()
#
            title = self.savepath + file_to_save_to + self.nickname + "_space_" + str(int(dat_pos[i,0,0]*1e-3))
            plt.savefig(title, format='png', frameon=None)
    
            if(i % 100 == 0):
                print("Plot_trajectories - ", i, " images plotted")
#
            plt.clf()
        plt.close()
        return
###############################################################################
#
#
###############################################################################
    def PlotPlasmaQuantities(self, t_start, t_end, t_freq, dat_erg, Enormalization, data, freq_b_points, file_to_save_to, plt_specification):
        """Function  """
#
        X_axis = np.arange(t_start,t_end)*t_freq
#
        plt.figure(figsize=(self.width, self.height))
        plt.subplots_adjust(top=self.top, bottom=self.bottom, left=self.left, right=self.right)
#
        for i in range(0,np.size(dat_erg[:,0])):
            plt.subplot(414)
            plt.plot(X_axis, data[i,0,:])
            for j in range(0,np.size(data[i,0,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], data[i,0,j+freq_b_points], 'go', markersize=8)
                except:
                    pass
            plt.xlim(X_axis[0], X_axis[-1])
            plt.ylabel(r"$\tau_{ci} \, \left( \, ns \, \right)$", fontsize=30)
            plt.xlabel(r"$Time \, \left( \, ns \, \right)$", fontsize=30)
#
            plt.subplot(413)
            plt.plot(X_axis,data[i,1,:])
            for j in range(0,np.size(data[i,1,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], data[i,1,j+freq_b_points], 'go', markersize=8)
                except:
                    pass
            plt.xlim(X_axis[0], X_axis[-1])
            plt.xticks([])
            plt.ylabel(r"$\rho_{ci} \, \left( \, \mu m \, \right)$", fontsize=30)
#
            plt.subplot(412)
            plt.plot(X_axis, data[i,2,:]/data[i,2,0])
            for j in range(0,np.size(data[i,2,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], data[i,2,j+freq_b_points]/data[i,2,0], 'go', markersize=8)
                except:
                    pass
            plt.xlim(X_axis[0], X_axis[-1])
            plt.xticks([])
            plt.ylabel(r'$\mu / \mu_0$', fontsize=30)
#
            plt.subplot(411)#, sharex=ax3)
            plt.plot(X_axis, dat_erg[i,1:]/Enormalization)
            for j in range(0,np.size(dat_erg[i,11:]),freq_b_points):
                try:
                    plt.plot(X_axis[j+freq_b_points], dat_erg[i,j+freq_b_points+1]/Enormalization, 'go', markersize=8)
                except:
                    pass
            plt.xticks([])
            plt.ylim([1.0e3, 2.0e7])
            plt.yscale('log')
            plt.ylabel(r'$\varepsilon \, / \,\varepsilon_{ref}$', fontsize=30)
#
            plt.xlim(X_axis[0], X_axis[-1])
            title = 'Particle: ' + str(int(dat_erg[i,0]))
            plt.suptitle(title, fontsize=18)
            plt.tight_layout()
#
            title = self.savepath + file_to_save_to + self.nickname + "_general_" + str(int(dat_erg[i,0]))
            plt.savefig(title, format='png', frameon=None)
#
            if(i % 100 == 0):
                print("PlotPlasmaQuantities - ", i, " images plotted")
#
            plt.clf()
        plt.close()
#
        average = np.zeros([np.size(data[:,0,0]),6])
        for j in range(0,np.size(data[:,0,0])):
            average[j,0] = np.nanmean(data[j,1,:])
            average[j,1] = np.nanstd(data[j,1,:])
            average[j,2] = np.nanmean(data[j,2,:]/data[j,2,0])
            average[j,3] = np.nanstd(data[j,2,:]/data[j,2,0])
            average[j,4] = np.mean(dat_erg[j,1:]/Enormalization)
            average[j,5] = np.std(dat_erg[j,1:]/Enormalization)
#    
        print("Energy: average, standard deviation and std/average - ", np.mean(average[:,4]), np.mean(average[:,5]), " and ", np.mean(average[:,5])/np.mean(average[:,4]))
        print("Larmor radius: average, std and std/average - ", np.mean(average[:,0]), np.mean(average[:,1]), "and ", np.mean(average[:,1])/np.mean(average[:,0]))
        print("First adiabatic invariant: avg, std and std/avg - ", np.average(average[:,2]), np.average(average[:,3]), "and ", np.average(average[:,3])/np.average(average[:,2]))
#
        return
###############################################################################
#
#
###############################################################################
#            def Plot_over_MHD_density(Positions, Energies, Eref, name_list, t_start, t_end, t_out_freq, MHD_dir, MHD_pref, nyt):
#                plt.rc('figure.subplot', top=0.85)
#                plt.rc('figure.subplot', bottom=0.25)
#                plt.rc('figure.subplot', left=0.10)
#                plt.rc('figure.subplot', right=0.90)
#                
#                plt.figure(figsize=(12,12))
#                gs = gridspec.GridSpec(2, 2)
#            
#                for i in range(t_start, t_end):
#                    a = PVTI.PVTI(MHD_dir,MHD_pref,i)
#                    print("MHD time:", i, " calling for PVTIclass")
#                    a.GetData()
#                    dl = a.CellSize()
#                
#                    for j in name_list:
#                        part = np.where(Positions[:,0,0] == j)[0]
#                        loc = int(Positions[part,2,i+1]/dl)
#                        tag = Positions[part,0,0]
#                        height = round(Positions[part,2,i+1]*1e3, 2)
#                        markersize = 5.0 #+ np.log((PICData_sub1[part,0,i+1]/Einit)[0])
#                
#                        plt.subplot(gs[0, 1])
#                        a.slice(a.d,loc,'Z',linlog=True,low=-4.,up=-1.,interp='bicubic',scale=1e3,colbar=False,clear=False,colmap=plt.cm.YlOrRd)
#                        acolbay = plt.colorbar()
#                        acolbay.set_label('log' r'$_{10} \, \, \rho$' ' (kg/m' r'$^{3}$' ')', size=20, color="black")
#                
#                        plt.plot(Positions[part,1,i+1]*1e3-7.5, Positions[part,0,i+1]*1e3-7.5, 'bo', markersize=markersize)
#                        #plt.text(-7.02, 7.0,str(i*PIC_time_out_frequency) + r'$\,ns$',size=20,color="black",horizontalalignment='left',verticalalignment='top')
#                        #plt.text(-7.20, 7.0,r"$z =$" + str(height) + r'$mm$',size=20,color="black",horizontalalignment='left',verticalalignment='top')
#                
#                        plt.xlim([-7.5,7.5])
#                        plt.ylim([-7.5,7.5])
#                        plt.xticks([-6,-3,0,3,6]) #-7,-5,-3,-1,1,3,5,7])
#                        plt.yticks([-6,-3,0,3,6]) #-7,-5,-3,-1,1,3,5,7])
#                        plt.xlabel("y-axis (mm)")
#                        plt.ylabel("x-axis (mm)")
#                
#                        plt.subplots_adjust(wspace=0.4)
#            
#                        plt.subplot(gs[0, 0])
#                        a.slice(a.d,int(nyt/2),'Y',linlog=True,low=-4.,up=-1.,interp='bicubic',scale=1e3,colbar=False,clear=False,colmap=plt.cm.YlOrRd)
#                        plt.plot([-7.5,7.5], [loc*dl*1e3,loc*dl*1e3], 'k-')
#                        if(Positions[part,1,i+1]*1e3-7.5 < 0.0) :
#                            plt.plot(Positions[part,0,i+1]*1e3-7.5, Positions[part,2,i+1]*1e3, 'wo')
#                        else :
#                            plt.plot(Positions[part,0,i+1]*1e3-7.5, Positions[part,2,i+1]*1e3, 'bo')
#                        
#                        plt.text(-7.02, 11.5,str(i*t_out_freq) + r'$\,ns$',size=20,color="black",horizontalalignment='left',verticalalignment='top')
#                        plt.text(1.70, 11.5,r"$z =$" + str(height) + r'$mm$',size=20,color="black",horizontalalignment='left',verticalalignment='top')
#                        
#                        plt.xlim([-7.5,7.5])
#                        plt.ylim([0.,12.])
#                        plt.xticks([-6,-3,0,3,6])
#                        plt.yticks([0,4,8,12])
#                        plt.xlabel(r"x-axis (mm)")
#                        plt.ylabel(r"z-axis (mm)")
#                        
#                        plt.subplots_adjust(hspace=0.4)
#            
#                        X_axis = np.arange(t_start,t_end)*t_out_freq
#                        plt.subplot(gs[1,:])
#                        plt.plot(X_axis[:i+1], Energies[part,0,1:i+2][0,:]/Eref, label=r'total')
#                        plt.plot(X_axis[i], Energies[part,0,i+1]/Eref, 'bo')
#                        plt.plot(X_axis[:i+1], Energies[part,1,1:i+2][0,:]/Eref, 'r-', label=r'$\perp$-part')
#                        plt.plot(X_axis[:i+1], Energies[part,2,1:i+2][0,:]/Eref, 'g-', label=r'$\parallel$-part')
#                        plt.legend(loc='upper left', fontsize=20)
#                        
#                        plt.xlim([X_axis[0],X_axis[-1]])
#                        plt.ylim([0.,np.max(Energies[part,0,1:])/Eref])
#                        plt.xlabel(r"Time (ns)")
#                        plt.ylabel(r"$Energy/E_{ref}$")
#                
#                        supertitle = r"$\rho_{MHD}$ with particle " + str(int(tag[0]))
#                        plt.suptitle(supertitle, y=0.95, fontsize=25)
#                
#                        SecurityCheckPaths(str(int(tag[0])) + "/")
#                        if i < 10 :
#                            title = savepath + str(int(tag[0])) + "/" + "PIC_over_MHD_00" + str(i)
#                        elif i >=10 and i < 100 :
#                            title = savepath + str(int(tag[0])) + "/" + "PIC_over_MHD_0" + str(i)
#                        else:
#                            title = savepath + str(int(tag[0])) + "/" + "PIC_over_MHD_" + str(i)
#                        plt.savefig(title, format='png', frameon=None)
#                        plt.clf()
#            
#                print("Particles plotted over MHD density: time " + str(i))
###############################################################################
#
#
###############################################################################







###############################################################################
#
#
###############################################################################



class PICreorganize:
    """This version of the PICclass works on the original data files produced by the MHD code GORGON (more specifically on the PIC files produced by it).
       The main difference with 'PICclass' is in the way it reads the data files; that is for one time it has to read as many files as the number of
       processors used in the simulation.
       Otherwise all the computation and transformation functions are the same... (dated 2018-08-22 10:24:35) """
    
    def __init__(self, wrkdir, pfx, tm, procs, svpath, chrge, atm_mass, fname=False):
        self.workdir = wrkdir
        self.prefix = pfx
        self.time = tm
        self.nb_proc = procs
        self.savepath = svpath
        self.charge = chrge
        self.atomic_weight = atm_mass
        self.filename=fname


###############################################################################
    def SecurityCheckPaths_OriginalFiles(self, dir_to_check):
        """Function checking the existence of the directories where the data are said to be stored and where the figures
        are said to be saved, that is the main directory and all the subdirectories of the tree according to the flags.
        -
        WARNING: For now it works without any arguments because it is all in the same script, but when going to a
                 modular way of working, it will be needed to add the arguments. """
        #
        if not os.path.exists(self.workdir):
            sys.exit("PICclass: directory to get the data does not exist... Stopping execution.")
        #
        if not os.path.exists(self.savepath + dir_to_check):
            os.makedirs(self.savepath + dir_to_check)
            print("SecurityCheckPaths function: creation of the directory: " + self.savepath + dir_to_check + ".")
        
        return
###############################################################################
#
#
###############################################################################
    def Header_OriginalFiles(self):
        """Gets the indices of the quantities by reading the header of a PIC data file.
        --> 3 for particle position, 3 for particle momentum components, 3 for E-field components, 3 for B-field components, 
            1 for the tag, 1 for the number of particles, 1 for PIC time step and 1 for the PIC time."""
#
        src = open(self.workdir+"PIC_output0.0.dat", "r")
        tmp_head = src.readline().rstrip('\n\r').split()
#
        idx_pos_x = tmp_head.index('#position_x')
        idx_pos_y = tmp_head.index('position_y')
        idx_pos_z = tmp_head.index('position_z')
#
        idx_imp_x = tmp_head.index('impulsion_x')
        idx_imp_y = tmp_head.index('impulsion_y')
        idx_imp_z = tmp_head.index('impulsion_z')
#
        idx_B_x = tmp_head.index('B_x')
        idx_B_y = tmp_head.index('B_y')
        idx_B_z = tmp_head.index('B_z')
#
        idx_E_x = tmp_head.index('E_x')
        idx_E_y = tmp_head.index('E_y')
        idx_E_z = tmp_head.index('E_z')
#
        idx_mass = tmp_head.index('mass')
        idx_charge = tmp_head.index('charge')
        idx_tag = tmp_head.index('tag')
#
        idx_MHD_rho = tmp_head.index('MHD_density')
        idx_MHD_Ti = tmp_head.index('MHD_Ti')
        idx_MHD_Te = tmp_head.index('MHD_Te')
        idx_MHD_Vx = tmp_head.index('MHD_Vx')
        idx_MHD_Vy = tmp_head.index('MHD_Vy')
        idx_MHD_Vz = tmp_head.index('MHD_Vz')
#
        idx_nb_part = idx_MHD_Vz + 2
        idx_PIC_dt = idx_MHD_Vz + 3
        idx_PIC_time = idx_MHD_Vz + 4

        return idx_pos_x,idx_pos_y,idx_pos_z,idx_imp_x,idx_imp_y,idx_imp_z,idx_B_x,idx_B_y,idx_B_z,idx_E_x,idx_E_y,idx_E_z,idx_mass,idx_charge,idx_tag, idx_MHD_rho,idx_MHD_Ti,idx_MHD_Te,idx_MHD_Vx,idx_MHD_Vy,idx_MHD_Vz, idx_nb_part,idx_PIC_dt,idx_PIC_time
###############################################################################
#
#
###############################################################################
    def CountParticles_OriginalFiles(self, ind_npart):
        """Function that counts the total number of particles in the whole domain at a given time
        --> 1st argument = total number of particles: is used to pre-allocate the arrays to save computation time (optimization purpose).
        --> 2nd argument = 2D array of the processors (1st column) containing a given number of particles (2nd column) (optimization purpose)."""
#
        n_proc = np.zeros([1,2])
        loc_ctr = 0
        n_part = 0
        for j in range(0,self.nb_proc):
            self.filename = self.workdir + "PIC_output" + str(j) + "." + str(self.time) + ".dat"
            dfile = open(self.filename, "r")
            tmp_head = dfile.readline().rstrip('\n\r').split()
#
            if int(tmp_head[ind_npart]) != 0 :
                n_part = n_part + int(tmp_head[ind_npart])
            
                truc = np.array([[j, int(tmp_head[ind_npart])]])           #Creates a 1D array storing the number of the processor and the number of particles it contains
            
                loc_ctr = loc_ctr + 1
                n_proc = np.append(n_proc, truc, 0)
#
        n_proc = np.delete(n_proc,0,0)
        #print self.time, n_part
            
        return n_part, np.int0(n_proc)
###############################################################################
#
#
###############################################################################
    def ExtractAllData_OriginalFiles(self, proc_list, columns):
        """Version 2 of the function that extracts all the data about all the particles without any selection criterion/criteria.
        The last argument is the list of columns to extract.
        --> WARNING 1: you must be really careful with the order of the indices you are giving because the output array will have
            the values in the same order!!!
        --> WARNING 2: must be used alongside CountParticles_V2 in order to have the list of processors containing particles."""
#
        Trsf_array = np.zeros([1,np.size(columns)])
#
        for j in proc_list[:,0] :
            self.filename = self.workdir + "PIC_output" + str(j) + "." + str(self.time) + ".dat"
            #print self.filename

            Test = np.loadtxt(self.filename, skiprows=1, usecols=columns)

            try :
                Trsf_array = np.append(Trsf_array, Test, axis=0)
            except :
                Trsf_array = np.append(Trsf_array, np.array([Test]), axis=0)
#
        Trsf_array = np.delete(Trsf_array, 0, 0)

        return Trsf_array
###############################################################################
#




class PIC_cst:
    
    def __init__(self):
        self.filename=None

    def Constant(self):
        self.elementary_charge = 1.6022e-19
        self.gamma = 5./3.
        self.K_Boltzmann = 1.3807e-23
        self.proton_mass = 1.6726e-27
        self.electron_mass = 9.1094e-31
        self.c_light = 2.9979e8/20.
        self.mu_zero = 4.0*np.pi*1e-7

    
###############################################################################








