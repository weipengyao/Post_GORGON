# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 14:40:58 2018

@author: kapichu
"""


import PICclass as PIC
import numpy as np
import time

###############################################################################
### Data ###
"""Here are all the stuff to declare before running the script. It consists in the name of the data files to be analyzed but also
in the details of the simulation box as well as how you want to analyze them """
#
##################################
###   Name of the simulation   ###
##################################
prefix = "PA7" # "PA_3_6"                                                              #"BWtest3"              #"PIC_results"

print("Studying simulation", prefix)
print("   ---------------------------------------------   ")


##################################
###    Path to get the data    ###
##################################
workdir = "/home/kapichu/mesopsl_work/Results/" + prefix + "/PIC_out_domain/"
##################################


##################################
###   Path to save the files   ###
##################################
#savepath = "/home/kapichu/Bureau/" + prefix + "/PIC_out_domain/"
#savepath = "/home/kapichu/Bureau/PIC/" + prefix + "/PIC_out_domain/"
#savepath = "/media/kapichu/data/Work_D2/Jets_collision/" + prefix + "/PIC_out_domain/"
savepath = "/media/kapichu/data/Work_D2/Weibul/" + prefix + "/PIC_out_domain/"


##################################
###       Simulation info      ###
##################################
nb_processor = 75#384
nb_particles = 0


#######################################
###   PICclass required arguments   ###
#######################################
Z = 1.                                                                         #Number of charge
atomic_weight = 1.                                                             #Ion mass expressed in units of the proton mass


###############################################################################
##           Function that rewrites all the original PIC files               ##
###############################################################################
pic = PIC.PICreorganize(workdir, prefix, 0, nb_processor, savepath, Z, atomic_weight)
pic.SecurityCheckPaths_OriginalFiles("")

print("Starting to rewrite the files - time:", time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))
timing_start = time.time()


Data = np.zeros([1,16])

for i in range(0,nb_processor):
    tmp_timing_start = time.time()
    loc_nb_particles = 0
#
    try:
        filename = workdir + "PIC_out_domain" + str(i) + ".dat"
        Test = np.loadtxt(filename, skiprows=1)
#
        try:
            Data = np.append(Data, Test, axis=0)
            loc_nb_particles = np.size(Test[:,0])
            nb_particles = nb_particles + loc_nb_particles
        except :                                                #Manages the case when there is only one particle in the file (exception...)
            Data = np.append(Data, np.array([Test]), axis=0)
            loc_nb_particles = 1
            nb_particles = nb_particles + loc_nb_particles
#
    except:
        print("No output files created by the processor:", i)
        pass
#
    tmp_timing_end = time.time()
    print("File from processor ", i, "read, nb of particles found:", loc_nb_particles, "time taken:", tmp_timing_end-tmp_timing_start)
    print(" ")
#
Data = np.delete(Data, 0, 0)
#

   
#        
"""Once all of the files have been read you can save everything in a global files """
file_header = "#position_x position_y position_z impulsion_x impulsion_y impulsion_z E_x E_y E_z B_x B_y B_z mass charge tag time"# + str(nb_particles)
#
"""Saving the particles in the global file"""
new_filename = savepath + "PIC_out_domain.dat"
np.savetxt(new_filename, Data, delimiter=' ', header=file_header)
#


#
timing_end = time.time()
print "Ended to rewrite the files - time:", time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
print "Number of particles:", nb_particles, "| Time taken:", timing_end-timing_start, " seconds"
print "End"
###############################################################################
##                                    End                                    ##
###############################################################################























