#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 12:49:26 2021

@author: evagmelichmeijling
"""

import numpy as np
from astropy.io import fits
from astropy.table import QTable
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import DAOStarFinder, CircularAperture, CircularAnnulus
from datetime import datetime
from photutils.datasets import make_100gaussians_image
from photutils import aperture_photometry
from photutils.background import Background2D, MedianBackground
from astropy.stats import sigma_clipped_stats, SigmaClip
from scipy import stats
import glob
import csv
import pandas as pd

start = datetime.now()

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'a')
DOUBLECLUSTERTXT.write("3_MAGCALIBRATE:\n")

###########################################NGC869A###########################################

#Read 1_Starfinder CSV to import instrumental magnitude values
StarfinderCSV_NGC869A = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869A.csv')

#Make lists of the instrumental magnitudes and their errors for every filter and exposure time
magB_inst = []
magB_error = []
aperture_sum_B = []
max_aperture_B = []
SNR_B = []

magRlong_inst = []
magRlong_error = []
aperture_sum_Rlong = []
max_aperture_Rlong = []
SNR_Rlong = []

magRshort_inst = []
magRshort_error = []
aperture_sum_Rshort = []
max_aperture_Rshort = []
SNR_Rshort = []

magVlong_inst = []
magVlong_error = []
aperture_sum_Vlong = []
max_aperture_Vlong = []
SNR_Vlong = []

magVshort_inst = []
magVshort_error = []
aperture_sum_Vshort = []
max_aperture_Vshort = []
SNR_Vshort = []

id_list = []
                  
#Load RA and DEC coordinates and put into list
NGC869_A_WCSpositions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC869_A_WCSpositions.csv', delimiter = ',')

RA = []
DEC = []
for i in range(len(NGC869_A_WCSpositions)):
    RA.append(NGC869_A_WCSpositions[i, 0])
    DEC.append(NGC869_A_WCSpositions[i, 1])

#Fill empty lists with appropriate magnitude values
for i in range(len(StarfinderCSV_NGC869A['filtertype'])):
    if StarfinderCSV_NGC869A['filtertype'][i] == 'B':
        starid = StarfinderCSV_NGC869A['id'][i]
        uniqueid = f"{starid}_NGC869A"
        id_list.append(uniqueid)
        magB_inst.append(StarfinderCSV_NGC869A['mag inst'][i])
        magB_error.append(StarfinderCSV_NGC869A['mag error'][i])
        aperture_sum_B.append(StarfinderCSV_NGC869A['aperture_sum'][i])
        max_aperture_B.append(StarfinderCSV_NGC869A['max aperture'][i])
        SNR_B.append(StarfinderCSV_NGC869A['SNR'][i])
    
    if StarfinderCSV_NGC869A['filtertype'][i] == 'R':
        if StarfinderCSV_NGC869A['exposuretime'][i] == 60.0:
            magRlong_inst.append(StarfinderCSV_NGC869A['mag inst'][i])
            magRlong_error.append(StarfinderCSV_NGC869A['mag error'][i])
            aperture_sum_Rlong.append(StarfinderCSV_NGC869A['aperture_sum'][i])
            max_aperture_Rlong.append(StarfinderCSV_NGC869A['max aperture'][i])
            SNR_Rlong.append(StarfinderCSV_NGC869A['SNR'][i])
    
    if StarfinderCSV_NGC869A['filtertype'][i] == 'R':
        if StarfinderCSV_NGC869A['exposuretime'][i] == 5.0:
            magRshort_inst.append(StarfinderCSV_NGC869A['mag inst'][i])
            magRshort_error.append(StarfinderCSV_NGC869A['mag error'][i])
            aperture_sum_Rshort.append(StarfinderCSV_NGC869A['aperture_sum'][i])
            max_aperture_Rshort.append(StarfinderCSV_NGC869A['max aperture'][i])
            SNR_Rshort.append(StarfinderCSV_NGC869A['SNR'][i])
         
    if StarfinderCSV_NGC869A['filtertype'][i] == 'V':
        if StarfinderCSV_NGC869A['exposuretime'][i] == 60.0:
            magVlong_inst.append(StarfinderCSV_NGC869A['mag inst'][i])
            magVlong_error.append(StarfinderCSV_NGC869A['mag error'][i])
            aperture_sum_Vlong.append(StarfinderCSV_NGC869A['aperture_sum'][i])
            max_aperture_Vlong.append(StarfinderCSV_NGC869A['max aperture'][i])
            SNR_Vlong.append(StarfinderCSV_NGC869A['SNR'][i])
    
    if StarfinderCSV_NGC869A['filtertype'][i] == 'V':
        if StarfinderCSV_NGC869A['exposuretime'][i] == 5.0:
            magVshort_inst.append(StarfinderCSV_NGC869A['mag inst'][i])
            magVshort_error.append(StarfinderCSV_NGC869A['mag error'][i])
            aperture_sum_Vshort.append(StarfinderCSV_NGC869A['aperture_sum'][i])
            max_aperture_Vshort.append(StarfinderCSV_NGC869A['max aperture'][i])
            SNR_Vshort.append(StarfinderCSV_NGC869A['SNR'][i])

#Create a table with one row for every star containing the magnitude values and errors for the various filters
magtable = QTable([id_list, RA, DEC, 
                   magRlong_inst, magRlong_error, 
                   magRshort_inst, magRshort_error, 
                   magVlong_inst, magVlong_error, 
                   magVshort_inst, magVshort_error, 
                   magB_inst, magB_error], 
                   names=('id', 'RA', 'DEC', 
                          'magRlong_inst', 'magRlong_error',
                          'magRshort_inst', 'magRshort_error', 
                          'magVlong_inst', 'magVlong_error',
                          'magVshort_inst', 'magVshort_error', 
                          'magB_inst', 'magB_error'))    

for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(magtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869A.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = magtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(magtable))

#Dictionary with literature reference star, and their RVB values & ID
lit_refstars_RVB_id= {
    'NGC 869 850 ra34.7034675246 dec+57.2029082369' : [13.610, 13.860, 14.690, 1371],
    'NGC 869 1166 ra34.8101488390 dec+57.2093657520' : [13.060, 13.125, 13.700, 1886],
    'NGC 869 1099 ra34.7832449900 dec+57.1727922525' : [13.090, 13.260, 13.790, 1740]} #1000 is verzonnen

#Plot B-V diagram w/ instrumental magnitudes
plt.scatter(magtable['magB_inst']-magtable['magVshort_inst'], magtable['magVshort_inst'], s=1)
plt.gca().invert_yaxis()
plt.xlabel("Uncalibrated B-V", fontsize=20)
plt.ylabel("V", fontsize=20)
plt.title('Color-magnitude diagram NGC869A', fontsize=20)
plt.show()

#Lists containing calibration factor per filter (one for every reference star)
C_Rlonglist = []
C_Rshortlist = []
C_Vlonglist = []
C_Vshortlist = []
C_Blist = []

#List of names referencestars taken from dictionary
lit_refstars_names = lit_refstars_RVB_id.keys()
#2D array of magnitudes&id of referencestars
lit_refstars_values = list(lit_refstars_RVB_id.values())


for i in range(len(lit_refstars_names)):
    #List of magnitudes&id of one referencestar
    lit_refstar_values = lit_refstars_values[i]
    #Find instrumental magnitudes of referencestar in magtable using star id
    referencestar = magtable[:][lit_refstar_values[3] - 1]
    
    #Calculate calibration constant: Minst,ref - Mlit,ref
    C_Rlong = referencestar[3] - lit_refstar_values[0]
    C_Rshort = referencestar[5] - lit_refstar_values[0]
    
    C_Vlong = referencestar[7] - lit_refstar_values[1]
    C_Vshort = referencestar[9] - lit_refstar_values[1]
    
    C_B = referencestar[11] - lit_refstar_values[2]
    
    #Append C_RVB values to list
    C_Rlonglist.append(C_Rlong)
    C_Rshortlist.append(C_Rshort)
    
    C_Vlonglist.append(C_Vlong)
    C_Vshortlist.append(C_Vshort)
    
    C_Blist.append(C_B)
    

#Print table with calibrated magnitudes per filter & exposure time (3 times, one for every 
#  literature reference star)
calmagtable = QTable([id_list , np.round(RA, 3), np.round(DEC, 3),
                      
                      aperture_sum_Rlong, max_aperture_Rlong, SNR_Rlong, magRlong_inst, magRlong_error,
                      magRlong_inst - C_Rlonglist[0], 
                      magRlong_inst - C_Rlonglist[1],
                      magRlong_inst - C_Rlonglist[2],
                      
                      aperture_sum_Rshort, max_aperture_Rshort, SNR_Rshort, magRshort_inst, magRshort_error,
                      magRshort_inst - C_Rshortlist[0],
                      magRshort_inst - C_Rshortlist[1],
                      magRshort_inst - C_Rshortlist[2],
                      
                      aperture_sum_Vlong, max_aperture_Vlong, SNR_Vlong, magVlong_inst, magVlong_error,
                      magVlong_inst - C_Vlonglist[0], 
                      magVlong_inst - C_Vlonglist[1],
                      magVlong_inst - C_Vlonglist[2],
                      
                      aperture_sum_Vshort, max_aperture_Vshort, SNR_Vshort, magVshort_inst, magVshort_error,
                      magVshort_inst - C_Vshortlist[0],
                      magVshort_inst - C_Vshortlist[1],
                      magVshort_inst - C_Vshortlist[2],
                      
                      aperture_sum_B, max_aperture_B, SNR_B, magB_inst, magB_error,
                      magB_inst - C_Blist[0], 
                      magB_inst - C_Blist[1],
                      magB_inst - C_Blist[2]],
                      names=('id', 'RA', 'DEC',
                             'aperture_sum_Rlong', 'max_aperture_Rlong', 'SNR_Rlong', 'magRlong_inst', 'magRlong_error',
                             'magRlong_ref1', 'magRlong_ref2', 'magRlong_ref3',
                             
                             'aperture_sum_Rshort', 'max_aperture_Rshort', 'SNR_Rshort', 'magRshort_inst', 'magRshort_error',
                             'magRshort_ref1', 'magRshort_ref2', 'magRshort_ref3',
                             
                             'aperture_sum_Vlong', 'max_aperture_Vlong', 'SNR_Vlong', 'magVlong_inst', 'magVlong_error',
                             'magVlong_ref1', 'magVlong_ref2', 'magVlong_ref3',
                             
                             'aperture_sum_Vshort', 'max_aperture_Vshort', 'SNR_Vshort', 'magVshort_inst', 'magVshort_error',
                             'magVshort_ref1', 'magVshort_ref2', 'magVshort_ref3',
                             
                             'aperture_sum_B', 'max_aperture_B', 'SNR_B', 'magB_inst', 'magB_error',
                             'magB_ref1', 'magB_ref2', 'magB_ref3'))
                                        
for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(calmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869A_calibrated.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = calmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(calmagtable))


###########################################NGC869B###########################################

#Read 1_Starfinder CSV to import instrumental magnitude values
StarfinderCSV_NGC869B = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869B.csv')

#Make lists of the instrumental magnitudes and their errors for every filter and exposure time
magB_inst = []
magB_error = []
aperture_sum_B = []
max_aperture_B = []
SNR_B = []

magRlong_inst = []
magRlong_error = []
aperture_sum_Rlong = []
max_aperture_Rlong = []
SNR_Rlong = []

magRshort_inst = []
magRshort_error = []
aperture_sum_Rshort = []
max_aperture_Rshort = []
SNR_Rshort = []

magVlong_inst = []
magVlong_error = []
aperture_sum_Vlong = []
max_aperture_Vlong = []
SNR_Vlong = []

magVshort_inst = []
magVshort_error = []
aperture_sum_Vshort = []
max_aperture_Vshort = []
SNR_Vshort = []

id_list = []
                  
#Load RA and DEC coordinates and put into list
NGC869_B_WCSpositions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC869_B_WCSpositions.csv', delimiter = ',')

RA = []
DEC = []
for i in range(len(NGC869_B_WCSpositions)):
    RA.append(NGC869_B_WCSpositions[i, 0])
    DEC.append(NGC869_B_WCSpositions[i, 1])

#Fill empty lists with appropriate magnitude values
for i in range(len(StarfinderCSV_NGC869B['filtertype'])):
    if StarfinderCSV_NGC869B['filtertype'][i] == 'B':
        starid = StarfinderCSV_NGC869B['id'][i]
        uniqueid = f"{starid}_NGC869B"
        id_list.append(uniqueid)
        magB_inst.append(StarfinderCSV_NGC869B['mag inst'][i])
        magB_error.append(StarfinderCSV_NGC869B['mag error'][i])
        aperture_sum_B.append(StarfinderCSV_NGC869B['aperture_sum'][i])
        max_aperture_B.append(StarfinderCSV_NGC869B['max aperture'][i])
        SNR_B.append(StarfinderCSV_NGC869B['SNR'][i])
    
    if StarfinderCSV_NGC869B['filtertype'][i] == 'R':
        if StarfinderCSV_NGC869B['exposuretime'][i] == 60.0:
            magRlong_inst.append(StarfinderCSV_NGC869B['mag inst'][i])
            magRlong_error.append(StarfinderCSV_NGC869B['mag error'][i])
            aperture_sum_Rlong.append(StarfinderCSV_NGC869B['aperture_sum'][i])
            max_aperture_Rlong.append(StarfinderCSV_NGC869B['max aperture'][i])
            SNR_Rlong.append(StarfinderCSV_NGC869B['SNR'][i])
    
    if StarfinderCSV_NGC869B['filtertype'][i] == 'R':
        if StarfinderCSV_NGC869B['exposuretime'][i] == 5.0:
            magRshort_inst.append(StarfinderCSV_NGC869B['mag inst'][i])
            magRshort_error.append(StarfinderCSV_NGC869B['mag error'][i])
            aperture_sum_Rshort.append(StarfinderCSV_NGC869B['aperture_sum'][i])
            max_aperture_Rshort.append(StarfinderCSV_NGC869B['max aperture'][i])
            SNR_Rshort.append(StarfinderCSV_NGC869B['SNR'][i])
         
    if StarfinderCSV_NGC869B['filtertype'][i] == 'V':
        if StarfinderCSV_NGC869B['exposuretime'][i] == 60.0:
            magVlong_inst.append(StarfinderCSV_NGC869B['mag inst'][i])
            magVlong_error.append(StarfinderCSV_NGC869B['mag error'][i])
            aperture_sum_Vlong.append(StarfinderCSV_NGC869B['aperture_sum'][i])
            max_aperture_Vlong.append(StarfinderCSV_NGC869B['max aperture'][i])
            SNR_Vlong.append(StarfinderCSV_NGC869B['SNR'][i])
    
    if StarfinderCSV_NGC869B['filtertype'][i] == 'V':
        if StarfinderCSV_NGC869B['exposuretime'][i] == 5.0:
            magVshort_inst.append(StarfinderCSV_NGC869B['mag inst'][i])
            magVshort_error.append(StarfinderCSV_NGC869B['mag error'][i])
            aperture_sum_Vshort.append(StarfinderCSV_NGC869B['aperture_sum'][i])
            max_aperture_Vshort.append(StarfinderCSV_NGC869B['max aperture'][i])
            SNR_Vshort.append(StarfinderCSV_NGC869B['SNR'][i])

#Create a table with one row for every star containing the magnitude values and errors for the various filters
magtable = QTable([id_list, RA, DEC, 
                   magRlong_inst, magRlong_error, 
                   magRshort_inst, magRshort_error, 
                   magVlong_inst, magVlong_error, 
                   magVshort_inst, magVshort_error, 
                   magB_inst, magB_error], 
                   names=('id', 'RA', 'DEC', 
                          'magRlong_inst', 'magRlong_error',
                          'magRshort_inst', 'magRshort_error', 
                          'magVlong_inst', 'magVlong_error',
                          'magVshort_inst', 'magVshort_error', 
                          'magB_inst', 'magB_error'))    

for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(magtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869B.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = magtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(magtable))

#Dictionary with literature reference star, and their RVB values & ID

# Sp.Type	  U-B	(B-V)	(V-R)C		
# B3.5	     -0.65	-0.19	-0.12

# So for NGC 869 1331 (spectral type B3.7V):
# V-R = ((B-V)*-0.12)/-0.19 = (0.353*-0.12)/-0.19 = 0.22294737
# R = V - (V-R) = 12.792 - 0.22294737 = 12.5690526

lit_refstars_RVB_id= {
    'NGC 869 593 ra34.5898315643 dec+56.9690227751' : [13.000, 13.050, 13.310, 704], #this star is a visual binary
    'BD+56 532 ra34.8497521072 dec+57.0044188540' : [10.228, 10.439, 10.937, 1393], #this star is a visual binary
    'NGC 869 1331 ra34.8966210059 dec+56.9960328997' : [12.5690526, 12.792, 13.145, 1527]} #B3.7V type MS star

#Plot B-V diagram w/ instrumental magnitudes
plt.scatter(magtable['magB_inst']-magtable['magVshort_inst'], magtable['magVshort_inst'], s=1)
plt.gca().invert_yaxis()
plt.xlabel("Uncalibrated B-V", fontsize=20)
plt.ylabel("V", fontsize=20)
plt.title('Color-magnitude diagram NGC869B', fontsize=20)
plt.show()

#Lists containing calibration factor per filter (one for every reference star)
C_Rlonglist = []
C_Rshortlist = []
C_Vlonglist = []
C_Vshortlist = []
C_Blist = []

#List of names referencestars taken from dictionary
lit_refstars_names = lit_refstars_RVB_id.keys()
#2D array of magnitudes&id of referencestars
lit_refstars_values = list(lit_refstars_RVB_id.values())


for i in range(len(lit_refstars_names)):
    #List of magnitudes&id of one referencestar
    lit_refstar_values = lit_refstars_values[i]
    #Find instrumental magnitudes of referencestar in magtable using star id
    referencestar = magtable[:][lit_refstar_values[3] - 1]
    
    #Calculate calibration constant: Minst,ref - Mlit,ref
    C_Rlong = referencestar[3] - lit_refstar_values[0]
    C_Rshort = referencestar[5] - lit_refstar_values[0]
    
    C_Vlong = referencestar[7] - lit_refstar_values[1]
    C_Vshort = referencestar[9] - lit_refstar_values[1]
    
    C_B = referencestar[11] - lit_refstar_values[2]
    
    #Append C_RVB values to list
    C_Rlonglist.append(C_Rlong)
    C_Rshortlist.append(C_Rshort)
    
    C_Vlonglist.append(C_Vlong)
    C_Vshortlist.append(C_Vshort)
    
    C_Blist.append(C_B)
    
#Print table with calibrated magnitudes per filter & exposure time (3 times, one for every 
#  literature reference star)
calmagtable = QTable([id_list , np.round(RA, 3), np.round(DEC, 3),
                      
                      aperture_sum_Rlong, max_aperture_Rlong, SNR_Rlong, magRlong_inst, magRlong_error,
                      magRlong_inst - C_Rlonglist[0], 
                      magRlong_inst - C_Rlonglist[1],
                      magRlong_inst - C_Rlonglist[2],
                      
                      aperture_sum_Rshort, max_aperture_Rshort, SNR_Rshort, magRshort_inst, magRshort_error,
                      magRshort_inst - C_Rshortlist[0],
                      magRshort_inst - C_Rshortlist[1],
                      magRshort_inst - C_Rshortlist[2],
                      
                      aperture_sum_Vlong, max_aperture_Vlong, SNR_Vlong, magVlong_inst, magVlong_error,
                      magVlong_inst - C_Vlonglist[0], 
                      magVlong_inst - C_Vlonglist[1],
                      magVlong_inst - C_Vlonglist[2],
                      
                      aperture_sum_Vshort, max_aperture_Vshort, SNR_Vshort, magVshort_inst, magVshort_error,
                      magVshort_inst - C_Vshortlist[0],
                      magVshort_inst - C_Vshortlist[1],
                      magVshort_inst - C_Vshortlist[2],
                      
                      aperture_sum_B, max_aperture_B, SNR_B, magB_inst, magB_error,
                      magB_inst - C_Blist[0], 
                      magB_inst - C_Blist[1],
                      magB_inst - C_Blist[2]],
                      names=('id', 'RA', 'DEC',
                             'aperture_sum_Rlong', 'max_aperture_Rlong', 'SNR_Rlong', 'magRlong_inst', 'magRlong_error',
                             'magRlong_ref1', 'magRlong_ref2', 'magRlong_ref3',
                             
                             'aperture_sum_Rshort', 'max_aperture_Rshort', 'SNR_Rshort', 'magRshort_inst', 'magRshort_error',
                             'magRshort_ref1', 'magRshort_ref2', 'magRshort_ref3',
                             
                             'aperture_sum_Vlong', 'max_aperture_Vlong', 'SNR_Vlong', 'magVlong_inst', 'magVlong_error',
                             'magVlong_ref1', 'magVlong_ref2', 'magVlong_ref3',
                             
                             'aperture_sum_Vshort', 'max_aperture_Vshort', 'SNR_Vshort', 'magVshort_inst', 'magVshort_error',
                             'magVshort_ref1', 'magVshort_ref2', 'magVshort_ref3',
                             
                             'aperture_sum_B', 'max_aperture_B', 'SNR_B', 'magB_inst', 'magB_error',
                             'magB_ref1', 'magB_ref2', 'magB_ref3'))
                                                                        
for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(calmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869B_calibrated.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = calmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(calmagtable))


###########################################NGC884A###########################################

#Read 1_Starfinder CSV to import instrumental magnitude values
StarfinderCSV_NGC884A = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884A.csv')

#Make lists of the instrumental magnitudes and their errors for every filter and exposure time
magB_inst = []
magB_error = []
aperture_sum_B = []
max_aperture_B = []
SNR_B = []

magRlong_inst = []
magRlong_error = []
aperture_sum_Rlong = []
max_aperture_Rlong = []
SNR_Rlong = []

magRshort_inst = []
magRshort_error = []
aperture_sum_Rshort = []
max_aperture_Rshort = []
SNR_Rshort = []

magVlong_inst = []
magVlong_error = []
aperture_sum_Vlong = []
max_aperture_Vlong = []
SNR_Vlong = []

magVshort_inst = []
magVshort_error = []
aperture_sum_Vshort = []
max_aperture_Vshort = []
SNR_Vshort = []

id_list = []

#Load RA and DEC coordinates and put into list
NGC884_A_WCSpositions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC884_A_WCSpositions.csv', delimiter = ',')

RA = []
DEC = []
for i in range(len(NGC884_A_WCSpositions)):
    RA.append(NGC884_A_WCSpositions[i, 0])
    DEC.append(NGC884_A_WCSpositions[i, 1])

#Fill empty lists with appropriate magnitude values
for i in range(len(StarfinderCSV_NGC884A['filtertype'])):
    if StarfinderCSV_NGC884A['filtertype'][i] == 'B':
        starid = StarfinderCSV_NGC884A['id'][i]
        uniqueid = f"{starid}_NGC884A"
        id_list.append(uniqueid)
        magB_inst.append(StarfinderCSV_NGC884A['mag inst'][i])
        magB_error.append(StarfinderCSV_NGC884A['mag error'][i])
        aperture_sum_B.append(StarfinderCSV_NGC884A['aperture_sum'][i])
        max_aperture_B.append(StarfinderCSV_NGC884A['max aperture'][i])
        SNR_B.append(StarfinderCSV_NGC884A['SNR'][i])
    
    if StarfinderCSV_NGC884A['filtertype'][i] == 'R':
        if StarfinderCSV_NGC884A['exposuretime'][i] == 60.0:
            magRlong_inst.append(StarfinderCSV_NGC884A['mag inst'][i])
            magRlong_error.append(StarfinderCSV_NGC884A['mag error'][i])
            aperture_sum_Rlong.append(StarfinderCSV_NGC884A['aperture_sum'][i])
            max_aperture_Rlong.append(StarfinderCSV_NGC884A['max aperture'][i])
            SNR_Rlong.append(StarfinderCSV_NGC884A['SNR'][i])
    
    if StarfinderCSV_NGC884A['filtertype'][i] == 'R':
        if StarfinderCSV_NGC884A['exposuretime'][i] == 5.0:
            magRshort_inst.append(StarfinderCSV_NGC884A['mag inst'][i])
            magRshort_error.append(StarfinderCSV_NGC884A['mag error'][i])
            aperture_sum_Rshort.append(StarfinderCSV_NGC884A['aperture_sum'][i])
            max_aperture_Rshort.append(StarfinderCSV_NGC884A['max aperture'][i])
            SNR_Rshort.append(StarfinderCSV_NGC884A['SNR'][i])
        
    if StarfinderCSV_NGC884A['filtertype'][i] == 'V':
        if StarfinderCSV_NGC884A['exposuretime'][i] == 60.0:
            magVlong_inst.append(StarfinderCSV_NGC884A['mag inst'][i])
            magVlong_error.append(StarfinderCSV_NGC884A['mag error'][i])
            aperture_sum_Vlong.append(StarfinderCSV_NGC884A['aperture_sum'][i])
            max_aperture_Vlong.append(StarfinderCSV_NGC884A['max aperture'][i])
            SNR_Vlong.append(StarfinderCSV_NGC884A['SNR'][i])
    
    if StarfinderCSV_NGC884A['filtertype'][i] == 'V':
        if StarfinderCSV_NGC884A['exposuretime'][i] == 5.0:
            magVshort_inst.append(StarfinderCSV_NGC884A['mag inst'][i])
            magVshort_error.append(StarfinderCSV_NGC884A['mag error'][i])
            aperture_sum_Vshort.append(StarfinderCSV_NGC884A['aperture_sum'][i])
            max_aperture_Vshort.append(StarfinderCSV_NGC884A['max aperture'][i])
            SNR_Vshort.append(StarfinderCSV_NGC884A['SNR'][i])

#Create a table with one row for every star containing the magnitude values and errors for the various filters
magtable = QTable([id_list, RA, DEC, 
                   magRlong_inst, magRlong_error, 
                   magRshort_inst, magRshort_error, 
                   magVlong_inst, magVlong_error, 
                   magVshort_inst, magVshort_error, 
                   magB_inst, magB_error], 
                   names=('id', 'RA', 'DEC', 
                          'magRlong_inst', 'magRlong_error',
                          'magRshort_inst', 'magRshort_error', 
                          'magVlong_inst', 'magVlong_error',
                          'magVshort_inst', 'magVshort_error', 
                          'magB_inst', 'magB_error'))    

for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(magtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884A.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = magtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(magtable))

# 'UCAC4 736-026265 ra35.6007227993 dec+57.0969145145' : [15.070, 15.290, 15.840, 1588]

#Dictionary with literature reference star, and their RVB values & ID
lit_refstars_RVB_id= {
    'UCAC4 737-025415 ra35.8470733222 dec+57.2048988721' : [14.710, 14.810, 15.390, 2111],
    'NGC 884 2413 ra35.5987214045 dec+57.0710730850' : [14.300, 14.311, 14.810, 1389],
    'NGC 884 2830 ra35.9116767562 dec+57.2508708742' : [13.190, 13.221, 13.801, 2277]}

#Plot B-V diagram w/ instrumental magnitudes
plt.scatter(magtable['magB_inst']-magtable['magVshort_inst'], magtable['magVshort_inst'], s=1)
plt.gca().invert_yaxis()
plt.xlabel("Uncalibrated B-V", fontsize=20)
plt.ylabel("V", fontsize=20)
plt.title('Color-magnitude diagram NGC884A', fontsize=20)
plt.show()

#Lists containing calibration factor per filter (one for every reference star)
C_Rlonglist = []
C_Rshortlist = []
C_Vlonglist = []
C_Vshortlist = []
C_Blist = []

#List of names referencestars taken from dictionary
lit_refstars_names = lit_refstars_RVB_id.keys()
#2D array of magnitudes&id of referencestars
lit_refstars_values = list(lit_refstars_RVB_id.values())


for i in range(len(lit_refstars_names)):
    #List of magnitudes&id of one referencestar
    lit_refstar_values = lit_refstars_values[i]
    #Find instrumental magnitudes of referencestar in magtable using star id
    referencestar = magtable[:][lit_refstar_values[3] - 1]
    
    #Calculate calibration constant: Minst,ref - Mlit,ref
    C_Rlong = referencestar[3] - lit_refstar_values[0]
    C_Rshort = referencestar[5] - lit_refstar_values[0]
    
    C_Vlong = referencestar[7] - lit_refstar_values[1]
    C_Vshort = referencestar[9] - lit_refstar_values[1]
    
    C_B = referencestar[11] - lit_refstar_values[2]
    
    #Append C_RVB values to list
    C_Rlonglist.append(C_Rlong)
    C_Rshortlist.append(C_Rshort)
    
    C_Vlonglist.append(C_Vlong)
    C_Vshortlist.append(C_Vshort)
    
    C_Blist.append(C_B)
    
#Print table with calibrated magnitudes per filter & exposure time (3 times, one for every 
#  literature reference star)
calmagtable = QTable([id_list , np.round(RA, 3), np.round(DEC, 3),
                      
                      aperture_sum_Rlong, max_aperture_Rlong, SNR_Rlong, magRlong_inst, magRlong_error,
                      magRlong_inst - C_Rlonglist[0], 
                      magRlong_inst - C_Rlonglist[1],
                      magRlong_inst - C_Rlonglist[2],
                      
                      aperture_sum_Rshort, max_aperture_Rshort, SNR_Rshort, magRshort_inst, magRshort_error,
                      magRshort_inst - C_Rshortlist[0],
                      magRshort_inst - C_Rshortlist[1],
                      magRshort_inst - C_Rshortlist[2],
                      
                      aperture_sum_Vlong, max_aperture_Vlong, SNR_Vlong, magVlong_inst, magVlong_error,
                      magVlong_inst - C_Vlonglist[0], 
                      magVlong_inst - C_Vlonglist[1],
                      magVlong_inst - C_Vlonglist[2],
                      
                      aperture_sum_Vshort, max_aperture_Vshort, SNR_Vshort, magVshort_inst, magVshort_error,
                      magVshort_inst - C_Vshortlist[0],
                      magVshort_inst - C_Vshortlist[1],
                      magVshort_inst - C_Vshortlist[2],
                      
                      aperture_sum_B, max_aperture_B, SNR_B, magB_inst, magB_error,
                      magB_inst - C_Blist[0], 
                      magB_inst - C_Blist[1],
                      magB_inst - C_Blist[2]],
                      names=('id', 'RA', 'DEC',
                             'aperture_sum_Rlong', 'max_aperture_Rlong', 'SNR_Rlong', 'magRlong_inst', 'magRlong_error',
                             'magRlong_ref1', 'magRlong_ref2', 'magRlong_ref3',
                             
                             'aperture_sum_Rshort', 'max_aperture_Rshort', 'SNR_Rshort', 'magRshort_inst', 'magRshort_error',
                             'magRshort_ref1', 'magRshort_ref2', 'magRshort_ref3',
                             
                             'aperture_sum_Vlong', 'max_aperture_Vlong', 'SNR_Vlong', 'magVlong_inst', 'magVlong_error',
                             'magVlong_ref1', 'magVlong_ref2', 'magVlong_ref3',
                             
                             'aperture_sum_Vshort', 'max_aperture_Vshort', 'SNR_Vshort', 'magVshort_inst', 'magVshort_error',
                             'magVshort_ref1', 'magVshort_ref2', 'magVshort_ref3',
                             
                             'aperture_sum_B', 'max_aperture_B', 'SNR_B', 'magB_inst', 'magB_error',
                             'magB_ref1', 'magB_ref2', 'magB_ref3'))
                                                                          
for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(calmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884A_calibrated.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = calmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(calmagtable))


###########################################NGC884B###########################################

#Read 1_Starfinder CSV to import instrumental magnitude values
StarfinderCSV_NGC884B = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884B.csv')

#Make lists of the instrumental magnitudes and their errors for every filter and exposure time
magB_inst = []
magB_error = []
aperture_sum_B = []
max_aperture_B = []
SNR_B = []

magRlong_inst = []
magRlong_error = []
aperture_sum_Rlong = []
max_aperture_Rlong = []
SNR_Rlong = []

magRshort_inst = []
magRshort_error = []
aperture_sum_Rshort = []
max_aperture_Rshort = []
SNR_Rshort = []

magVlong_inst = []
magVlong_error = []
aperture_sum_Vlong = []
max_aperture_Vlong = []
SNR_Vlong = []

magVshort_inst = []
magVshort_error = []
aperture_sum_Vshort = []
max_aperture_Vshort = []
SNR_Vshort = []

id_list = []
                  
#Load RA and DEC coordinates and put into list
NGC884_B_WCSpositions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC884_B_WCSpositions.csv', delimiter = ',')

RA = []
DEC = []
for i in range(len(NGC884_B_WCSpositions)):
    RA.append(NGC884_B_WCSpositions[i, 0])
    DEC.append(NGC884_B_WCSpositions[i, 1])

#Fill empty lists with appropriate magnitude values
for i in range(len(StarfinderCSV_NGC884B['filtertype'])):
    if StarfinderCSV_NGC884B['filtertype'][i] == 'B':
        starid = StarfinderCSV_NGC884B['id'][i]
        uniqueid = f"{starid}_NGC884B"
        id_list.append(uniqueid)
        magB_inst.append(StarfinderCSV_NGC884B['mag inst'][i])
        magB_error.append(StarfinderCSV_NGC884B['mag error'][i])
        aperture_sum_B.append(StarfinderCSV_NGC884B['aperture_sum'][i])
        max_aperture_B.append(StarfinderCSV_NGC884B['max aperture'][i])
        SNR_B.append(StarfinderCSV_NGC884B['SNR'][i])
    
    if StarfinderCSV_NGC884B['filtertype'][i] == 'R':
        if StarfinderCSV_NGC884B['exposuretime'][i] == 60.0:
            magRlong_inst.append(StarfinderCSV_NGC884B['mag inst'][i])
            magRlong_error.append(StarfinderCSV_NGC884B['mag error'][i])
            aperture_sum_Rlong.append(StarfinderCSV_NGC884B['aperture_sum'][i])
            max_aperture_Rlong.append(StarfinderCSV_NGC884B['max aperture'][i])
            SNR_Rlong.append(StarfinderCSV_NGC884B['SNR'][i])
    
    if StarfinderCSV_NGC884B['filtertype'][i] == 'R':
        if StarfinderCSV_NGC884B['exposuretime'][i] == 5.0:
            magRshort_inst.append(StarfinderCSV_NGC884B['mag inst'][i])
            magRshort_error.append(StarfinderCSV_NGC884B['mag error'][i])
            aperture_sum_Rshort.append(StarfinderCSV_NGC884B['aperture_sum'][i])
            max_aperture_Rshort.append(StarfinderCSV_NGC884B['max aperture'][i])
            SNR_Rshort.append(StarfinderCSV_NGC884B['SNR'][i])
         
    if StarfinderCSV_NGC884B['filtertype'][i] == 'V':
        if StarfinderCSV_NGC884B['exposuretime'][i] == 60.0:
            magVlong_inst.append(StarfinderCSV_NGC884B['mag inst'][i])
            magVlong_error.append(StarfinderCSV_NGC884B['mag error'][i])
            aperture_sum_Vlong.append(StarfinderCSV_NGC884B['aperture_sum'][i])
            max_aperture_Vlong.append(StarfinderCSV_NGC884B['max aperture'][i])
            SNR_Vlong.append(StarfinderCSV_NGC884B['SNR'][i])
    
    if StarfinderCSV_NGC884B['filtertype'][i] == 'V':
        if StarfinderCSV_NGC884B['exposuretime'][i] == 5.0:
            magVshort_inst.append(StarfinderCSV_NGC884B['mag inst'][i])
            magVshort_error.append(StarfinderCSV_NGC884B['mag error'][i])
            aperture_sum_Vshort.append(StarfinderCSV_NGC884B['aperture_sum'][i])
            max_aperture_Vshort.append(StarfinderCSV_NGC884B['max aperture'][i])
            SNR_Vshort.append(StarfinderCSV_NGC884B['SNR'][i])

#Create a table with one row for every star containing the magnitude values and errors for the various filters
magtable = QTable([id_list, RA, DEC, 
                    magRlong_inst, magRlong_error, 
                    magRshort_inst, magRshort_error, 
                    magVlong_inst, magVlong_error, 
                    magVshort_inst, magVshort_error, 
                    magB_inst, magB_error], 
                    names=('id', 'RA', 'DEC', 
                          'magRlong_inst', 'magRlong_error',
                          'magRshort_inst', 'magRshort_error', 
                          'magVlong_inst', 'magVlong_error',
                          'magVshort_inst', 'magVshort_error', 
                          'magB_inst', 'magB_error'))    

for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(magtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884B.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = magtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(magtable))

#Dictionary with literature reference star, and their RVB values & ID
lit_refstars_RVB_id= {
    'NGC 884 2615 ra35.7423023923 dec+57.0037244441' : [14.240, 14.344, 14.843, 1620],
    'NGC 884 1869 ra35.2764041132 dec+56.8866202266' : [12.492, 12.791, 13.643, 416], #visual binary
    'NGC 884 2413 ra35.5987214045 dec+57.0710730850' : [14.300, 14.311, 14.810, 1275]}

#Plot B-V diagram w/ instrumental magnitudes
plt.scatter(magtable['magB_inst']-magtable['magVshort_inst'], magtable['magVshort_inst'], s=1)
plt.gca().invert_yaxis()
plt.xlabel("Uncalibrated B-V", fontsize=20)
plt.ylabel("V", fontsize=20)
plt.title('Color-magnitude diagram NGC884B', fontsize=20)
plt.show()

#Lists containing calibration factor per filter (one for every reference star)
C_Rlonglist = []
C_Rshortlist = []
C_Vlonglist = []
C_Vshortlist = []
C_Blist = []

#List of names referencestars taken from dictionary
lit_refstars_names = lit_refstars_RVB_id.keys()
#2D array of magnitudes&id of referencestars
lit_refstars_values = list(lit_refstars_RVB_id.values())


for i in range(len(lit_refstars_names)):
    #List of magnitudes&id of one referencestar
    lit_refstar_values = lit_refstars_values[i]
    #Find instrumental magnitudes of referencestar in magtable using star id
    referencestar = magtable[:][lit_refstar_values[3] - 1]
    
    #Calculate calibration constant: Minst,ref - Mlit,ref
    C_Rlong = referencestar[3] - lit_refstar_values[0]
    C_Rshort = referencestar[5] - lit_refstar_values[0]
    
    C_Vlong = referencestar[7] - lit_refstar_values[1]
    C_Vshort = referencestar[9] - lit_refstar_values[1]
    
    C_B = referencestar[11] - lit_refstar_values[2]
    
    #Append C_RVB values to list
    C_Rlonglist.append(C_Rlong)
    C_Rshortlist.append(C_Rshort)
    
    C_Vlonglist.append(C_Vlong)
    C_Vshortlist.append(C_Vshort)
    
    C_Blist.append(C_B)
    
#Print table with calibrated magnitudes per filter & exposure time (3 times, one for every 
#  literature reference star)
calmagtable = QTable([id_list , np.round(RA, 3), np.round(DEC, 3),
                      
                      aperture_sum_Rlong, max_aperture_Rlong, SNR_Rlong, magRlong_inst, magRlong_error,
                      magRlong_inst - C_Rlonglist[0], 
                      magRlong_inst - C_Rlonglist[1],
                      magRlong_inst - C_Rlonglist[2],
                      
                      aperture_sum_Rshort, max_aperture_Rshort, SNR_Rshort, magRshort_inst, magRshort_error,
                      magRshort_inst - C_Rshortlist[0],
                      magRshort_inst - C_Rshortlist[1],
                      magRshort_inst - C_Rshortlist[2],
                      
                      aperture_sum_Vlong, max_aperture_Vlong, SNR_Vlong, magVlong_inst, magVlong_error,
                      magVlong_inst - C_Vlonglist[0], 
                      magVlong_inst - C_Vlonglist[1],
                      magVlong_inst - C_Vlonglist[2],
                      
                      aperture_sum_Vshort, max_aperture_Vshort, SNR_Vshort, magVshort_inst, magVshort_error,
                      magVshort_inst - C_Vshortlist[0],
                      magVshort_inst - C_Vshortlist[1],
                      magVshort_inst - C_Vshortlist[2],
                      
                      aperture_sum_B, max_aperture_B, SNR_B, magB_inst, magB_error,
                      magB_inst - C_Blist[0], 
                      magB_inst - C_Blist[1],
                      magB_inst - C_Blist[2]],
                      names=('id', 'RA', 'DEC',
                             'aperture_sum_Rlong', 'max_aperture_Rlong', 'SNR_Rlong', 'magRlong_inst', 'magRlong_error',
                             'magRlong_ref1', 'magRlong_ref2', 'magRlong_ref3',
                             
                             'aperture_sum_Rshort', 'max_aperture_Rshort', 'SNR_Rshort', 'magRshort_inst', 'magRshort_error',
                             'magRshort_ref1', 'magRshort_ref2', 'magRshort_ref3',
                             
                             'aperture_sum_Vlong', 'max_aperture_Vlong', 'SNR_Vlong', 'magVlong_inst', 'magVlong_error',
                             'magVlong_ref1', 'magVlong_ref2', 'magVlong_ref3',
                             
                             'aperture_sum_Vshort', 'max_aperture_Vshort', 'SNR_Vshort', 'magVshort_inst', 'magVshort_error',
                             'magVshort_ref1', 'magVshort_ref2', 'magVshort_ref3',
                             
                             'aperture_sum_B', 'max_aperture_B', 'SNR_B', 'magB_inst', 'magB_error',
                             'magB_ref1', 'magB_ref2', 'magB_ref3'))
                                                                         
for col in magtable.colnames:
    if not col == 'id':
        magtable[col].info.format = '%.6g' 
    
print(calmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884B_calibrated.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = calmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(calmagtable))

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)





