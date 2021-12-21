#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 13:05:48 2021

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
import seaborn as sns
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from mpl_toolkits.mplot3d import Axes3D

start = datetime.now()

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'a')
DOUBLECLUSTERTXT.write("5_FILTER:\n")

#Create Qtable from CSV, selecting either short or long exposure for R and V filter depending on the max aperture value
MergedMagCalCSV = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/4_Combined/MergedMagCalCSV.csv')

print("Number of stars before filtering:", len(MergedMagCalCSV))
DOUBLECLUSTERTXT.write(f"Total amount of stars before filtering long and short exposures: {len(MergedMagCalCSV)}\n")

# Count how many stars were overexposed in the R an V filter
counter_overexposed = 0
counter_nonoverexposed = 0
for i in range(len(MergedMagCalCSV)):
    if MergedMagCalCSV['max_aperture_Rlong'][i] > 50000:
        counter_overexposed = counter_overexposed + 1
    if MergedMagCalCSV['max_aperture_Rlong'][i] < 50000:
        counter_nonoverexposed = counter_nonoverexposed + 1
    
counter_overexposed = 0
counter_nonoverexposed = 0
for i in range(len(MergedMagCalCSV)):
    if MergedMagCalCSV['max_aperture_Vlong'][i] > 50000:
        counter_overexposed = counter_overexposed + 1
    if MergedMagCalCSV['max_aperture_Vlong'][i] < 50000:
        counter_nonoverexposed = counter_nonoverexposed + 1
    
ID = []
RA = []
DEC = []
parallax = []
parallax_err = []
disttoearth = []

magR_inst = []
aperture_sum_R = []
magR_error = []
SNR_R = []
magR_ref1 = []
magR_ref2 = []
magR_ref3 = []

magV_inst = []
aperture_sum_V = []
magV_error = []
SNR_V = []
magV_ref1 = []
magV_ref2 = []
magV_ref3 = []

magB_inst = []
aperture_sum_B = []
magB_error = []
SNR_B = []
magB_ref1 = []
magB_ref2 = []
magB_ref3 = []

for i in range(len(MergedMagCalCSV)):
    ID.append(MergedMagCalCSV['id'][i])
    RA.append(MergedMagCalCSV['RA'][i])
    DEC.append(MergedMagCalCSV['DEC'][i])
    parallax.append(None)
    parallax_err.append(None)
    disttoearth.append(None)
    
    if MergedMagCalCSV['max_aperture_Rlong'][i] > 50000:
        magR_inst.append(MergedMagCalCSV['magRshort_inst'][i])
        aperture_sum_R.append(MergedMagCalCSV['aperture_sum_Rshort'][i])
        magR_error.append(MergedMagCalCSV['magRshort_error'][i])
        SNR_R.append(MergedMagCalCSV['SNR_Rshort'][i])
        magR_ref1.append(MergedMagCalCSV['magRshort_ref1'][i])
        magR_ref2.append(MergedMagCalCSV['magRshort_ref2'][i])
        magR_ref3.append(MergedMagCalCSV['magRshort_ref3'][i])
    else: 
        magR_inst.append(MergedMagCalCSV['magRlong_inst'][i])
        aperture_sum_R.append(MergedMagCalCSV['aperture_sum_Rlong'][i])
        magR_error.append(MergedMagCalCSV['magRlong_error'][i])
        SNR_R.append(MergedMagCalCSV['SNR_Rlong'][i])
        magR_ref1.append(MergedMagCalCSV['magRlong_ref1'][i])
        magR_ref2.append(MergedMagCalCSV['magRlong_ref2'][i])
        magR_ref3.append(MergedMagCalCSV['magRlong_ref3'][i])
        
    if MergedMagCalCSV['max_aperture_Vlong'][i] > 50000:
        magV_inst.append(MergedMagCalCSV['magVshort_inst'][i])
        aperture_sum_V.append(MergedMagCalCSV['aperture_sum_Vshort'][i])
        magV_error.append(MergedMagCalCSV['magVshort_error'][i])
        SNR_V.append(MergedMagCalCSV['SNR_Vshort'][i])
        magV_ref1.append(MergedMagCalCSV['magVshort_ref1'][i])
        magV_ref2.append(MergedMagCalCSV['magVshort_ref2'][i])
        magV_ref3.append(MergedMagCalCSV['magVshort_ref3'][i])
    else: 
        magV_inst.append(MergedMagCalCSV['magVlong_inst'][i])
        aperture_sum_V.append(MergedMagCalCSV['aperture_sum_Vlong'][i])
        magV_error.append(MergedMagCalCSV['magVlong_error'][i])
        SNR_V.append(MergedMagCalCSV['SNR_Vlong'][i])
        magV_ref1.append(MergedMagCalCSV['magVlong_ref1'][i])
        magV_ref2.append(MergedMagCalCSV['magVlong_ref2'][i])
        magV_ref3.append(MergedMagCalCSV['magVlong_ref3'][i])
        
    magB_inst.append(MergedMagCalCSV['magB_inst'][i])
    aperture_sum_B.append(MergedMagCalCSV['aperture_sum_B'][i])
    magB_error.append(MergedMagCalCSV['magB_error'][i])
    SNR_B.append(MergedMagCalCSV['SNR_B'][i])
    magB_ref1.append(MergedMagCalCSV['magB_ref1'][i])
    magB_ref2.append(MergedMagCalCSV['magB_ref2'][i])
    magB_ref3.append(MergedMagCalCSV['magB_ref3'][i])
    
RVBcalmagtable = QTable([ID, RA, DEC, 
                         parallax, parallax_err, disttoearth,
                         magR_inst, aperture_sum_R, magR_error, SNR_R, magR_ref1, magR_ref2, magR_ref3,
                         magV_inst, aperture_sum_V, magV_error, SNR_V, magV_ref1, magV_ref2, magV_ref3,
                         magB_inst, aperture_sum_B, magB_error, SNR_B, magB_ref1, magB_ref2, magB_ref3],
                         names=('id', 'RA', 'DEC', 'parallax', 'parallax_err', 'disttoearth',
                                'magR_inst', 'aperture_sum_R', 'magR_error', 'SNR_R', 'magR_ref1', 'magR_ref2', 'magR_ref3',
                                'magV_inst', 'aperture_sum_V', 'magV_error', 'SNR_V', 'magV_ref1', 'magV_ref2', 'magV_ref3',
                                'magB_inst', 'aperture_sum_B', 'magB_error', 'SNR_B', 'magB_ref1', 'magB_ref2', 'magB_ref3'))
           
print("Number of stars after filtering long and short exposures:", len(RVBcalmagtable))
         
DOUBLECLUSTERTXT.write(f"Total amount of stars after filtering long and short: {len(RVBcalmagtable)}\n")
        
for col in RVBcalmagtable.colnames:
    if not col == 'id':
        if not col == 'parallax':
            if not col == 'parallax_err':   
                if not col == 'disttoearth':
                    RVBcalmagtable[col].info.format = '%.6g' 
    
print(RVBcalmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/5_Filter/RVBMergedMagCalCSV.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = RVBcalmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(RVBcalmagtable))
    
#Remove all stars with negative magnitude error value
negerrorlist = []

index = 0
for star in RVBcalmagtable:
    if star['magR_error'] < 0:
        negerrorlist.append(index)
    if star['magV_error'] < 0:
        negerrorlist.append(index)
    if star['magB_error'] < 0:
        negerrorlist.append(index)
    index = index + 1

RVBcalmagtable.remove_rows(negerrorlist)

DOUBLECLUSTERTXT.write(f"Number of stars after filtering neg magerror: {len(RVBcalmagtable)}\n")

print("Number of stars after filtering neg magerror:", len(RVBcalmagtable))

#Remove all stars from R and V filter which have a SNR value below 3
smallSNR = []

index = 0
for star in RVBcalmagtable:
    if star['SNR_R'] < 3:
        smallSNR.append(index)
    if star['SNR_V'] < 3:
        negerrorlist.append(index)
    index = index + 1

RVBcalmagtable.remove_rows(smallSNR)

DOUBLECLUSTERTXT.write(f"Number of stars after filtering R and V SNR under the 3: {len(RVBcalmagtable)}\n")

print("Number of stars after filtering R and V SNR under the 3:", len(RVBcalmagtable))

#Remove stars from overlapping planes depending on their SNR value 
overlappingstars = []

for i in range(len(RVBcalmagtable)-1):
    for j in range(i+1, len(RVBcalmagtable)):
        if RVBcalmagtable['RA'][i] == RVBcalmagtable['RA'][j]:
            if RVBcalmagtable['DEC'][i] == RVBcalmagtable['DEC'][j]:
                if RVBcalmagtable['SNR_R'][i] < RVBcalmagtable['SNR_R'][j]: 
                    RVBcalmagtable['magR_inst'][i] = RVBcalmagtable['magR_inst'][j]
                    RVBcalmagtable['aperture_sum_R'][i] = RVBcalmagtable['aperture_sum_R'][j]
                    RVBcalmagtable['magR_error'][i] = RVBcalmagtable['magR_error'][j]
                    RVBcalmagtable['SNR_R'][i] = RVBcalmagtable['SNR_R'][j]
                    RVBcalmagtable['magR_ref1'][i] = RVBcalmagtable['magR_ref1'][j]
                    RVBcalmagtable['magR_ref2'][i] = RVBcalmagtable['magR_ref2'][j]
                    RVBcalmagtable['magR_ref3'][i] = RVBcalmagtable['magR_ref3'][j]
                if RVBcalmagtable['SNR_V'][i] < RVBcalmagtable['SNR_V'][j]:
                    RVBcalmagtable['magV_inst'][i] = RVBcalmagtable['magV_inst'][j]
                    RVBcalmagtable['aperture_sum_V'][i] = RVBcalmagtable['aperture_sum_V'][j]
                    RVBcalmagtable['magV_error'][i] = RVBcalmagtable['magV_error'][j]
                    RVBcalmagtable['SNR_V'][i] = RVBcalmagtable['SNR_V'][j]
                    RVBcalmagtable['magV_ref1'][i] = RVBcalmagtable['magV_ref1'][j]
                    RVBcalmagtable['magV_ref2'][i] = RVBcalmagtable['magV_ref2'][j]
                    RVBcalmagtable['magV_ref3'][i] = RVBcalmagtable['magV_ref3'][j]
                if RVBcalmagtable['SNR_B'][i] < RVBcalmagtable['SNR_B'][j]:
                    RVBcalmagtable['magB_inst'][i] = RVBcalmagtable['magB_inst'][j]
                    RVBcalmagtable['aperture_sum_B'][i] = RVBcalmagtable['aperture_sum_B'][j]
                    RVBcalmagtable['magB_error'][i] = RVBcalmagtable['magB_error'][j]
                    RVBcalmagtable['SNR_B'][i] = RVBcalmagtable['SNR_B'][j]
                    RVBcalmagtable['magB_ref1'][i] = RVBcalmagtable['magB_ref1'][j]
                    RVBcalmagtable['magB_ref2'][i] = RVBcalmagtable['magB_ref2'][j]
                    RVBcalmagtable['magB_ref3'][i] = RVBcalmagtable['magB_ref3'][j]
                overlappingstars.append(j)

RVBcalmagtable.remove_rows(overlappingstars)

DOUBLECLUSTERTXT.write(f"Number of stars after removing overlap from planes: {len(RVBcalmagtable)}\n")

print("Number of stars after removing overlap:", len(RVBcalmagtable))

#Compare data to gaia data to add column with parallax value 
DoubleClusterGaia = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/5_Filter/GaiaData/DoubleCluster_0.7deg-result.csv')

DOUBLECLUSTERTXT.write(f"Number of stars found 0.7 degrees from center double cluster Gaia Archive: {len(DoubleClusterGaia)}\n")

RA = []
DEC = []
parallax = []
parallax_err = []
mag = []
for i in range(len(DoubleClusterGaia)):
    if DoubleClusterGaia['phot_g_mean_mag'][i] < 17:
        RA.append(DoubleClusterGaia['ra'][i]- 0.0057)
        DEC.append(DoubleClusterGaia['dec'][i] + 0.0029)
        parallax.append(DoubleClusterGaia['parallax'][i])
        parallax_err.append(DoubleClusterGaia['parallax_error'][i])
        mag.append(DoubleClusterGaia['phot_g_mean_mag'][i])

plt.scatter(RA, DEC, label='Gaia', s=0.5)
plt.scatter(RVBcalmagtable['RA'], RVBcalmagtable['DEC'], label='Observed data', s=0.5)
plt.legend(prop={'size': 7})
plt.show()

Gaiadata = QTable([RA, DEC, parallax, parallax_err, mag],
                   names=('RA', 'DEC', 'parallax', 'parallax_err', 'mag'))

for i in range(len(RVBcalmagtable)):
    print(i)
    r_list = []
    for j in range(len(Gaiadata)):
        # r = np.sqrt(((DoubleClusterGaia['ra'][j] - 0.0057) - RVBcalmagtable['RA'][i])**2 + ((DoubleClusterGaia['dec'][j] + 0.0029) - RVBcalmagtable['DEC'][i])**2)
        r = ((Gaiadata['RA'][j] - 0.0057) - RVBcalmagtable['RA'][i])**2 + ((Gaiadata['DEC'][j] + 0.0029) - RVBcalmagtable['DEC'][i])**2
        r_list.append(r)
    RVBcalmagtable['parallax'][i] = Gaiadata['parallax'][np.argmin(r_list)]
    RVBcalmagtable['parallax_err'][i] = Gaiadata['parallax_err'][np.argmin(r_list)]
    RVBcalmagtable['disttoearth'][i] = 1/(RVBcalmagtable['parallax'][i]/1000)

#Plot a 3d graph of the Double Cluster 
# fig = plt.figure()
fig = plt.figure(figsize = (25,25))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(RVBcalmagtable['RA'], RVBcalmagtable['DEC'], RVBcalmagtable['disttoearth'], c = 'b', marker='o', s=5*(15-RVBcalmagtable['magV_ref1']))
# for star in range(len(RVBcalmagtable)):
#     plt.text(RVBcalmagtable['RA'][star], RVBcalmagtable['DEC'][star], RVBcalmagtable['disttoearth'][star], RVBcalmagtable['id'][star], fontweight = 'bold')
ax.axes.set_zlim3d(bottom=0, top=3000) 
ax.set_xlabel('RA')
ax.set_ylabel('DEC')
ax.set_zlabel('Distance (pc)')
plt.show()

for col in RVBcalmagtable.colnames:
    if not col == 'id':
        if not col == 'parallax':
            if not col == 'parallax_err':   
                if not col == 'disttoearth':
                    RVBcalmagtable[col].info.format = '%.6g'  

posdisttoearth = []
for disttoearth in RVBcalmagtable['disttoearth']:
    if disttoearth > 0:
        posdisttoearth.append(disttoearth)

DOUBLECLUSTERTXT.write(f"Minimum distance found from Gaia Data: {np.min(posdisttoearth)} parsec\n")
DOUBLECLUSTERTXT.write(f"Maximum distance found from Gaia Data: {np.max(RVBcalmagtable['disttoearth'])} parsec\n")

print(RVBcalmagtable)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/5_Filter/FilteredRVBMergedMagCalCSV.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = RVBcalmagtable.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(RVBcalmagtable))
    
DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)

#### ALLE STERREN WAARBIJ DE APERTURE OVERLAPT ERUIT ####

   