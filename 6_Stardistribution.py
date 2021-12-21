#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 10:12:08 2021

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

# %matplotlib auto
# %matplotlib inline 

start = datetime.now()

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'a')
DOUBLECLUSTERTXT.write("6_STARDISTRIBUTION:\n")

#Create Qtable from CSV and combine 3 literature reference stars to 1 & calculate final error 
DoubleClusterCSV = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/5_Filter/FilteredRVBMergedMagCalCSV.csv')

ID = []
RA = []
DEC = []
parallax = []
parallax_err = []
disttoearth = []

magR_inst = []
aperture_sum_R = []
Rerror = []
SNR_R = []
m_R = []

magV_inst = []
aperture_sum_V = []
Verror = []
SNR_V = []
m_V = []

magB_inst = []
aperture_sum_B = []
Berror = []
SNR_B = []
m_B = []

for i in range(len(DoubleClusterCSV)):
    ID.append(DoubleClusterCSV['id'][i])
    RA.append(DoubleClusterCSV['RA'][i])
    DEC.append(DoubleClusterCSV['DEC'][i])
    parallax.append(DoubleClusterCSV['parallax'][i])
    parallax_err.append(DoubleClusterCSV['parallax_err'][i])
    disttoearth.append(DoubleClusterCSV['disttoearth'][i])
    
    magR_inst.append(DoubleClusterCSV['magR_inst'][i])
    aperture_sum_R.append(DoubleClusterCSV['aperture_sum_R'][i])
    Rerror.append(DoubleClusterCSV['magR_error'][i]/np.sqrt(3))
    SNR_R.append(DoubleClusterCSV['SNR_R'][i])
    m_R.append((DoubleClusterCSV['magR_ref1'][i] + DoubleClusterCSV['magR_ref2'][i] + DoubleClusterCSV['magR_ref3'][i])/3)
    
    magV_inst.append(DoubleClusterCSV['magV_inst'][i])
    aperture_sum_V.append(DoubleClusterCSV['aperture_sum_V'][i])
    Verror.append(DoubleClusterCSV['magV_error'][i]/np.sqrt(3))
    SNR_V.append(DoubleClusterCSV['SNR_V'][i])
    m_V.append((DoubleClusterCSV['magV_ref1'][i] + DoubleClusterCSV['magV_ref2'][i] + DoubleClusterCSV['magV_ref3'][i])/3)

    magB_inst.append(DoubleClusterCSV['magB_inst'][i])
    aperture_sum_B.append(DoubleClusterCSV['aperture_sum_B'][i])
    Berror.append(DoubleClusterCSV['magB_error'][i]/np.sqrt(3))
    SNR_B.append(DoubleClusterCSV['SNR_B'][i])
    m_B.append((DoubleClusterCSV['magB_ref1'][i] + DoubleClusterCSV['magB_ref2'][i] + DoubleClusterCSV['magB_ref3'][i])/3)

DoubleCluster = QTable([ID, RA, DEC, 
                        parallax, parallax_err, disttoearth,
                        magR_inst, aperture_sum_R, Rerror, SNR_R, m_R,
                        magV_inst, aperture_sum_V, Verror, SNR_V, m_V,
                        magB_inst, aperture_sum_B, Berror, SNR_B, m_B],
                        names=('id', 'RA', 'DEC', 'parallax', 'parallax_err', 'disttoearth',
                               'magR_inst', 'aperture_sum_R', 'Rerror', 'SNR_R', 'm_R',
                               'magV_inst', 'aperture_sum_V', 'Verror', 'SNR_V', 'm_V',
                               'magB_inst', 'aperture_sum_B', 'Berror', 'SNR_B', 'm_B'))

for col in DoubleCluster.colnames:
    if not col == 'id':
        if not col == 'parallax':
            if not col == 'parallax_err':   
                if not col == 'disttoearth':
                    DoubleCluster[col].info.format = '%.6g' 

DOUBLECLUSTERTXT.write(f"Lowest m_R: {np.max(DoubleCluster['m_R'])}\n")    
DOUBLECLUSTERTXT.write(f"Lowest m_V: {np.max(DoubleCluster['m_V'])}\n")    
DOUBLECLUSTERTXT.write(f"Lowest m_B: {np.max(DoubleCluster['m_B'])}\n")    

print(DoubleCluster)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/6_Stardistribution/DoubleCluster.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = DoubleCluster.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(DoubleCluster))


#Create Qtable of only NGC869
NGC869 = QTable([ID, RA, DEC, parallax, parallax_err, disttoearth,
                    magR_inst, aperture_sum_R, Rerror, SNR_R, m_R,
                    magV_inst, aperture_sum_V, Verror, SNR_V, m_V,
                    magB_inst, aperture_sum_B, Berror, SNR_B, m_B],
                    names=('id', 'RA', 'DEC', 'parallax', 'parallax_err', 'disttoearth',
                          'magR_inst', 'aperture_sum_R', 'Rerror', 'SNR_R', 'm_R',
                          'magV_inst', 'aperture_sum_V', 'Verror', 'SNR_V', 'm_V',
                          'magB_inst', 'aperture_sum_B', 'Berror', 'SNR_B', 'm_B'))
deletelist_NGC869 = []
for i in range(len(NGC869)):
    if 'NGC884' in NGC869['id'][i]:
        deletelist_NGC869.append(i)
NGC869.remove_rows(deletelist_NGC869)

#Create heatmap of stellar distribution NGC869
sns.set_style("white")
sns.kdeplot(x=NGC869['RA'], y=NGC869['DEC'], n_levels=30, cmap="Blues_d")
sns.kdeplot(x=NGC869['RA'], y=NGC869['DEC'],
            color='b', shade=True,
            cmap="Blues", thresh=0, n_levels=200)
plt.xlabel("RA")
plt.ylabel("DEC")
plt.title('Stellar distribution NGC 869')
plt.show()

#Create Qtable of only NGC884
NGC884 = QTable([ID, RA, DEC, parallax, parallax_err, disttoearth,
                    magR_inst, aperture_sum_R, Rerror, SNR_R, m_R,
                    magV_inst, aperture_sum_V, Verror, SNR_V, m_V,
                    magB_inst, aperture_sum_B, Berror, SNR_B, m_B],
                    names=('id', 'RA', 'DEC', 'parallax', 'parallax_err', 'disttoearth',
                          'magR_inst', 'aperture_sum_R', 'Rerror', 'SNR_R', 'm_R',
                          'magV_inst', 'aperture_sum_V', 'Verror', 'SNR_V', 'm_V',
                          'magB_inst', 'aperture_sum_B', 'Berror', 'SNR_B', 'm_B'))
deletelist_NGC884 = []
for i in range(len(NGC884)):
    if 'NGC869' in NGC884['id'][i]:
        deletelist_NGC884.append(i)
NGC884.remove_rows(deletelist_NGC884)

#Create heatmap of stellar distribution NGC869
sns.set_style("white")
sns.kdeplot(x=NGC884['RA'], y=NGC884['DEC'], n_levels=30, cmap='Purples_d')
sns.kdeplot(x=NGC884['RA'], y=NGC884['DEC'],
            color='m', shade=True,
            cmap="Purples", thresh=0, n_levels=200)
plt.xlabel("RA")
plt.ylabel("DEC")
plt.title('Stellar distribution NGC 884')
plt.show()

#Create heatmap of stellar distribution Double Cluster
sns.kdeplot(x=NGC869['RA'], y=NGC869['DEC'], n_levels=30, cmap="Blues_d")
sns.kdeplot(x=NGC884['RA'], y=NGC884['DEC'], n_levels=30, cmap='Purples_d')
sns.kdeplot(x=NGC869['RA'], y=NGC869['DEC'],
            color='b', shade=True,
            cmap="Blues", shade_lowest=False)
sns.kdeplot(x=NGC884['RA'], y=NGC884['DEC'],
            color='m', shade=True,
            cmap="Purples", shade_lowest=False) 
plt.title('Stellar distribution of the Double Cluster')
plt.show()

#Manually find the center of the two heatmaps & extract central star from Qtable
for i in range(len(DoubleCluster)):
    if '1641_NGC869A' in DoubleCluster['id'][i]:
        NGC869_center_index = i

for i in range(len(DoubleCluster)):
    if '1188_NGC884A' in DoubleCluster['id'][i]:
        NGC884_center_index = i
        
NGC869_center = DoubleCluster[NGC869_center_index]
NGC884_center = DoubleCluster[NGC884_center_index]

#Convert RA and DEC to SkyCoordinates 
NGC869_center_coord = SkyCoord(NGC869_center['RA'], NGC869_center['DEC'], frame='icrs', unit="deg")
NGC884_center_coord = SkyCoord(NGC884_center['RA'], NGC884_center['DEC'], frame='icrs', unit="deg")

DOUBLECLUSTERTXT.write(f"Center coordinates of NGC869: {NGC869_center['RA']} {NGC869_center['DEC']}\n")   
DOUBLECLUSTERTXT.write(f"Center coordinates of NGC884: {NGC884_center['RA']} {NGC884_center['DEC']}\n")   

#Determine distance from every star to the center & add to Qtable
disttocenter = []
for i in range(len(DoubleCluster)):
    if 'NGC869' in DoubleCluster['id'][i]:
        disttocenter.append((NGC869_center_coord.separation(SkyCoord(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], frame='icrs', unit="deg"))).degree)
    if 'NGC884' in DoubleCluster['id'][i]:
        disttocenter.append((NGC884_center_coord.separation(SkyCoord(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], frame='icrs', unit="deg"))).degree)

DoubleCluster = QTable([ID, RA, DEC, 
                        parallax, parallax_err, disttoearth, disttocenter,
                        magR_inst, aperture_sum_R, Rerror, SNR_R, m_R,
                        magV_inst, aperture_sum_V, Verror, SNR_V, m_V,
                        magB_inst, aperture_sum_B, Berror, SNR_B, m_B],
                        names=('id', 'RA', 'DEC', 
                               'parallax', 'parallax_err', 'disttoearth', 'disttocenter',
                               'magR_inst', 'aperture_sum_R', 'Rerror', 'SNR_R', 'm_R',
                               'magV_inst', 'aperture_sum_V', 'Verror', 'SNR_V', 'm_V',
                               'magB_inst', 'aperture_sum_B', 'Berror', 'SNR_B', 'm_B'))

for col in DoubleCluster.colnames:
    if not col == 'id':
        if not col == 'parallax':
            if not col == 'parallax_err':   
                if not col == 'disttoearth':
                    DoubleCluster[col].info.format = '%.6g' 
    
print(DoubleCluster)

#Export magtable information to CSV file
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/6_Stardistribution/FinalDoubleCluster.csv', 'w') as file:
    writer = csv.writer(file)

    headerlist = DoubleCluster.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()

    writer.writerows(np.array(DoubleCluster))
    
#Create scatter plot of the stellar distribution of double cluster, with different shade of color for various radii from center
for i in range(len(DoubleCluster)):
    if 'NGC869' in DoubleCluster['id'][i]:
        if DoubleCluster['disttocenter'][i] < 0.1:
            first869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='darkblue', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.1 and DoubleCluster['disttocenter'][i] < 0.2:
            second869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='royalblue',  s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.2 and DoubleCluster['disttocenter'][i] < 0.3:
            third869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='steelblue', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.3 and DoubleCluster['disttocenter'][i] < 0.4:
            fourth869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='lightskyblue', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.4 and DoubleCluster['disttocenter'][i] < 0.5:
            fifth869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='lightblue', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.5 and DoubleCluster['disttocenter'][i] < 0.6:
            sixth869 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='lightcyan', s=1.5, marker="*")
for i in range(len(DoubleCluster)):
    if 'NGC884' in DoubleCluster['id'][i]:
        if DoubleCluster['disttocenter'][i] < 0.1:
            first884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='indigo', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.1 and DoubleCluster['disttocenter'][i] < 0.2:
            second884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='rebeccapurple',  s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.2 and DoubleCluster['disttocenter'][i] < 0.3:
            third884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='darkviolet', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.3 and DoubleCluster['disttocenter'][i] < 0.4:
            fourth884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='mediumorchid', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.4 and DoubleCluster['disttocenter'][i] < 0.5:
            fifth884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='plum', s=1.5, marker="*")
        if DoubleCluster['disttocenter'][i] > 0.5 and DoubleCluster['disttocenter'][i] < 0.6:
            sixth884 = plt.scatter(DoubleCluster['RA'][i], DoubleCluster['DEC'][i], color='thistle', s=1.5, marker="*")

first_legend = plt.legend((first869, second869, third869, fourth869, fifth869, sixth869), ('0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6'), 
               markerscale=3, fontsize=6, title='Degrees from centre', title_fontsize=6, loc='lower left')
plt.gca().add_artist(first_legend)
plt.legend((first884, second884, third884, fourth884, fifth884, sixth884), ('0.0-0.1', '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6'), 
               markerscale=3, fontsize=6, title='Degrees from centre', title_fontsize=6, loc='lower right')

plt.xlabel("RA")
plt.ylabel("DEC")
plt.title('Stellar distribution of The Double Cluster')
plt.show()

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)
