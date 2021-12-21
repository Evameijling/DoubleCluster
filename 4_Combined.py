#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 12:34:15 2021

@author: evagmelichmeijling
"""

import numpy as np
from astropy.io import fits
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
from astropy.wcs import WCS
import astropy.units as u

start = datetime.now()

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'a')
DOUBLECLUSTERTXT.write("4_COMBINED:\n")

#Import master .fit files from each object&plane with B, Vlong, Vshort, Rlong, Rshort filter&exptime
alldata = glob.glob('/Users/evagmelichmeijling/OneDrive/Capstone/DoubleCluster/Platesolved/*')

#Create list with array for each masterframe
data = []
for masterframe in alldata:
    data.append(fits.getdata(masterframe)) 
    
header = []
for masterframe in alldata:
    header.append(fits.getheader(masterframe))

filters = []
for fitfile in alldata:
    parts = fitfile.split('/')
    name = parts[len(parts)-1].split('.')[0]
    filters.append(name)

MagCalCSV_NGC869A = '/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869A_calibrated.csv'
MagCalCSV_NGC869B = '/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC869B_calibrated.csv'
MagCalCSV_NGC884A = '/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884A_calibrated.csv'
MagCalCSV_NGC884B = '/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/3_Magcalibrate/NGC884B_calibrated.csv'

planes = [MagCalCSV_NGC869A, MagCalCSV_NGC869B, MagCalCSV_NGC884A, MagCalCSV_NGC884B]

MergedMagCalCSV = pd.DataFrame()

for plane in planes:
    MergedMagCalCSV = MergedMagCalCSV.append(pd.read_csv(plane))

MergedMagCalCSV.to_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/4_Combined/MergedMagCalCSV.csv', index=False)

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)