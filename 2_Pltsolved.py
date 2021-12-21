#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 15:03:26 2021

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
DOUBLECLUSTERTXT.write("2_PLTSOLVED:\n")

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
        
NGC869_A_positions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC869_A_positions.csv', delimiter = ',')
NGC869_B_positions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC869_B_positions.csv', delimiter = ',')
NGC884_A_positions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC884_A_positions.csv', delimiter = ',')
NGC884_B_positions = np.loadtxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC884_B_positions.csv', delimiter = ',')

NGC869_A_WCSpositions = []
NGC869_B_WCSpositions = []
NGC884_A_WCSpositions = []
NGC884_B_WCSpositions = []

for i in range(len(data)):
    if 'R' in fits.getheader(alldata[i])['FILTER']: 
        if fits.getheader(alldata[i])['EXPTIME'] == 60.0:
            if 'NGC869_A' in filters[i]:
                for star in range(len(NGC869_A_positions)):
                    f = fits.open(alldata[i])
                    w = WCS(f[0].header)
                    skycoordinates = []
                    sky = w.pixel_to_world(NGC869_A_positions[star,0], NGC869_A_positions[star,1])
                    skycoordinates.append(sky.ra.degree)
                    skycoordinates.append(sky.dec.degree)
                    NGC869_A_WCSpositions.append(skycoordinates)
                DOUBLECLUSTERTXT.write(f"Centercoordinate of {filters[i]}: {fits.getheader(alldata[i])['OBJCTRA']} {fits.getheader(alldata[i])['OBJCTDEC']} \n")
            if 'NGC869_B' in filters[i]:
                for star in range(len(NGC869_B_positions)):
                    f = fits.open(alldata[i])
                    w = WCS(f[0].header)
                    skycoordinates = []
                    sky = w.pixel_to_world(NGC869_B_positions[star,0], NGC869_B_positions[star,1])
                    skycoordinates.append(sky.ra.degree)
                    skycoordinates.append(sky.dec.degree)
                    NGC869_B_WCSpositions.append(skycoordinates)
                DOUBLECLUSTERTXT.write(f"Centercoordinate of {filters[i]}: {fits.getheader(alldata[i])['OBJCTRA']} {fits.getheader(alldata[i])['OBJCTDEC']} \n")
            if 'NGC884_A' in filters[i]:
                for star in range(len(NGC884_A_positions)):
                    f = fits.open(alldata[i])
                    w = WCS(f[0].header)
                    skycoordinates = []
                    sky = w.pixel_to_world(NGC884_A_positions[star,0], NGC884_A_positions[star,1])
                    skycoordinates.append(sky.ra.degree)
                    skycoordinates.append(sky.dec.degree)
                    NGC884_A_WCSpositions.append(skycoordinates)
                DOUBLECLUSTERTXT.write(f"Centercoordinate of {filters[i]}: {fits.getheader(alldata[i])['OBJCTRA']} {fits.getheader(alldata[i])['OBJCTDEC']} \n")
            if 'NGC884_B' in filters[i]:
                for star in range(len(NGC884_B_positions)):
                    f = fits.open(alldata[i])
                    w = WCS(f[0].header)
                    skycoordinates = []
                    sky = w.pixel_to_world(NGC884_B_positions[star,0], NGC884_B_positions[star,1])
                    skycoordinates.append(sky.ra.degree)
                    skycoordinates.append(sky.dec.degree)
                    NGC884_B_WCSpositions.append(skycoordinates)
                DOUBLECLUSTERTXT.write(f"Centercoordinate of {filters[i]}: {fits.getheader(alldata[i])['OBJCTRA']} {fits.getheader(alldata[i])['OBJCTDEC']} \n")
                    
np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC869_A_WCSpositions.csv', NGC869_A_WCSpositions, delimiter = ',')
np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC869_B_WCSpositions.csv', NGC869_B_WCSpositions, delimiter = ',')
np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC884_A_WCSpositions.csv', NGC884_A_WCSpositions, delimiter = ',')
np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/2_Pltsolved/NGC884_B_WCSpositions.csv', NGC884_B_WCSpositions, delimiter = ',')

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)

    