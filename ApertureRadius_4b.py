#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 11:02:43 2021

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

##CHOOSE FILTER TYPE & LONG/SHORT##
filtertype = 'Vlong'
filtertype_compare = 'B'

MergedMagCalCSV_20 = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/Hulpprogrammas/ApertureRadius/MergedMagCalCSV_20.csv')
MergedMagCalCSV_15 = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/Hulpprogrammas/ApertureRadius/MergedMagCalCSV_15.csv')
MergedMagCalCSV_10 = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/Hulpprogrammas/ApertureRadius/MergedMagCalCSV_10.csv')
MergedMagCalCSV_07 = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/Hulpprogrammas/ApertureRadius/MergedMagCalCSV_07.csv')

SNR_20 = []
SNR_15 = []
SNR_10 = []
SNR_07 = []

instmag_20 = []
instmag_15 = []
instmag_10 = []
instmag_07 = []

instmag_compare_20 = []
instmag_compare_15 = []
instmag_compare_10 = []
instmag_compare_07 = []

SNR_20_min_SNR_15 = []
SNR_20_min_SNR_10 = []
SNR_20_min_SNR_07 = []

for i in range(len(MergedMagCalCSV_20['id'])):
    SNR_20.append(MergedMagCalCSV_20[f"SNR_{filtertype}"][i])
    SNR_15.append(MergedMagCalCSV_15[f"SNR_{filtertype}"][i])
    SNR_10.append(MergedMagCalCSV_10[f"SNR_{filtertype}"][i])
    SNR_07.append(MergedMagCalCSV_07[f"SNR_{filtertype}"][i])

    instmag_20.append(MergedMagCalCSV_20[f"mag{filtertype}_inst"][i])
    instmag_15.append(MergedMagCalCSV_15[f"mag{filtertype}_inst"][i])
    instmag_10.append(MergedMagCalCSV_10[f"mag{filtertype}_inst"][i])
    instmag_07.append(MergedMagCalCSV_07[f"mag{filtertype}_inst"][i])
    
    instmag_compare_20.append(MergedMagCalCSV_20[f"mag{filtertype_compare}_inst"][i])
    instmag_compare_15.append(MergedMagCalCSV_15[f"mag{filtertype_compare}_inst"][i])
    instmag_compare_10.append(MergedMagCalCSV_10[f"mag{filtertype_compare}_inst"][i])
    instmag_compare_07.append(MergedMagCalCSV_07[f"mag{filtertype_compare}_inst"][i])
   
    SNR_20_min_SNR_15.append(MergedMagCalCSV_20[f"SNR_{filtertype}"][i] - MergedMagCalCSV_15[f"SNR_{filtertype}"][i])
    SNR_20_min_SNR_10.append(MergedMagCalCSV_20[f"SNR_{filtertype}"][i] - MergedMagCalCSV_10[f"SNR_{filtertype}"][i])
    SNR_20_min_SNR_07.append(MergedMagCalCSV_20[f"SNR_{filtertype}"][i] - MergedMagCalCSV_07[f"SNR_{filtertype}"][i])


#SNR voor verschillende radius grootte
plt.hist(SNR_20, bins=10000, label = 'radius 20')
plt.hist(SNR_15, bins=10000, label = 'radius 15')
plt.hist(SNR_10, bins=10000, label = 'radius 10')
plt.hist(SNR_07, bins=10000, label = 'radius 07')
plt.title('The SNR distribution')
plt.xlabel('SNR')
plt.ylabel('Amount of stars with this SNR')
plt.xlim(0, 100)
plt.legend()
plt.show()

#SNR voor verschillende radius grootte vergelijken
plt.scatter(SNR_20, SNR_15)
plt.title('SNR met radius 20 vs SNR met radius 15')
plt.xlabel('SNR_20')
plt.ylabel('SNR_15')
plt.show()

#Toename in SNR vs helderheid(flux)
plt.scatter(instmag_20, SNR_20_min_SNR_15, label = 'SNR 20 - SNR 15', s = 5)
plt.scatter(instmag_20, SNR_20_min_SNR_10, label = 'SNR 20 - SNR 10', s = 5)
plt.scatter(instmag_20, SNR_20_min_SNR_07, label = 'SNR 20 - SNR 07', s = 5)
plt.title('Ik weet niet zo goed wat dit moet voorstellen?')
plt.xlabel('instmag_20')
plt.ylabel('SNR_2 minus SNR_15/10/07')
plt.legend()
plt.show()

#Mag-Mag diagram --> the bigger the fan, the less accurate the magnitude determination
plt.scatter(instmag_20, instmag_compare_20, label = 'radius 20', s=1)
plt.scatter(instmag_15, instmag_compare_15, label = 'radius 15', s=1)
plt.scatter(instmag_10, instmag_compare_10, label = 'radius 10', s=1)
plt.scatter(instmag_07, instmag_compare_07, label = 'radius 07', s=1)
plt.gca().invert_yaxis()
plt.xlabel(f"mag{filtertype}_inst")
plt.ylabel(f"mag{filtertype_compare}_inst")
plt.title('Mag-Mag diagram')
plt.legend()
plt.show()


