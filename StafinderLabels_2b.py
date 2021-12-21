#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:10:09 2021

@author: evagmelichmeijling
"""

# %matplotlib auto --> interactieve plot window
# %matplotlib inline --> terug naar normaal

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
from VariablesDC import apertureradius, annulusradius_in, annulusradius_out

start = datetime.now()

#Import master .fit files from each object&plane with B, Vlong, Vshort, Rlong, Rshort filter&exptime
alldata = glob.glob('/Users/evagmelichmeijling/OneDrive/Capstone/DoubleCluster/Non-Platesolved/*')
alldata_pltsolved = glob.glob('/Users/evagmelichmeijling/OneDrive/Capstone/DoubleCluster/Platesolved/*')

#Create list with array for each masterframe
data = []
for masterframe in alldata:
    data.append(fits.getdata(masterframe)) 
    
#Create list with header information of each masterframe
header = []
for masterframe in alldata:
    header.append(fits.getheader(masterframe))

#List of all the included masterfiles
filters = []
for fitfile in alldata:
    parts = fitfile.split('/')
    name = parts[len(parts)-1].split('.')[0]
    filters.append(name)
    
def plot_image(image, title):
    vmin, vmax = np.percentile(image, [5, 95])
    #print(vmin, vmax)
    fig, ax1 = plt.subplots(1,1, figsize=(15,15))
    plt.imshow(image, cmap='gray' , norm=LogNorm(vmin=1*vmin,vmax=1.03*vmax))
    plt.title(title)
    plt.show()

#Apply starfinder on the Rlong filter of every plane. Save locations & plot the aperture and annulus rings around found objects
background = []
for i in range(len(data)):
    if 'R' in fits.getheader(alldata[i])['FILTER']: 
        if fits.getheader(alldata[i])['EXPTIME'] == 60.0:
            sigma_clip = SigmaClip(sigma=3.)
            bkg_estimator = MedianBackground()
            bkg = Background2D(data[i], (75, 75), filter_size=(3, 3),
                                sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
            background.append(bkg.background_median) 
            
            pedestrial = np.min(data[i] - bkg.background)
            bkg_sub = data[i] - bkg.background - pedestrial + 100
            
            # plot_image(data[i], ('data:', filters[i]))
            # plot_image(bkg.background, ('background', filters[i]))
            # plot_image(bkg_sub, ('bkg_sub:', filters[i]))
            
            #Save the x and y coordinate of the found objects for every plane
            mean, median, std = sigma_clipped_stats(bkg_sub, sigma=3.0)
            FWHM = fits.getheader(alldata_pltsolved[i])['FWHM']
            Sterrenvinder = DAOStarFinder(threshold=100, fwhm=FWHM+6, exclude_border=True) (bkg_sub[25:4000, 25:4000]) 
            Sterrenvinder['xcentroid'] #the xcoordinate of all the found stars
            Sterrenvinder['ycentroid'] #the ycoordinate of all the found stars
            for col in Sterrenvinder.colnames:
                Sterrenvinder[col].info.format = '%.8g'
            positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid'])) 
        
            print('Amount of stars found in', filters[i], ':', len(Sterrenvinder))
            print("mean:", mean, "median:", median, "std:", std)
            
            #Plot the Rlong masterfile with aperture and annulus rings
            aperture = CircularAperture(positions, r=apertureradius) 
            annulus_aperture = CircularAnnulus(positions, r_in=annulusradius_in, r_out=annulusradius_out)
            annulus_masks = annulus_aperture.to_mask(method='center')
            plt.figure(figsize=(20,20))
            aperture.plot(color='red', lw=1., alpha=1)
            annulus_aperture.plot(color='red', lw=0.5, alpha=0.5)
            vmin, vmax = np.percentile(bkg_sub, [5, 95])
            plt.imshow(bkg_sub[25:4000, 25:4000], cmap='gray', norm=LogNorm(vmin=1*vmin,vmax=1.03*vmax))
            for star in range(len(np.array(Sterrenvinder['id']))):
                    plt.text(Sterrenvinder['xcentroid'][star],Sterrenvinder['ycentroid'][star],Sterrenvinder['id'][star], fontweight = 'bold')
            plt.title(filters[i])
            plt.show()
        
print(datetime.now() - start)    
            
            