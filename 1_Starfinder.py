#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  3 11:24:38 2021

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

DOUBLECLUSTERTXT = open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/AnalysisResults.txt', 'a')
DOUBLECLUSTERTXT.write("1_STARFINDER:\n")

#Import master .fit files from each object&plane with B, Vlong, Vshort, Rlong, Rshort filter&exptime
alldata = glob.glob('/Users/evagmelichmeijling/OneDrive/Capstone/DoubleCluster/Non-Platesolved/*')
alldata_pltsolved = glob.glob('/Users/evagmelichmeijling/OneDrive/Capstone/DoubleCluster/Platesolved/*')

DOUBLECLUSTERTXT.write(f"Amount of masterframes: {len(alldata)}\n")

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

#Print information and plot of every masterframe to check anormalities

background = []
for i in range(len(data)):
    mean, median, std = sigma_clipped_stats(data[i], sigma=3.0)
    print("")
    print(filters[i])
    print("Measurements", np.shape(data[i]))
    print("Mean data", mean) 
    print("Standard deviation", std) 
    plot_image(data[i], filters[i])
    
    #Background and background noise estimation
    sigma_clip = SigmaClip(sigma=3.)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data[i], (75, 75), filter_size=(3, 3), 
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    background.append(bkg.background_median) 
    plot_image(bkg.background, ('background', filters[i]))

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
            
            plot_image(data[i], ('data:', filters[i]))
            plot_image(bkg.background, ('background', filters[i]))
            plot_image(bkg_sub, ('bkg_sub:', filters[i]))
            
            #Save the x and y coordinate of the found objects for every plane
            mean, median, std = sigma_clipped_stats(bkg_sub, sigma=3.0)
            FWHM = fits.getheader(alldata_pltsolved[i])['FWHM']
            Sterrenvinder = DAOStarFinder(threshold=100, fwhm=FWHM+6, exclude_border=True) (bkg_sub[25:4000, 25:4000]) 
            Sterrenvinder['xcentroid'] #the xcoordinate of all the found stars
            Sterrenvinder['ycentroid'] #the ycoordinate of all the found stars
            for col in Sterrenvinder.colnames:
                Sterrenvinder[col].info.format = '%.8g'
            positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid'])) 
        
            if 'NGC869_A' in filters[i]:
                NGC869_A_positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid']))
                np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC869_A_positions.csv', NGC869_A_positions, delimiter = ',')
                DOUBLECLUSTERTXT.write(f"Starfinder found {len(NGC869_A_positions)} stars for {filters[i]}\n")
            if 'NGC869_B' in filters[i]:
                NGC869_B_positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid']))
                np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC869_B_positions.csv', NGC869_B_positions, delimiter = ',')
                DOUBLECLUSTERTXT.write(f"Starfinder found {len(NGC869_B_positions)} stars for {filters[i]}\n")
            if 'NGC884_A' in filters[i]:
                NGC884_A_positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid']))
                np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC884_A_positions.csv', NGC884_A_positions, delimiter = ',')
                DOUBLECLUSTERTXT.write(f"Starfinder found {len(NGC884_A_positions)} stars for {filters[i]}\n")
            if 'NGC884_B' in filters[i]:
                NGC884_B_positions = np.transpose((Sterrenvinder['xcentroid'], Sterrenvinder['ycentroid'])) 
                np.savetxt('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/StarfinderPositionsPerPlane/NGC884_B_positions.csv', NGC884_B_positions, delimiter = ',')
                DOUBLECLUSTERTXT.write(f"Starfinder found {len(NGC884_B_positions)} stars for {filters[i]}\n")
            
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
            plt.title(filters[i])
            plt.show()

###########################################NGC869A###########################################

#Apply photometry on all filters of every plane, using the location of objects found in the Rlong filter 
phot869A_allfilters = []
for i in range(len(filters)):
    if 'NGC869_A' in filters[i]:
        #Aperture and annulus rings of masterfiles
        aperture = CircularAperture(NGC869_A_positions, r=apertureradius) 
        annulus_aperture = CircularAnnulus(NGC869_A_positions, r_in=annulusradius_in, r_out=annulusradius_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        ##bkg_mode = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data[i][25:4000, 25:4000])   
            annulus_data_1d = annulus_data[mask.data > 0] 
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            ##mode = stats.mode(annulus_data_1d)
            ##bkg_mode.append(mode[1])
        bkg_median = np.array(bkg_median)
        ##bkg_mode = np.array(bkg_mode)
        
        max_aperture = []
        aperture_masks = aperture.to_mask(method='center')
        for mask in aperture_masks:
            aperture_data = mask.multiply(data[i][25:4000, 25:4000])
            max_aperture_value = np.amax(aperture_data)
            max_aperture.append(max_aperture_value)
         
        #Create phottable for every filter of object plane 
        phot869A = aperture_photometry(data[i][25:4000, 25:4000], aperture)
        phot869A['object'] = filters[i]
        filtertype = fits.getheader(alldata[i])['FILTER']
        phot869A['filtertype'] = filtertype
        exposuretime = fits.getheader(alldata[i])['EXPTIME']
        phot869A['exposuretime'] = exposuretime
        phot869A['max aperture'] = max_aperture
        phot869A['annulus_median'] = bkg_median
        phot869A['aperture_bkg'] = bkg_median * aperture.area
        phot869A['aperture_sum_bkgsub'] = phot869A['aperture_sum'] - phot869A['aperture_bkg']
        phot869A['SNR'] = phot869A['aperture_sum_bkgsub'] / (phot869A['aperture_sum_bkgsub'] + phot869A['aperture_bkg'])**0.5
        phot869A['Absolute error'] = phot869A['aperture_sum_bkgsub'] / phot869A['SNR']
        phot869A['mag inst'] = -2.5 * np.log10(np.abs(phot869A['aperture_sum_bkgsub']))
        phot869A['mag error'] = 2.5 * 0.434 * (1/phot869A['SNR']) #alleen afhankelijk van de SNR
        for col in phot869A.colnames:
            if not col == 'filtertype':
                if not col == 'exposuretime':
                    if not col == 'object':
                        phot869A[col].info.format = '%.6g'  #for consistent table output
        print(phot869A)
        print("")
        
        phot869A_allfilters.append(phot869A)

#Create a CSV file and add headers (corresponding to the phottable headers)           
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869A.csv', 'w') as file:
    writer = csv.writer(file)
    
    headerlist = phot869A.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()
    
    #Add phottables to CSV file
    for filter in phot869A_allfilters:
        phot869A_array = np.array(filter)
        writer.writerows(phot869A_array)

#Sort the csv file so every object is grouped with each a different filter & exposure time 
sortedNGC869A = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869A.csv')
sortedNGC869A.sort_values(['id'], axis=0, inplace=True)
sortedNGC869A.to_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869A.csv', index=False)


###########################################NGC869B###########################################

#Apply photometry on all filters of every plane, using the location of objects found in the Rlong filter 
phot869B_allfilters = []
for i in range(len(filters)):
    if 'NGC869_B' in filters[i]:
        #Aperture and annulus rings of masterfiles        
        aperture = CircularAperture(NGC869_B_positions, r=apertureradius) 
        annulus_aperture = CircularAnnulus(NGC869_B_positions, r_in=annulusradius_in, r_out=annulusradius_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        ##bkg_mode = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data[i][25:4000, 25:4000])   
            annulus_data_1d = annulus_data[mask.data > 0] ##DUS HIER WEL DATA GEBRUIKEN IPV BKG_SUB?
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            ##mode = stats.mode(annulus_data_1d)
            ##bkg_mode.append(mode[1])
        bkg_median = np.array(bkg_median)
        ##bkg_mode = np.array(bkg_mode)
        
        max_aperture = []
        aperture_masks = aperture.to_mask(method='center')
        for mask in aperture_masks:
            aperture_data = mask.multiply(data[i][25:4000, 25:4000])
            max_aperture_value = np.amax(aperture_data)
            max_aperture.append(max_aperture_value)
        
        #Create phottable for every filter of object plane 
        phot869B = aperture_photometry(data[i][25:4000, 25:4000], aperture)
        phot869B['object'] = filters[i]
        filtertype = fits.getheader(alldata[i])['FILTER']
        phot869B['filtertype'] = filtertype
        exposuretime = fits.getheader(alldata[i])['EXPTIME']
        phot869B['exposuretime'] = exposuretime
        phot869B['max aperture'] = max_aperture
        phot869B['annulus_median'] = bkg_median
        phot869B['aperture_bkg'] = bkg_median * aperture.area
        phot869B['aperture_sum_bkgsub'] = phot869B['aperture_sum'] - phot869B['aperture_bkg']
        phot869B['SNR'] = phot869B['aperture_sum_bkgsub'] / (phot869B['aperture_sum_bkgsub'] + phot869B['aperture_bkg'])**0.5
        phot869B['Absolute error'] = phot869B['aperture_sum_bkgsub'] / phot869B['SNR']
        phot869B['mag inst'] = -2.5 * np.log10(np.abs(phot869B['aperture_sum_bkgsub']))
        phot869B['mag error'] = 2.5 * 0.434 * (1/phot869B['SNR']) #alleen afhankelijk van de SNR
        for col in phot869B.colnames:
            if not col == 'filtertype':
                if not col == 'exposuretime':
                    if not col == 'object':
                        phot869B[col].info.format = '%.6g'  #for consistent table output
        print(phot869B)
        print("")
        
        phot869B_allfilters.append(phot869B)

#Create a CSV file and add headers (corresponding to the phottable headers)           
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869B.csv', 'w') as file:
    writer = csv.writer(file)
    
    headerlist = phot869B.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()
    
    #Add phottables to CSV file
    for filter in phot869B_allfilters:
        phot869B_array = np.array(filter)
        writer.writerows(phot869B_array)

#Sort the csv file so every object is grouped with each a different filter & exposure time 
sortedNGC869B = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869B.csv')
sortedNGC869B.sort_values(['id'], axis=0, inplace=True)
sortedNGC869B.to_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC869B.csv', index=False)


###########################################NGC884A###########################################

#Apply photometry on all filters of every plane, using the location of objects found in the Rlong filter 
phot884A_allfilters = []
for i in range(len(filters)):
    if 'NGC884_A' in filters[i]:
        #Aperture and annulus rings of masterfiles
        aperture = CircularAperture(NGC884_A_positions, r=apertureradius) 
        annulus_aperture = CircularAnnulus(NGC884_A_positions, r_in=annulusradius_in, r_out=annulusradius_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        ##bkg_mode = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data[i][25:4000, 25:4000])   
            annulus_data_1d = annulus_data[mask.data > 0] ##DUS HIER WEL DATA GEBRUIKEN IPV BKG_SUB?
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            ##mode = stats.mode(annulus_data_1d)
            ##bkg_mode.append(mode[1])
        bkg_median = np.array(bkg_median)
        ##bkg_mode = np.array(bkg_mode)
        
        max_aperture = []
        aperture_masks = aperture.to_mask(method='center')
        for mask in aperture_masks:
            aperture_data = mask.multiply(data[i][25:4000, 25:4000])
            max_aperture_value = np.amax(aperture_data)
            max_aperture.append(max_aperture_value)

        #Create phottable for every filter of object plane 
        phot884A = aperture_photometry(data[i][25:4000, 25:4000], aperture)
        phot884A['object'] = filters[i]
        filtertype = fits.getheader(alldata[i])['FILTER']
        phot884A['filtertype'] = filtertype
        exposuretime = fits.getheader(alldata[i])['EXPTIME']
        phot884A['exposuretime'] = exposuretime
        phot884A['max aperture'] = max_aperture
        phot884A['annulus_median'] = bkg_median
        phot884A['aperture_bkg'] = bkg_median * aperture.area
        phot884A['aperture_sum_bkgsub'] = phot884A['aperture_sum'] - phot884A['aperture_bkg']
        phot884A['SNR'] = phot884A['aperture_sum_bkgsub'] / (phot884A['aperture_sum_bkgsub'] + phot884A['aperture_bkg'])**0.5
        phot884A['Absolute error'] = phot884A['aperture_sum_bkgsub'] / phot884A['SNR']
        phot884A['mag inst'] = -2.5 * np.log10(np.abs(phot884A['aperture_sum_bkgsub']))
        phot884A['mag error'] = 2.5 * 0.434 * (1/phot884A['SNR']) #alleen afhankelijk van de SNR
        for col in phot884A.colnames:
            if not col == 'filtertype':
                if not col == 'exposuretime':
                    if not col == 'object':
                        phot884A[col].info.format = '%.6g'  #for consistent table output
        id = phot884A['id'] 
        print(phot884A)
        print("")
        
        phot884A_allfilters.append(phot884A)

#Create a CSV file and add headers (corresponding to the phottable headers)           
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884A.csv', 'w') as file:
    writer = csv.writer(file)
    
    headerlist = phot884A.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()
    
    #Add phottables to CSV file
    for filter in phot884A_allfilters:
        phot884A_array = np.array(filter)
        writer.writerows(phot884A_array)

#Sort the csv file so every object is grouped with each a different filter & exposure time 
sortedNGC884A = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884A.csv')
sortedNGC884A.sort_values(['id'], axis=0, inplace=True)
sortedNGC884A.to_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884A.csv', index=False)


###########################################NGC884B###########################################

#Apply photometry on all filters of every plane, using the location of objects found in the Rlong filter 
phot884B_allfilters = []
for i in range(len(filters)):
    if 'NGC884_B' in filters[i]:
        #Aperture and annulus rings of masterfiles
        aperture = CircularAperture(NGC884_B_positions, r=apertureradius)
        annulus_aperture = CircularAnnulus(NGC884_B_positions, r_in=annulusradius_in, r_out=annulusradius_out)
        annulus_masks = annulus_aperture.to_mask(method='center')

        bkg_median = []
        ##bkg_mode = []
        for mask in annulus_masks:
            annulus_data = mask.multiply(data[i][25:4000, 25:4000])   
            annulus_data_1d = annulus_data[mask.data > 0] ##DUS HIER WEL DATA GEBRUIKEN IPV BKG_SUB?
            _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
            bkg_median.append(median_sigclip)
            ##mode = stats.mode(annulus_data_1d)
            ##bkg_mode.append(mode[1])
        bkg_median = np.array(bkg_median)
        ##bkg_mode = np.array(bkg_mode)
        
        max_aperture = []
        aperture_masks = aperture.to_mask(method='center')
        for mask in aperture_masks:
            aperture_data = mask.multiply(data[i][25:4000, 25:4000])
            max_aperture_value = np.amax(aperture_data)
            max_aperture.append(max_aperture_value)
            
        #Create phottable for every filter of object plane 
        phot884B = aperture_photometry(data[i][25:4000, 25:4000], aperture)
        phot884B['object'] = filters[i]
        filtertype = fits.getheader(alldata[i])['FILTER']
        phot884B['filtertype'] = filtertype
        exposuretime = fits.getheader(alldata[i])['EXPTIME']
        phot884B['exposuretime'] = exposuretime
        phot884B['max aperture'] = max_aperture
        phot884B['annulus_median'] = bkg_median
        phot884B['aperture_bkg'] = bkg_median * aperture.area
        phot884B['aperture_sum_bkgsub'] = phot884B['aperture_sum'] - phot884B['aperture_bkg']
        phot884B['SNR'] = phot884B['aperture_sum_bkgsub'] / (phot884B['aperture_sum_bkgsub'] + phot884B['aperture_bkg'])**0.5
        phot884B['Absolute error'] = phot884B['aperture_sum_bkgsub'] / phot884B['SNR']
        phot884B['mag inst'] = -2.5 * np.log10(np.abs(phot884B['aperture_sum_bkgsub']))
        phot884B['mag error'] = 2.5 * 0.434 * (1/phot884B['SNR']) #alleen afhankelijk van de SNR
        for col in phot884B.colnames:
            if not col == 'filtertype':
                if not col == 'exposuretime':
                    if not col == 'object':
                        phot884B[col].info.format = '%.6g'  #for consistent table output
        id = phot884B['id'] 
        print(phot884B)
        print("")
        
        phot884B_allfilters.append(phot884B)

#Create a CSV file and add headers (corresponding to the phottable headers)           
with open('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884B.csv', 'w') as file:
    writer = csv.writer(file)
    
    headerlist = phot884B.colnames
    headers = csv.DictWriter(file, delimiter = ',', fieldnames = headerlist)
    headers.writeheader()
    
    #Add phottables to CSV file
    for filter in phot884B_allfilters:
        phot884B_array = np.array(filter)
        writer.writerows(phot884B_array)

#Sort the csv file so every object is grouped with each a different filter & exposure time 
sortedNGC884B = pd.read_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884B.csv')
sortedNGC884B.sort_values(['id'], axis=0, inplace=True)
sortedNGC884B.to_csv('/Users/evagmelichmeijling/OneDrive/Capstone/BachelorProjectProgram/1_Starfinder/NGC884B.csv', index=False)

DOUBLECLUSTERTXT.write("\n")
DOUBLECLUSTERTXT.close()

print(datetime.now() - start)