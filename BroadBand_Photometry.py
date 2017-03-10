# -*- coding: utf-8 -*-
"""
Updated on Tue Mar 07 09:41:47 2017
###############################################################################
NAME:       BroadBand_Photometry.py

PURPOSE:    To extract stellar photometric information from FITS image files
            and plot the results
            
INPUTS:     A single parameter, "Target", points to a configuration file
            that provides the stellar targets and plotting information. That
            configuration file also points to secondary configuraiton files
            that contain the FITS data file lists for each observation.
            
LIBRARIES:  This code calls the SpecPhotLibNew.py library. It is an updated
            subset of the SpecPhotLibV006.py library that had grown cumbersome.
                    

###############################################################################
@author: Steven Hill
"""
import sys
drive='f:'
sys.path.append(drive+'/Astronomy/Python Play')
sys.path.append(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry')
sys.path.append(drive+'/Astronomy/Python Play/FITSImageStuff')

#import numpy as np
#from astropy.io import fits
#from photutils import CircularAperture
#from photutils import aperture_photometry
#from photutils import CircularAnnulus
#from astropy.table import Table, hstack
import pylab as pl
import SpecPhotLibNew as SPL

#Target="Photometry Project 2011"
#Target="16-17-Dra-Photometry"
#Target="Vega Photometry"
#Target="Antares"
Target="70-Oph"
#DateObs="20110814UT"
width=0.5
clrs=SPL.StarRainbow()
#Observations=SPL.observation_list('Photometry Project 2011.txt')
Observations=SPL.observation_list(Target+'.txt')
print "Observations=",Observations
print Observations.DateUT
print "drive,Target=",drive,Target
plotparams=SPL.spec_plot_params(drive,Target)
#plotparams.DateTimeKey=Observations.DateUT[1]
plotparams.DateTimeKey=""

first=True
print Observations.FileList

for Obsindex in range(0,Observations.NObs):
    Centroid=[Observations.Xcen[Obsindex],Observations.Ycen[Obsindex]]                
    Radii=[Observations.R1[Obsindex],Observations.R2[Obsindex],Observations.R3[Obsindex]]
    PathName='f:/Astronomy/Projects/Stars/'+Observations.StarIdentifier[Obsindex]+\
        '/Imaging Data/'+Observations.DateUT[Obsindex]+'/'
    FNArray=SPL.GetStarObsFileNames(PathName,Observations.FileList[Obsindex])

    WavelengthCenters,NetCountsArray=SPL.BroadBandSpectrum(PathName,FNArray,Centroid,Radii)
    print "TEST:",Observations.Target[Obsindex][-2:]
    if Observations.Target[Obsindex][-3:]==" ND":
        print "Hi"
        NetCountsArray=NetCountsArray/0.04

    Label=Observations.StarIdentifier[Obsindex]+Observations.Target[Obsindex]+\
        ' '+Observations.DateUT[Obsindex]
    clr=clrs.c1[Obsindex % 6,:]
    SPL.PlotBroadBand(WavelengthCenters,NetCountsArray,Label,clr,first,plotparams,width)
    first=False

pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.90,
            wspace=None, hspace=None)

pl.savefig(Target+'.png',dpi=300)

###Testing Area

