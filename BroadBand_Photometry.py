# -*- coding: utf-8 -*-
"""
Created on Tue Jun 24 09:41:47 2014

A simple example of using ASTROPY to open and disply a FITS file of Jupiter
along with a color bar.

Not sure why, but when I ran this deep into a session, it gave errors the
first couple of times, then worked. Maybe not importing all I need to?

@author: Steven
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

Target="Vega"
DateObs="20110814UT"
width=0.5
clrs=SPL.StarRainbow()
Observations=SPL.observation_list('BroadBandPhotometry.txt')
print "Observations=",Observations
print Observations.DateUT
print "drive,Target=",drive,Target
plotparams=SPL.spec_plot_params(drive,Target)
plotparams.DateTimeKey=Observations.DateUT[1]

first=True
print Observations.FileList

for Obsindex in range(0,Observations.NObs):
    Centroid=[Observations.Xcen[Obsindex],Observations.Ycen[Obsindex]]                
    Radii=[Observations.R1[Obsindex],Observations.R2[Obsindex],Observations.R3[Obsindex]]
    PathName='f:/Astronomy/Projects/Photometry Project 2011/'+Observations.DateUT[Obsindex]+'/'
    FNArray=SPL.GetStarObsFileNames(PathName,Observations.FileList[Obsindex])

    ############################### Antares 20110809UT                
    
    WavelengthCenters,NetCountsArray=SPL.BroadBandSpectrum(PathName,FNArray,Centroid,Radii)

    Label=Observations.StarIdentifierDD[Obsindex]+' '+Observations.DateUT[Obsindex]
    clr=clrs.c1[Obsindex % 6,:]
    #PlotBroadBand(WavelengthCenters,NetCountsArray,LBL,clr,first,plotparams,width):
    SPL.PlotBroadBand(WavelengthCenters,NetCountsArray,Label,clr,first,plotparams,width)
    first=False


pl.savefig('Photometry2011.png',dpi=300)

###Testing Area

