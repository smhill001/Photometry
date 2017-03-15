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

import pylab as pl
import Meta_and_Control_Data_Operations as Meta
import SpecPhotPlot as SPP

#Plot="Photometry Project 2011"
#Plot="16-17-Dra"
#Plot="Vega"
#Plot="Antares"
#Plot="70-Oph"
#Plot="Alp-Her"
#Plot="BetLyr"
#Plot="S-ORI"
Plot="EpsLyr"
width=0.5
clrs=SPP.StarRainbow()

plotparams=Meta.spec_plot_params(drive,Plot)
#plotparams.DateTimeKey=Measurements.DateUT[1]
plotparams.DateTimeKey=""
print "plotparams.MeasList=",plotparams.ID,plotparams.MeasList
Measurements=Meta.measurement_list(plotparams.MeasList)
print "Measurements=",Measurements
print Measurements.DateUT
print "drive,Target=",drive,Plot
print "Measurements.FileList=",Measurements.FileList

first=True
for MeasIndex in range(0,Measurements.NObs):
    Centroid=[Measurements.Xcen[MeasIndex],Measurements.Ycen[MeasIndex]]                
    Radii=[Measurements.R1[MeasIndex],Measurements.R2[MeasIndex],Measurements.R3[MeasIndex]]
    PathName='f:/Astronomy/Projects/'+Measurements.DataType[MeasIndex]+'/'+Measurements.DataTarget[MeasIndex]+\
        '/Imaging Data/'+Measurements.DateUT[MeasIndex]+'/'
    print Measurements.FileList[MeasIndex]
    FNArray=Meta.GetObsFileNames(PathName,Measurements.FileList[MeasIndex])

    WavelengthCenters,NetCountsArray=SPP.BroadBandSpectrum(PathName,FNArray,Centroid,Radii)
    print "TEST:",Measurements.MeasTarget[MeasIndex][-2:]
    if Measurements.MeasTarget[MeasIndex][-3:]==" ND":
        print "Hi"
        NetCountsArray=NetCountsArray/0.04

    Label=Measurements.MeasTarget[MeasIndex]+'; '+Measurements.DataTarget[MeasIndex]+\
        '; '+Measurements.DateUT[MeasIndex]
    clr=clrs.c1[MeasIndex % 6,:]
    SPP.PlotBroadBand(WavelengthCenters,NetCountsArray,Label,clr,first,plotparams,width)
    first=False

pl.subplots_adjust(left=0.08, bottom=0.15, right=0.98, top=0.90,
            wspace=None, hspace=None)

pl.savefig(Plot+'.png',dpi=300)

###Testing Area

