# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 09:50:38 2017

MODULES:    StarRainbow
            BroadBandSpectrum
            ComputeNetRate
            ComputeNetRate1
            PlotBroadBand
            
@author: Astronomy

New baseline version of SPL Library 3/5/2017

"""


class StarRainbow:
    def __init__(self):
        import numpy as np
        self.c1=np.array([[0.7,0.,0.],[0.7,0.7,0.0],[0.,0.7,0.], \
                             [0.7,0.,0.7],[0.,0.7,0.7],[0.,0.,0.7]])
                             

def BroadBandSpectrum(PathName,FNArray,Centroid,Radii):
    import numpy as np
    ArrayLen=len(FNArray)
    WavelengthCenters=np.zeros(ArrayLen)
    NetCountsArray=np.zeros(ArrayLen)
    
    for FNindex in range(0,ArrayLen):
        
        rate,WVCenter=ComputeNetRate(PathName+FNArray[FNindex],
                                     Centroid,Radii)
        WavelengthCenters[FNindex]=WVCenter
        NetCountsArray[FNindex]=rate
    return WavelengthCenters,NetCountsArray

def ComputeNetRate(FN,positions,radii):
    from photutils import CircularAperture
    from photutils import aperture_photometry
    from photutils import CircularAnnulus
    import Meta_and_Control_Data_Operations as Meta
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=Meta.FilterParameters(hdulist[0].header['FILTER'])
    scidata=hdulist[0].data

    #Create aperture objects    
    apertures = CircularAperture(positions, r=radii[0])
    annulus_apertures = CircularAnnulus(positions, r_in=radii[1], r_out=radii[2])

    #Compute raw fluxes
    rawflux_table = aperture_photometry(scidata, apertures)
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()

    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_raw'] - bkg_sum
    rate=((final_sum/hdulist[0].header['EXPTIME'])/Filter.Aperture)/Filter.EW
    phot_table['net_count_rate'] = rate

    WVCenter=Filter.CenterWV
    
    return rate,WVCenter

def ComputeNetRate1(FN,positions,radii):
    """
    Alternative test version of ComputeNetRate. This version is based on
    the image segmentation and morphology capabilities in photutils instead
    of on simple aperture photometry.  Ultimately, it may combine both
    
    SMH 6/6/2016
    """
    from photutils import CircularAperture
    from photutils import aperture_photometry
    from photutils import CircularAnnulus
    import Meta_and_Control_Data_Operations as Meta
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=Meta.FilterParameters(hdulist[0].header['FILTER'])
    scidata=hdulist[0].data

    #Create aperture objects    
    apertures = CircularAperture(positions, r=radii[0])
    annulus_apertures = CircularAnnulus(positions, r_in=radii[1], r_out=radii[2])

    #Compute raw fluxes
    rawflux_table = aperture_photometry(scidata, apertures)
    bkgflux_table = aperture_photometry(scidata, annulus_apertures)
    
    phot_table = hstack([rawflux_table, bkgflux_table], table_names=['raw', 'bkg'])
    bkg_mean = phot_table['aperture_sum_bkg'] / annulus_apertures.area()

    bkg_sum = bkg_mean * apertures.area()
    final_sum = phot_table['aperture_sum_raw'] - bkg_sum
    rate=((final_sum/hdulist[0].header['EXPTIME'])/Filter.Aperture)/Filter.EW
    phot_table['net_count_rate'] = rate

    WVCenter=Filter.CenterWV
    
    return rate,WVCenter

def PlotBroadBand(WavelengthCenters,NetCountsArray,LBL,clr,first,plotparams,width):
    import pylab as pl
    import numpy as np
    print "Hi. In PlotBroadBand", first
    if first==True:
        pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")
        pl.subplot(1, 1, 1)
        #Plot Layout Configuration

        # Set x limits
        pl.xlim(plotparams.X0,plotparams.X1)
        # Set x ticks
        pl.xticks(np.linspace(plotparams.X0,plotparams.X1,plotparams.NX, endpoint=True))
        # Set y limits
        pl.ylim(plotparams.Y0,plotparams.Y1)
        # Set y ticks
        if plotparams.Ytype=='log':
            pl.yscale(plotparams.Ytype)
        
        pl.grid()
        pl.tick_params(axis='both', which='major', labelsize=7)
        pl.ylabel(r"$Counts-s^{-1}$-$m^{-2}$-$nm^{-1}$",fontsize=7)
        pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
        
        pl.title(plotparams.ID+' '+plotparams.DateTimeKey,fontsize=9) #Should change this to include objectID and common name
        #print plotparams.Type
        if plotparams.ObsType=='broadband':
            print "in branch"
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=1,label=LBL,color=clr,marker='o',s=2)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':4})    
    else:
        if plotparams.ObsType=='broadband':
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=1,label=LBL,color=clr,marker='o',s=2)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=3, borderaxespad=0.,prop={'size':5})
                             
