<<<<<<< HEAD
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 09:50:38 2017

@author: Astronomy

New baseline version of SPL Library 3/5/2017

"""


class spec_plot_params:
    def __init__(self,drive,ObjIdentifierDD):
        print "In spec_plot_params"
        #View has two options: raw or flux?
        self.ID=ObjIdentifierDD
        self.TargetType=''  #PN, ST, etc.
        self.X0=0.
        self.X1=0.
        self.NX=0.
        self.Y0=0.
        self.Y1=0.
        self.NY=0.
        self.Ytype=''  #log or linear
        self.ObsType=''   #broadband or spectrum
        self.ModelFile='' #spectral model file.observation_list

        CfgFile=open(drive+'/Astronomy/Python Play/SpectroPhotometryLibraries/SpecPlotConfig.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            #print fields[0], fields[1]
            if fields[0] == ObjIdentifierDD:
                #print "In first if, fields[1]",fields[1],fields[0]
                self.TargetType=str(fields[1])
                self.X0=float(fields[2])
                self.X1=float(fields[3])
                self.NX=float(fields[4])
                self.Y0=float(fields[5])
                self.Y1=float(fields[6])
                self.NY=float(fields[7])
                self.Ytype=str(fields[8])
                self.ObsType=str(fields[9])

class StarRainbow:
    def __init__(self):
        import numpy as np
        self.c1=np.array([[0.7,0.,0.],[0.7,0.7,0.0],[0.,0.7,0.], \
                             [0.7,0.,0.7],[0.,0.7,0.7],[0.,0.,0.7]])
                             
class observation_list:
    def __init__(self,ObservationListFile):
        #The initial plan is to read ALL records in the observation list
        self.StarIdentifierDD=[''] #Keyword for star identification
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Instrument=['']       #Instrument code, to be used for aperture
        self.FileList=['']         #List of observation image files (FITS)
        self.Xcen=[0.]             #X center coordinate in FITS image (pixels)
        self.Ycen=[0.]             #Y center coordinate in FITS image (pixels)
        self.R1=[0.]               #Radius of photometry aperture (pixels)
        self.R2=[0.]               #Inner radius of background annulus (pixels)
        self.R3=[0.]               #Outer radius of background annulus (pixels)
        self.NObs=0                #Number of observatinos
        self.FirstTime=True
        
        CfgFile=open(ObservationListFile,'r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            if self.FirstTime:
                self.StarIdentifierDD[0]=str(fields[0])
                self.DateUT[0]=str(fields[1])
                self.Instrument[0]=str(fields[2])
                self.FileList[0]=str(fields[3])
                self.Xcen[0]=float(fields[4])
                self.Ycen[0]=float(fields[5])
                self.R1[0]=float(fields[6])
                self.R2[0]=float(fields[7])
                self.R3[0]=float(fields[8])
                self.FirstTime=False
                self.NObs=1
            else:
                self.StarIdentifierDD.extend([str(fields[0])])
                self.DateUT.extend([str(fields[1])])
                self.Instrument.extend([str(fields[2])])
                self.FileList.extend([str(fields[3])])
                self.Xcen.extend([float(fields[4])])
                self.Ycen.extend([float(fields[5])])
                self.R1.extend([float(fields[6])])
                self.R2.extend([float(fields[7])])
                self.R3.extend([float(fields[8])])
                self.NObs=self.NObs+1

def GetStarObsFileNames(Path,IndexFile):
    """
    Used by *BroadBand_Photometry_V1.py* to obtain a list of FITS files
    contributing to a single broadband spectrum. Each FITS file
    self-identifies it's filter from the header meta-data and the 
    parameters of that filter with regard to photometry are then
    obtained by a table look up routine (list here)
    """
    CfgFile=open(IndexFile,'r')
    CfgLines=CfgFile.readlines()
    CfgFile.close()
    nrecords=len(CfgLines)
    #print CfgLines
    FNArray=['']
    FirstTime=True
    for recordindex in range(0,nrecords):
        fields=CfgLines[recordindex].split(',')
        if FirstTime:
            FNArray[0]=str(fields[0])
            FirstTime=False
        else:
            FNArray.extend([str(fields[0])])
            
    return FNArray

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
    import PhotometryLib as phot
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=phot.FilterParameters(hdulist[0].header['FILTER'])
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
    import PhotometryLib as phot
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=phot.FilterParameters(hdulist[0].header['FILTER'])
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
    #import pylab as pl
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
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=0,label=LBL,color=clr)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=1, borderaxespad=0.,prop={'size':4})    
    else:
        if plotparams.ObsType=='broadband':
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=0,label=LBL,color=clr)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=1, borderaxespad=0.,prop={'size':5})
                             
=======
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 09:50:38 2017

@author: Astronomy

New baseline version of SPL Library 3/5/2017

"""


class spec_plot_params:
    def __init__(self,drive,ObjIdentifierDD):
        print "In spec_plot_params"
        #View has two options: raw or flux?
        self.ID=ObjIdentifierDD
        self.TargetType=''  #PN, ST, etc.
        self.X0=0.
        self.X1=0.
        self.NX=0.
        self.Y0=0.
        self.Y1=0.
        self.NY=0.
        self.Ytype=''  #log or linear
        self.ObsType=''   #broadband or spectrum
        self.ModelFile='' #spectral model file.observation_list

        CfgFile=open(drive+'/Astronomy/Python Play/SpectroPhotometryLibraries/SpecPlotConfig.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            #print fields[0], fields[1]
            if fields[0] == ObjIdentifierDD:
                #print "In first if, fields[1]",fields[1],fields[0]
                self.TargetType=str(fields[1])
                self.X0=float(fields[2])
                self.X1=float(fields[3])
                self.NX=float(fields[4])
                self.Y0=float(fields[5])
                self.Y1=float(fields[6])
                self.NY=float(fields[7])
                self.Ytype=str(fields[8])
                self.ObsType=str(fields[9])

class StarRainbow:
    def __init__(self):
        import numpy as np
        self.c1=np.array([[0.7,0.,0.],[0.7,0.7,0.0],[0.,0.7,0.], \
                             [0.7,0.,0.7],[0.,0.7,0.7],[0.,0.,0.7]])
                             
class observation_list:
    def __init__(self,ObservationListFile):
        #The initial plan is to read ALL records in the observation list
        self.StarIdentifierDD=[''] #Keyword for star identification
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Instrument=['']       #Instrument code, to be used for aperture
        self.FileList=['']         #List of observation image files (FITS)
        self.Xcen=[0.]             #X center coordinate in FITS image (pixels)
        self.Ycen=[0.]             #Y center coordinate in FITS image (pixels)
        self.R1=[0.]               #Radius of photometry aperture (pixels)
        self.R2=[0.]               #Inner radius of background annulus (pixels)
        self.R3=[0.]               #Outer radius of background annulus (pixels)
        self.NObs=0                #Number of observatinos
        self.FirstTime=True
        
        CfgFile=open(ObservationListFile,'r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            if self.FirstTime:
                self.StarIdentifierDD[0]=str(fields[0])
                self.DateUT[0]=str(fields[1])
                self.Instrument[0]=str(fields[2])
                self.FileList[0]=str(fields[3])
                self.Xcen[0]=float(fields[4])
                self.Ycen[0]=float(fields[5])
                self.R1[0]=float(fields[6])
                self.R2[0]=float(fields[7])
                self.R3[0]=float(fields[8])
                self.FirstTime=False
                self.NObs=1
            else:
                self.StarIdentifierDD.extend([str(fields[0])])
                self.DateUT.extend([str(fields[1])])
                self.Instrument.extend([str(fields[2])])
                self.FileList.extend([str(fields[3])])
                self.Xcen.extend([float(fields[4])])
                self.Ycen.extend([float(fields[5])])
                self.R1.extend([float(fields[6])])
                self.R2.extend([float(fields[7])])
                self.R3.extend([float(fields[8])])
                self.NObs=self.NObs+1

def GetStarObsFileNames(Path,IndexFile):
    """
    Used by *BroadBand_Photometry_V1.py* to obtain a list of FITS files
    contributing to a single broadband spectrum. Each FITS file
    self-identifies it's filter from the header meta-data and the 
    parameters of that filter with regard to photometry are then
    obtained by a table look up routine (list here)
    """
    CfgFile=open(IndexFile,'r')
    CfgLines=CfgFile.readlines()
    CfgFile.close()
    nrecords=len(CfgLines)
    #print CfgLines
    FNArray=['']
    FirstTime=True
    for recordindex in range(0,nrecords):
        fields=CfgLines[recordindex].split(',')
        if FirstTime:
            FNArray[0]=str(fields[0])
            FirstTime=False
        else:
            FNArray.extend([str(fields[0])])
            
    return FNArray

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
    import PhotometryLib as phot
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=phot.FilterParameters(hdulist[0].header['FILTER'])
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
    import PhotometryLib as phot
    from astropy.table import Table, hstack
    from astropy.io import fits

    hdulist=fits.open(FN)
    Filter=phot.FilterParameters(hdulist[0].header['FILTER'])
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
    #import pylab as pl
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
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=0,label=LBL,color=clr)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=1, borderaxespad=0.,prop={'size':4})    
    else:
        if plotparams.ObsType=='broadband':
            pl.scatter(WavelengthCenters,NetCountsArray,linewidth=0,label=LBL,color=clr)
        elif plotparams.ObsType=='spectrum':
            pl.plot(WavelengthCenters,NetCountsArray,linewidth=width,label=LBL,color=clr)
        pl.legend(loc=1,ncol=1, borderaxespad=0.,prop={'size':5})
                             
>>>>>>> origin/master
