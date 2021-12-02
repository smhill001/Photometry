# -*- coding: utf-8 -*-
"""
Created on Sun Mar 05 09:50:38 2017

MODULES:    FilterParameters
            spec_plot_params
            measurement_list
            GetStarObsFileNames
            
@author: Astronomy

New baseline version of SPL Library 3/5/2017

"""

class FilterParameters:
    def __init__(self,FilterName):
        import numpy as np
        #Extracts the EW parameters for a specific absorption band
        #from an EW file for a specific observation (target, time, 
        #processing level).
        
        #NEED TO CREATE A CLASS THAT WILL READ IN ALL THE EWS FOR
        #ALL BANDS IN AN EW FILE AND HANDLE THE SELECTION OF INDIVIDUAL
        #BAND PARAMETERS IN A SEPARATE METHOD OR CLASS.
        
        self.CenterWV=0.
        self.EW=0.
        self.Aperture=0.

        #FilterNames=["Blue","Green","Green","NIR","NUV","Red"]
        #WavelengthCenters=[450,550.,550.,750.,380.,650.]
        #WavelengthWidths=[110.,95.4,95.4,130.,40.,62.8]
        FilterNames=["380NUV","450BLU","486HIB","501OIII","550GRN","650RED","632OI",
                     "647CNT","656HIA","658NII","672SII","685NIR","730OII","742NIR",
                     "807NIR","889CH4","940NIR","1000NIR"]
        WavelengthCenters=[379.5,460.0,486.0,499.0,525.0,647.0,632.0,
                           647.0,656.0,658.0,672.0,842.5,730.0,871.0,
                           903.5,889.0,940.0,1000.0]
        WavelengthWidths=[15.4,100.0,10.0,10.0,90.0,75.0,10.0,
                          10.0,11.0,10.0,10.0,315.0,258.0,10.0,
                          193.0,11.2,10.0,30.0]
        #Widths are from integration bounds in "FluxCalibration.xlsx" for 
        #broadband filters. FWHM from analysis for NUV, HIA and CH4.  OIII,
        #HIB and SII are 10nm placeholders.
                          
        FilterIndex = [k for k, x in enumerate(FilterNames) if x == FilterName] #what does this do!?
        FI=np.int(FilterIndex[0])
        self.CenterWV=WavelengthCenters[FI]
        self.EW=WavelengthWidths[FI]
        self.Aperture= 0.2**2.-0.07**2. #meters^2
        
class spec_plot_params:
    def __init__(self,drive,PlotIdentifier):
        print "In spec_plot_params"
        #View has two options: raw or flux?
        self.ID=PlotIdentifier
        self.X0=0.
        self.X1=0.
        self.NX=0.
        self.Y0=0.
        self.Y1=0.
        self.NY=0.
        self.Ytype=''  #log or linear
        self.ObsType=''   #broadband or spectrum
        self.MeasList='' #spectral model file.measurement_list

        CfgFile=open(drive+'/Astronomy/Python Play/SpectroPhotometry/Photometry/PlotConfig.txt','r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            #print fields[0], fields[1]
            if fields[0] == PlotIdentifier:
                #print "In first if, fields[1]",fields[1],fields[0]
                self.X0=float(fields[1])
                self.X1=float(fields[2])
                self.NX=float(fields[3])
                self.Y0=float(fields[4])
                self.Y1=float(fields[5])
                self.NY=float(fields[6])
                self.Ytype=str(fields[7])
                self.ObsType=str(fields[8])
                self.MeasList=str(fields[9])

                             
class measurement_list:
    def __init__(self,MeasurementListFile):
        #The initial plan is to read ALL records in the observation list
        self.MeasTarget=['']   #Keyword for star identification
        self.DataType=['']           #Target, e.g., component of a multiple star
        self.DataTarget=['']           #Target, e.g., component of a multiple star
        self.DateUT=['']           #UT Date of observation: YYYYMMDDUT
        self.Optics=['']       #Instrument code, to be used for aperture
        self.Camera=['']       #Instrument code, to be used for aperture
        self.FileList=['']         #List of observation image files (FITS)
        self.Xcen=[0.]             #X center coordinate in FITS image (pixels)
        self.Ycen=[0.]             #Y center coordinate in FITS image (pixels)
        self.R1=[0.]               #Radius of photometry aperture (pixels)
        self.R2=[0.]               #Inner radius of background annulus (pixels)
        self.R3=[0.]               #Outer radius of background annulus (pixels)
        self.NObs=0                #Number of observatinos
        self.FirstTime=True
        
        CfgFile=open(MeasurementListFile,'r')
        CfgLines=CfgFile.readlines()
        CfgFile.close()
        nrecords=len(CfgLines)
        #print CfgLines

        for recordindex in range(1,nrecords):
            fields=CfgLines[recordindex].split(',')
            if self.FirstTime:
                self.MeasTarget[0]=str(fields[0])
                self.DataType[0]=str(fields[1])
                self.DataTarget[0]=str(fields[2])
                self.DateUT[0]=str(fields[3])
                self.Optics[0]=str(fields[4])
                self.Camera[0]=str(fields[5])
                self.FileList[0]=str(fields[6])
                self.Xcen[0]=float(fields[7])
                self.Ycen[0]=float(fields[8])
                self.R1[0]=float(fields[9])
                self.R2[0]=float(fields[10])
                self.R3[0]=float(fields[11])
                self.FirstTime=False
                self.NObs=1
            else:
                self.MeasTarget.extend([str(fields[0])])
                self.DataType.extend([str(fields[1])])
                self.DataTarget.extend([str(fields[2])])
                self.DateUT.extend([str(fields[3])])
                self.Optics.extend([str(fields[4])])
                self.Camera.extend([str(fields[5])])
                self.FileList.extend([str(fields[6])])
                self.Xcen.extend([float(fields[7])])
                self.Ycen.extend([float(fields[8])])
                self.R1.extend([float(fields[9])])
                self.R2.extend([float(fields[10])])
                self.R3.extend([float(fields[11])])
                self.NObs=self.NObs+1

def GetObsFileNames(Path,IndexFile):
    """
    Used by *BroadBand_Photometry_V1.py* to obtain a list of FITS files
    contributing to a single broadband spectrum. Each FITS file
    self-identifies it's filter from the header meta-data and the 
    parameters of that filter with regard to photometry are then
    obtained by a table look up routine (list here)
    """
    print "GetStarObsFileNames.IndexFile=",IndexFile
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

