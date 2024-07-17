from ccdproc import Combiner
from ccdproc import wcs_project
import astropy.units as u
import numpy as np
from astropy.nddata import CCDData
from astropy.wcs import WCS
from reproject import reproject_interp
from astropy.table import Table
import astropy.io.fits as fits
import os
from astropy.wcs import WCS, Sip
import astropy.stats as astats
import astroalign as aa

import warnings
warnings.filterwarnings("ignore")


def combine_kernal(file):
    image1=np.array(fits.getdata(file[0]),dtype='float32')
    header=fits.open(file[0])[0].header
    for i in range(len(file)-1):
        image2=np.array(fits.getdata(file[i+1]),dtype='float32')
        registered_image, footprint = aa.register(image2, image1)
        ccd_1=CCDData(image1,unit=u.adu)
        ccd_2=CCDData(registered_image,unit=u.adu)
        combiner=Combiner([ccd_1,ccd_2])
        combiner.sigma_clipping(low_thresh=3,high_thresh=3,func=np.ma.median)
        image1=combiner.median_combine()
    return image1,header
def combine_fits_withoutwcs_auto(listfile,num,output):
    list_fits=[f.strip() for f in open(listfile).readlines() if f.strip()]
    count=0
    while count<len(list_fits):
        #print(list_fits[count:count+num])
        combine_ma,header=combine_kernal(list_fits[count:count+num])
        header['COMfitsNUM']=num
        fits.writeto(output+'-%04d.fits'%(count+1),\
                     data=combine_ma,header=header,overwrite=True)
        count+=num
        print(count)
    if len(list_fits)-count>0:
        combine_ma,header=combine_kernal(list_fits[count:len(list_fits)])
        header['COMfitsNUM']=len(list_fits)-count
        fits.writeto(output+'-%04d.fits'%(count+1),\
                     data=combine_ma,header=header,overwrite=True)

def combine_fits_withoutwcs(listfile,red_filepath,output):
    image1=np.array(fits.getdata(listfile[0]),dtype='float32')
    exptime=fits.getval(listfile[0],'EXPTIME')
    image1=image1/exptime
    header=fits.open(listfile[0])[0].header
    header['COMfitsNUM']=len(listfile)
    ccd_1=CCDData(image1,unit=u.adu)
    tempcom=[]
    tempcom.append(ccd_1)
    for i in range(len(listfile)-1):
        image2=np.array(fits.getdata(listfile[i+1]),dtype='float32')
        image2=image2/exptime
        registered_image, footprint = aa.register(image2, image1)
        ccd_2=CCDData(registered_image,unit=u.adu)
        tempcom.append(ccd_2)
        
    combiner=Combiner(tempcom)
    combiner.sigma_clipping(low_thresh=3,high_thresh=3,func=np.ma.median)
    image1=combiner.median_combine()
    
    fits.writeto(f'{red_filepath}/{output}',data=image1,header=header,overwrite=True)