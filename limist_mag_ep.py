import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import astropy.io.fits as fits
import os
import pandas as pd
import astropy.stats as astats
import scipy.stats as stats
from scipy.signal import savgol_filter 
from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d


import warnings
warnings.filterwarnings("ignore")

def run_cmd_Popen_fileno(cmd_string):
    import subprocess
    return subprocess.Popen(cmd_string, shell=True, stdout=None, stderr=None).wait()

def run_xmarch(file_path_cat,ref_cata,ra,dec):
    run_cmd_Popen_fileno("stilts cdsskymatch cdstable=%s \
                find=best in=%s \
                ra=%s dec=%s radius=10 out=%s-xm.csv"\
                %(ref_cata,file_path_cat,ra,dec,file_path_cat[:-4]))
    file_path_csv=file_path_cat[:-4]+'-xm.csv'
    return file_path_csv

def linear_func(x, a, b):
        return a * x + b
def func(x, a, b,c,d,e):
    return c*np.log(a * x) + b+np.e**(d*x)+e*x**4

def get_Zmag_limitmag(filename,starname,inst_mag,inst_magerr,refmag):
    exptime=fits.getval(filename,'EXPTIME')
    data=pd.read_csv(starname)
    data=data.sort_values(inst_mag)
    data_mag=np.array([data[inst_mag],data[inst_magerr],data[refmag]])
    max_mag=np.nanmax(data[inst_mag])
    min_mag=np.nanmin(data[inst_mag])
    mag_err_median=np.nanmedian(data[refmag]-data[inst_mag])
    mag_err_std=np.nanstd(data[refmag]-data[inst_mag])
    #print(len(data[refmag]))

    err_=np.arange(0.02,0.1,step=0.01)
    std_=np.arange(0.5,5,step=0.5)
    esaz_list=[]
    for i in err_:
        for j in std_:
            index=np.where((data_mag[0]>min_mag) & (data_mag[0]<max_mag) & (np.isfinite(data_mag[0])) & (np.isfinite(data_mag[2])) & \
                (abs(data_mag[2]-data_mag[0])<=mag_err_median+j*mag_err_std) & \
                (abs(data_mag[2]-data_mag[0])>=mag_err_median-j*mag_err_std) & (data_mag[1]<=i))
    #print(index)
            if len(index[0])>=3:
                params, _ = curve_fit(linear_func,data_mag[0][index],data_mag[2][index])
                a, Zmag = params

                esaz_list.append([i,j,abs(a-1),Zmag])

    esaz_arr=np.array(esaz_list).T
    Zmag=esaz_arr[3][np.nanargmin(esaz_arr[2])]
    mag_c=data_mag[0]+Zmag
    
    data_mag_new=np.array([mag_c,data_mag[1]])
    data_mag_new_sort=np.sort(data_mag_new)
    
    index_sig5=np.where((data_mag_new_sort[1]<0.22)&\
                       (data_mag_new_sort[1]>0.18))
    limit_mag_mean=np.nanmean(data_mag_new_sort[0][index_sig5])
    limit_mag_median=np.nanmedian(data_mag_new_sort[0][index_sig5])
    print(limit_mag_mean)
    if np.isnan(limit_mag_mean):
        popt, pcov = curve_fit(func,data_mag_new[1], data_mag_new[0],maxfev = 100000)
        limit_mag=func(0.2,*popt)
        limit_mag_mean=limit_mag
        limit_mag_median = limit_mag
    
    return filename,exptime,Zmag,limit_mag_mean,limit_mag_median



