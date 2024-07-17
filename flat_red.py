import numpy as np
import astropy.io.fits as fits
from tqdm import tqdm


def flatcomcpu(flat,biascom,file,objfilter):
    count=0
    nx = fits.getval(flat[0], 'NAXIS1')
    ny = fits.getval(flat[0], 'NAXIS2')
    flat_cube = np.empty((len(flat), ny, nx), dtype='float32')
    # 加载数据
    for i in range(len(flat)):
        flat_tmp = fits.getdata(flat[i])
        flat_tmp = flat_tmp - biascom
        flat_med = np.median(flat_tmp)
        flat_cube[i] = flat_tmp / flat_med
    flatcom = np.median(flat_cube, axis=0)

    
    #print(flatcom.shape)
    fits.writeto(f'{file}/Flat{objfilter}com.fit',data=flatcom,overwrite=True)
    return(flatcom)