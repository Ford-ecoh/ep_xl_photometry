import numpy as np
import astropy.io.fits as fits
from tqdm import tqdm


def biascomcpu(bias,cambineway,red_filepath):
    nx = fits.getval(bias[0], 'NAXIS1')
    ny = fits.getval(bias[0], 'NAXIS2')
    bias_cube = np.empty((len(bias), ny, nx), dtype='float32')
    # 加载数据
    for i in range(len(bias)):
        bias_cube[i] = fits.getdata(bias[i])
    if cambineway=="median":
        biascom = np.float32(np.median(bias_cube, axis=0))
    elif cambineway=="mean":
        biascom = np.float32(np.mean(bias_cube, axis=0))
    else:
        print("NO way")
    
    print(biascom.shape)

    fits.writeto(f'%s/Bias_com.fit'%(red_filepath),data=biascom,overwrite=True)
    return biascom