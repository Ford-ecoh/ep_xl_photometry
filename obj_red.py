import numpy as np
import astropy.io.fits as fits
from tqdm import tqdm

def obj_flat_bias(obj,flatcom,biascom,filepath,red_filepath):
    for i in tqdm(range(len(obj))):
        obji_g=np.array(fits.open(obj[i])[0].data,dtype=np.float32)
        obji_gz=np.subtract(obji_g,biascom)
        obji_gzf=np.divide(obji_gz,flatcom)

        objheader=fits.open(obj[i])[0].header
        corhdr = objheader.copy()
        fits.writeto(f'{red_filepath}/{obj[i][len(filepath):-4]}zf.fit',\
                     data=obji_gzf,header=corhdr,overwrite=True)
