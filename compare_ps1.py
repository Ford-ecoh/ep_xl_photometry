from ps1_single import *
from matplotlib.patches import Circle
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord 
from astropy.coordinates import ICRS, Galactic, FK4, FK5 
import astropy.units as u
import matplotlib.pyplot as plt
import astropy.stats as astats
import pandas as pd

def show_fits(fitsname,torget_name,red_filepath,outpng,Zmag):

    data=fits.open(f'{fitsname}.fit')[0].data
    head=fits.open(f'{fitsname}.fit')[0].header

    img_x=fits.getval(f'{fitsname}.fit','NAXIS1')
    img_y=fits.getval(f'{fitsname}.fit','NAXIS2')

    ep_xm=pd.read_csv(f'{fitsname}-xm.csv')
    data_mag=np.array([ep_xm['NUMBER'],\
                  ep_xm['X_IMAGE'],ep_xm['Y_IMAGE'],\
                  ep_xm['X_WORLD'],ep_xm['Y_WORLD'],\
                  ep_xm['MAG_BEST']+Zmag,ep_xm['MAGERR_BEST'],ep_xm['gmag'],])
    max_mag=np.nanmax(ep_xm['MAG_BEST']+Zmag)
    min_mag=np.nanmin(ep_xm['MAG_BEST']+Zmag)
    mag_err_median=np.nanmedian(ep_xm['gmag']-ep_xm['MAG_BEST']-Zmag)
    mag_err_std=np.nanstd(ep_xm['gmag']-ep_xm['MAG_BEST']-Zmag)

    index_may_1=np.where((data_mag[6]<=0.04)&((abs(data_mag[7]-data_mag[5])>=mag_err_median+2*mag_err_std) | \
                  (abs(data_mag[7]-data_mag[5])<=mag_err_median-2*mag_err_std)) \
                & (abs(img_y-data_mag[1])>100) & (abs(img_x-data_mag[2])>100)\
                & (data_mag[1]>100) & (data_mag[2]>100))


    im_cli=astats.sigma_clip(data,masked=False)
    im_med=np.median(im_cli)
    im_std=np.std(im_cli)
    
    fig = plt.figure(figsize=(10,10))
    plt.set_cmap('gray')
    plt.imshow(data,vmin=im_med-5*im_std,vmax=im_med+5*im_std)
    plt.axis('off')
    for i, x, y  in np.array([data_mag[0][index_may_1[0]],data_mag[1][index_may_1[0]],data_mag[2][index_may_1[0]]]).T:
        plt.scatter(x, y, color='red', s=5)  # Draw a red circle
        plt.text(x+20, y+20, str(int(i)),\
                 color='white', fontsize=12, ha='right', va='top')
    plt.savefig(f'{red_filepath}/{outpng}_maybe_VS.pdf')
    plt.close()

    col_name=['NUMBER','X_IMAGE',"Y_IMAGE","RA","DEC",'MAG_INST','MAGERR_INST','MAG_PS']
    lzfl=pd.DataFrame(data=data_mag[:,index_may_1[0]].T,columns=col_name)
    lzfl.to_csv(f'{red_filepath}/{outpng}_maybe_VS.csv')

    Ra=fits.getval(f'{fitsname}.fit','RA')
    Dec=fits.getval(f'{fitsname}.fit','DEC')
    size=fits.getval(f'{fitsname}.fit','NAXIS1')
    w = WCS(head)
    coords = [f"{Ra} {Dec}"]

    c = SkyCoord(coords, frame=ICRS, unit=(u.hour, u.deg), obstime="J2000")

    gfit,ghead= ps1_data(c.ra.deg[0], c.dec.deg[0], size, 'g')
    rfit,rhead= ps1_data(c.ra.deg[0], c.dec.deg[0], size, 'r')
    ifit,ihead= ps1_data(c.ra.deg[0], c.dec.deg[0], size, 'i')

    x, y = w.world_to_pixel(c)
    circle = Circle((x, y), 500, edgecolor='r', facecolor='none') 

    
    im_cli=astats.sigma_clip(data,masked=False)
    im_med=np.median(im_cli)
    im_std=np.std(im_cli)
    
    fig, ([ax1,ax2],[ax3,ax4]) = plt.subplots(2, 2,figsize=(50,50))
    fig.suptitle(outpng)
    #ax1=plt.subplot(2,2,1,projection=WCS(head))
    ax1.set_title(torget_name,fontsize=25)
    ax1.imshow(data,vmin=im_med-5*im_std,vmax=im_med+5*im_std,cmap='gray')
    #ax1.coords['ra'].set_axislabel('Right Ascension')
    #ax1.coords['dec'].set_axislabel('Declination')
    #ax1.add_patch(circle)
    ax1.set_axis_off()


    #ax2=plt.subplot(2,2,2,projection=WCS(ghead))
    ax2.set_title(f'PanStar g',fontsize=25)
    ax2.imshow(gfit,cmap="gray",origin="lower")
    ax2.set_axis_off()

    #ax3=plt.subplot(2,2,3,projection=WCS(rhead))
    ax3.set_title(f'PanStar r',fontsize=25)
    ax3.imshow(rfit,cmap="gray",origin="lower")
    ax3.set_axis_off()

    #ax4=plt.subplot(2,2,4,projection=WCS(ihead))
    ax4.set_title(f'PanStar i',fontsize=25)
    ax4.imshow(ifit,cmap="gray",origin="lower")
    ax4.set_axis_off()

    plt.savefig(f'{red_filepath}/{outpng}.pdf')
    plt.close()

