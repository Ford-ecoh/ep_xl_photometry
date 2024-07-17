from bias_red import *
from flat_red import *
from obj_red import *
from wcs_sex import *
from limist_mag_ep import *
from ps1_single import *
import os
from tqdm import tqdm
import astropy.io.fits as fits
from combine import *
from compare_ps1 import *



def run_cmd_Popen_fileno(cmd_string):
    import subprocess
    return subprocess.Popen(cmd_string, shell=True, stdout=None, stderr=None).wait()



#obj_cooder

def ep_proresss(filepath,torget_name,filter):

    file_name = os.path.basename(filepath)
    file_path=os.path.dirname(filepath)
    red_filepath=filepath+"_red"

    run_cmd_Popen_fileno('mkdir %s'%(red_filepath))

    run_cmd_Popen_fileno('ls %s/*fit* > %s/all_list'%(filepath,red_filepath))
    file_list=[f.strip() for f in open("%s/all_list"%(red_filepath)).readlines() if f.strip()]

    bias_word = "bias"
    bias=[]
    flat_word = "flat"
    flat=[]
    for i in file_list:
        if bias_word.lower() in i.lower():
            bias.append(i)
        elif flat_word.lower() in i.lower():
            flat.append(i)



    print('Bias Begin')
    biascom=biascomcpu(bias,'median',red_filepath)
    print('Bias Done')



    print('Flat Begin')
    flat_f=[]
    flatcom_i=[]
    for i in filter:
        for j in flat:
            if i.lower() in j.lower():
                flat_f.append(j)

        flatcom_temp=flatcomcpu(flat_f,biascom,red_filepath,i)
        flatcom_i.append(flatcom_temp)
        flat_f.clear()
    print('Flat done')



    #obj
    flj=0
    for j in filter:
        run_cmd_Popen_fileno(f'ls {filepath}/*{torget_name}*{j}* > {red_filepath}/{torget_name}_{j}')
        obj_f=[f.strip() for f in open(f'{red_filepath}/{torget_name}_{j}').readlines() if f.strip()]
        if len(obj_f)!=0:
            obj_flat_bias(obj_f,flatcom_i[flj],biascom,filepath,red_filepath)
        flj+=1
    print(f'{torget_name} done')


    list_zmag_fwhm_lim=[]
    for j in filter:
        run_cmd_Popen_fileno(f'ls {red_filepath}/*{torget_name}*{j}*zf.fit* > {red_filepath}/{torget_name}_{j}zf')
        obj_zf=[f.strip() for f in open(f'{red_filepath}/{torget_name}_{j}zf').readlines() if f.strip()]
        
        if len(obj_f)!=0:
            RA=fits.getval(obj_zf[0], 'RA')
            DEC=fits.getval(obj_zf[0], 'DEC')
            for i in range(len(obj_zf)):
                wcs_solve(obj_zf[i],RA,DEC,red_filepath)
                temp_name=obj_zf[i][:-4]+'w'
                sex_phot(temp_name)

                starname=run_xmarch(temp_name+'.cat','II/349/ps1','X_WORLD',"Y_WORLD")
                if 'r'.lower() in j.lower():
                    lsfl=get_Zmag_limitmag(obj_zf[i],starname,'MAG_BEST','MAGERR_BEST','rmag')
                    list_zmag_fwhm_lim.append(lsfl)
                elif 'i'.lower() in j.lower():
                    lsfl=get_Zmag_limitmag(obj_zf[i],starname,'MAG_BEST','MAGERR_BEST','imag')
                    list_zmag_fwhm_lim.append(lsfl)
                else:
                    lsfl=get_Zmag_limitmag(obj_zf[i],starname,'MAG_BEST','MAGERR_BEST','gmag')
                    list_zmag_fwhm_lim.append(lsfl)

                show_fits(f'{temp_name}',torget_name,red_filepath,f'{os.path.basename(temp_name)}_comparePS1',lsfl[2])

    print("WCS and SEX done")



    
    for j in filter:
        run_cmd_Popen_fileno(f'ls {red_filepath}/*{torget_name}*{j}*zfw.fit* > {red_filepath}/{torget_name}_{j}zfw')
        obj_zfw=[f.strip() for f in open(f'{red_filepath}/{torget_name}_{j}zfw').readlines() if f.strip()]
        if len(obj_f)>1 :
            combine_fits_withoutwcs(obj_zfw,red_filepath,f'{file_name}_{torget_name}_{j}_com.fit')
            temp_name=f'{red_filepath}/{file_name}_{torget_name}_{j}_com'
            sex_phot(temp_name)

            starname=run_xmarch(temp_name+'.cat','II/349/ps1','X_WORLD',"Y_WORLD")

            if 'r'.lower() in j.lower():
                fzl_temp=get_Zmag_limitmag(temp_name+'.fit',starname,'MAG_BEST','MAGERR_BEST','rmag')
                list_zmag_fwhm_lim.append(fzl_temp)
            elif 'i'.lower() in j.lower():
                fzl_temp=get_Zmag_limitmag(temp_name+'.fit',starname,'MAG_BEST','MAGERR_BEST','imag')
                list_zmag_fwhm_lim.append(fzl_temp)
            else:
                fzl_temp=get_Zmag_limitmag(temp_name+'.fit',starname,'MAG_BEST','MAGERR_BEST','gmag')
                list_zmag_fwhm_lim.append(fzl_temp)

            show_fits(f'{red_filepath}/{file_name}_{torget_name}_{j}_com',\
                torget_name,red_filepath,f'{file_name}_{torget_name}_{j}_comparePS1',fzl_temp[2])
    print("COMBINE done")

    list_zmag_fwhm_lim_a=np.array(list_zmag_fwhm_lim)
    col_name=['NAME','EXPTIME',"ZMAG","LIMIT_MAG_MEAN","LIMIT_MAG_MEDIAN"]
    lzfl=pd.DataFrame(data=list_zmag_fwhm_lim_a,columns=col_name)
    lzfl.to_csv(f'{red_filepath}/{torget_name}-JZFL.csv')




    print("ALL done")