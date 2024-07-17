
def run_cmd_Popen_fileno(cmd_string):
    import subprocess
    return subprocess.Popen(cmd_string, shell=True, stdout=None, stderr=None).wait()

def sex_phot(sf):
        run_cmd_Popen_fileno(f'sex {sf}.fit -c tfind.sex -PARAMETERS_NAME findstar.param -CATALOG_NAME {sf}.cat')


def wcs_solve(sf,RA,DEC,red_filepath):
    run_cmd_Popen_fileno(f'solve-field {sf[:-4]}.fit -N {sf[:-4]}w.fit -O -p -3 {RA} -4 {DEC} -5 1')

    run_cmd_Popen_fileno('rm %s/*.match %s/*-indx.xyls %s/*.corr %s/*.axy %s/*.wcs %s/*.solved  %s/*.rdls'\
                     %(red_filepath,red_filepath,red_filepath,\
                      red_filepath,red_filepath,red_filepath,\
                      red_filepath))
