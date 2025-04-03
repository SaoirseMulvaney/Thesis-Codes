# %% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import glob
import os
from iris_fitting.fit_iris_lines import fit_raster
from iris_fitting import iris_get_mg_features_lv2 as get_mg
import asdf
from astropy.io import fits
import datetime as dt

if os.uname().sysname == 'Darwin':
    IRIS_data_loc = '/Users/dml/Data/IRIS/'
    output_loc = '/Users/dml/python_output/IRIS_output/'
else:
    IRIS_data_loc = '/mnt/scratch/data/mulvans8/Data/'
    output_loc = '/home/ug/mulvans8/python_output/'


# %% [markdown]
# Fit the different IRIS spectral lines

# %%
# Fit the Si IV spectral lines
def fit_iris(file, iris_window, event):
    a = fit_raster(file, iris_window)
    results_array, int_map, dopp_map, width_map, vnt_map, asym_map = a.fit_iris_data(fulldisk=True)

    plot_time = int_map.date.strftime('%Y%m%d_%H%M%S')

# Save the outputs as an asdf file for each run
    tree = {'results_array':results_array, 'int_map':int_map, 'dopp_map':dopp_map, 'width_map':width_map, 
            'vnt_map':vnt_map, 'asym_map':asym_map}
    with asdf.AsdfFile(tree) as asdf_file:  
        asdf_file.write_to(output_loc+event+'/IRIS_fitting_'+iris_window.replace(' ', '_')+'_'+plot_time+'.asdf',
                           all_array_compression='zlib')

    return

# %%
# Get the Mg II line properties
def fit_iris_mgii(file, iris_window, event, plot_time, onlyk=False, onlyh=False):

    lc, rp, bp = get_mg.iris_get_mg_features_lev2(file,onlyh=onlyh,onlyk=onlyk,fulldisk=True)

    lc_v = lc[0,:,:,0]
    lc_i = lc[0,:,:,1]

    rp_v = rp[0,:,:,0]
    rp_i = rp[0,:,:,1]

    bp_v = bp[0,:,:,0]
    bp_i = bp[0,:,:,1]

    a = fit_raster(file, iris_window)
    sp = fits.open(file, memmap=False)
    header = sp[0].header
 
    if onlyk:
        dv_k3_map = a.mk_iris_map(lc_v, header, header, fulldisk=True) # Upper chromosphere velocity
        k2_diff = rp_v - bp_v
        k2_sep_map = a.mk_iris_map(k2_diff, header, header, fulldisk=True) # Mid chromosphere velocity
        tree = {'dv_k3_map':dv_k3_map, 'k2_sep_map':k2_sep_map}
        with asdf.AsdfFile(tree) as asdf_file:  
            asdf_file.write_to(output_loc+event+'/IRIS_fitting_MgIIk_'+plot_time+'.asdf', all_array_compression='zlib')

    if onlyh:
        dv_h3_map = a.mk_iris_map(lc_v, header, header, fulldisk=True) # Upper chromosphere velocity
        h2_diff = rp_v - bp_v
        h2_sep_map = a.mk_iris_map(h2_diff, header, header, fulldisk=True) # Mid chromosphere velocity
        tree = {'dv_h3_map':dv_h3_map, 'h2_sep_map':h2_sep_map}
        with asdf.AsdfFile(tree) as asdf_file:  
            asdf_file.write_to(output_loc+event+'/IRIS_fitting_MgIIh_'+plot_time+'.asdf', all_array_compression='zlib')

    return

# %%

def fit_iris_data_fulldisk():

    iris_evts = ['20140317']

    for event in iris_evts:

        os.makedirs(output_loc+event+'/', exist_ok='True')

# Si IV 1394
        window = 'Si IV 1394'
        f_iris_raster = glob.glob(IRIS_data_loc+event+"/"+window.replace(' ', '_')+"/*.fits")
        f_iris_raster.sort()
        file = f_iris_raster[0]

        sp = fits.open(file, memmap=False)
        header = sp[0].header
        plot_time = dt.datetime.strftime(dt.datetime.strptime(header['DATE_OBS'],'%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

        files = os.path.exists(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_'+plot_time+'.asdf')
        if files == False:
            fit_iris(file, window, event)

# Si IV 1403
        window = 'Si IV 1403'
        f_iris_raster = glob.glob(IRIS_data_loc+event+"/"+window.replace(' ', '_')+"/*.fits")
        f_iris_raster.sort()
        file = f_iris_raster[0]

        sp = fits.open(file, memmap=False)
        header = sp[0].header
        plot_time = dt.datetime.strftime(dt.datetime.strptime(header['DATE_OBS'],'%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

        files = os.path.exists(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_'+plot_time+'.asdf')
        if files == False:
            fit_iris(file, window, event)

# Mg II h
        window = 'Mg II h'
        f_iris_raster = glob.glob(IRIS_data_loc+event+"/"+window.replace(' ', '_')+"/*.fits")
        f_iris_raster.sort()
        file = f_iris_raster[0]

        sp = fits.open(file, memmap=False)
        header = sp[0].header
        plot_time = dt.datetime.strftime(dt.datetime.strptime(header['DATE_OBS'],'%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

        files = os.path.exists(output_loc+event+'/IRIS_fitting_MgIIh_'+plot_time+'.asdf')
        if files == False:
            fit_iris_mgii(file, 'Mg II h', event, plot_time, onlyh=True)

# Mg II k
        window = 'Mg II k'
        f_iris_raster = glob.glob(IRIS_data_loc+event+"/"+window.replace(' ', '_')+"/*.fits")
        f_iris_raster.sort()
        file = f_iris_raster[0]

        sp = fits.open(file, memmap=False)
        header = sp[0].header
        plot_time = dt.datetime.strftime(dt.datetime.strptime(header['DATE_OBS'],'%Y-%m-%dT%H:%M:%S.%f'), '%Y%m%d_%H%M%S')

        files = os.path.exists(output_loc+event+'/IRIS_fitting_MgIIk_'+plot_time+'.asdf')
        if files == False:
            fit_iris_mgii(file, 'Mg II k', event, plot_time, onlyk=True)

if __name__ == "__main__":
    fit_iris_data_fulldisk()

