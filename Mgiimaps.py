
Original Author: David Long
Contact: david.long@dcu.ie

Updates Author: Saoirse Mulvaney
Contact: saoirse.mulvaney8@mail.dcu.ie

# %% [markdown]
# **Notebook to fit IRIS data using multiple cores**

# %%
import os
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import numpy as np
import datetime as dt
import asdf
from sunpy.map import Map
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib as mpl


if os.uname().sysname == 'Darwin':
    IRIS_data_loc = '/Users/dml/Data/IRIS/'
    output_loc = '/Users/dml/python_output/IRIS_output/'
else:
    IRIS_data_loc = '/home/ug/mulvans8/python_output/'
    output_loc = '/home/ug/mulvans8/python_output/'


iris_evts = ['20241016']
#iris_evts = ['20230329_111458']
#iris_evts = ['20231024_165953']

for event in iris_evts:
    os.makedirs(output_loc+event+'/', exist_ok='True')

# Mg II K
    window = 'Mg II k'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIk_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
             dv_k3_map = asdf_file['dv_k3_map']
             k2_sep_map = asdf_file['k2_sep_map']

             aspect_ratio_dv_k3 = dv_k3_map.meta['CDELT1'] / np.abs(dv_k3_map.meta['CDELT2'])
             aspect_ratio_k2_sep = k2_sep_map.meta['CDELT1'] / np.abs(k2_sep_map.meta['CDELT2'])
             

             fig = plt.figure(figsize=(15, 5))

             ax1 = fig.add_subplot(121, projection=dv_k3_map)
             dv_k3_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             dv_k3_map.plot_settings['norm'] = Normalize(-10, 10)
             fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-10, 10), cmap='coolwarm'),
                           ax=ax1, orientation='vertical')
             dv_k3_map.plot(axes=ax1, aspect=aspect_ratio_dv_k3)

             ax1.set_xlabel('Solar X (arsec)')
             ax1.set_ylabel('Solar Y (arsec)')

             ax2 = fig.add_subplot(122, projection=k2_sep_map)
             k2_sep_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             k2_sep_map.plot_settings['norm'] = Normalize(0, 60)
             fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0, 60), cmap='coolwarm'),
                           ax=ax2, orientation='vertical')
             k2_sep_map.plot(axes=ax2, aspect=aspect_ratio_k2_sep)

             ax2.set_xlabel('Solar X (arsec)')
             ax2.set_ylabel('Solar Y (arsec)')

             ax1.set_title("K Asymmetry Map")
             ax2.set_title("K-H Separation Map")
             fig.tight_layout()
             plt.savefig(output_loc+event+'/20241016.mgk.f.png')
             
# Mg II H
    window = 'Mg II h'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIh_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
             dv_h3_map = asdf_file['dv_h3_map']
             h2_sep_map = asdf_file['h2_sep_map']

             aspect_ratio_dv_h3 = dv_h3_map.meta['CDELT1'] / np.abs(dv_h3_map.meta['CDELT2'])
             aspect_ratio_h2_sep = h2_sep_map.meta['CDELT1'] / np.abs(h2_sep_map.meta['CDELT2'])
             

             fig = plt.figure(figsize=(15, 5))

             ax1 = fig.add_subplot(121, projection=dv_h3_map)
             dv_h3_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             dv_h3_map.plot_settings['norm'] = Normalize(-10, 10)
             fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(-10, 10), cmap='coolwarm'),
                           ax=ax1, orientation='vertical')
             dv_h3_map.plot(axes=ax1, aspect=aspect_ratio_dv_h3)

             ax1.set_xlabel('Solar X (arsec)')
             ax1.set_ylabel('Solar Y (arsec)')

             ax2 = fig.add_subplot(122, projection=h2_sep_map)
             h2_sep_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             h2_sep_map.plot_settings['norm'] = Normalize(0, 60)
             fig.colorbar(mpl.cm.ScalarMappable(norm=Normalize(0, 60), cmap='coolwarm'),
                           ax=ax2, orientation='vertical')
             h2_sep_map.plot(axes=ax2, aspect=aspect_ratio_h2_sep)

             ax2.set_xlabel('Solar X (arsec)')
             ax2.set_ylabel('Solar Y (arsec)')

             ax1.set_title("H Asymmetry Map")
             ax2.set_title("H-K Separation Map")
             fig.tight_layout()
             plt.savefig(output_loc+event+'/20241016.mgh.f.png')
