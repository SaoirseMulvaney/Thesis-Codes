# %% [markdown]
# **Notebook to fit IRIS data using multiple cores**
Original Author: David Long
Contact: david.long@dcu.ie

Updates Author: Saoirse Mulvaney
Contact: saoirse.mulvaney8@mail.dcu.ie

# %%
import os
import glob
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.colors as colors
import numpy as np
import datetime as dt
import asdf
from sunpy.map import Map
import astropy. units as u
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib as mpl
from astropy.visualization import ImageNormalize, AsinhStretch

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

# Si IV 1393
    window = 'Si IV 1393'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
             doppler_map = asdf_file['dopp_map']
             int_map = asdf_file['int_map']
             vnt_map = asdf_file['vnt_map']

#             incorrect_map= doppler_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             doppler_map = Map(r_data, r_map.meta)

#             incorrect_map= int_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             int_map = Map(r_data, r_map.meta)

#             incorrect_map= vnt_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             vnt_map = Map(r_data, r_map.meta)

             aspect_ratio_dopp = doppler_map.meta['CDELT2'] / np.abs(doppler_map.meta['CDELT1'])
             aspect_ratio_int = int_map.meta['CDELT2'] / np.abs(int_map.meta['CDELT1'])
             aspect_ratio_vnt = vnt_map.meta['CDELT2'] / np.abs(vnt_map.meta['CDELT1'])
             

             fig = plt.figure(figsize=(15, 5))

             ax1 = fig.add_subplot(131, projection=doppler_map)
             doppler_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             doppler_map.plot_settings['norm'] = Normalize(-20, 20)
             doppler_map.plot(axes=ax1, aspect=aspect_ratio_dopp)
             plt.colorbar(cmap='coolwarm', ax=ax1, orientation='vertical')

             ax1.set_xlabel('Solar X (arsec)')
             ax1.set_ylabel('Solar Y (arsec)')

             ax2 = fig.add_subplot(132, projection=int_map)
             int_map.plot_settings['cmap'] = mpl.colormaps['Reds']
             norm = ImageNormalize(vmin = 0, vmax = 300, stretch=AsinhStretch())
             int_map.plot_settings['norm'] = norm
             int_map.plot(axes=ax2, aspect=aspect_ratio_int)
             plt.colorbar(cmap='Reds', ax=ax2, orientation='vertical')

             ax2.set_xlabel('Solar X (arsecs)')
             ax2.set_ylabel('Solar Y (arsecs)')

             ax3 = fig.add_subplot(133, projection=vnt_map)
             vnt_map.plot_settings['cmap'] = mpl.colormaps['viridis']
             vnt_map.plot_settings['norm'] = Normalize(0, 0.1)
             vnt_map.plot(axes=ax3, aspect=aspect_ratio_vnt)
             plt.colorbar(cmap='viridis', ax=ax3, orientation='vertical')

             ax3.set_xlabel('Solar X (arsec)')
             ax3.set_ylabel('Solar Y (arsec)')

             ax1.set_title("Doppler Map")
             ax2.set_title("Intensity Map")
             ax3.set_title("Line Width Map")
             fig.tight_layout()
             plt.savefig(output_loc+event+'/20241016.1394.f.png')


# Si IV 1403
    window = 'Si IV 1403'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
             doppler_map = asdf_file['dopp_map']
             int_map = asdf_file['int_map']
             vnt_map = asdf_file['vnt_map']

#             incorrect_map= doppler_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             doppler_map = Map(r_data, r_map.meta)
#
#             incorrect_map= int_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             int_map = Map(r_data, r_map.meta)
#
#             incorrect_map= vnt_map
#             r_map = incorrect_map.rotate(angle = 90*u.deg)
#             r_data = np.fliplr(r_map.data)
#             vnt_map = Map(r_data, r_map.meta)

             aspect_ratio_dopp = doppler_map.meta['CDELT2'] / np.abs(doppler_map.meta['CDELT1'])
             aspect_ratio_int = int_map.meta['CDELT2'] / np.abs(int_map.meta['CDELT1'])
             aspect_ratio_vnt = vnt_map.meta['CDELT2'] / np.abs(vnt_map.meta['CDELT1'])
             

             fig = plt.figure(figsize=(15, 5))

             ax1 = fig.add_subplot(131, projection=doppler_map)
             doppler_map.plot_settings['cmap'] = mpl.colormaps['coolwarm']
             doppler_map.plot_settings['norm'] = Normalize(-20, 20)
             doppler_map.plot(axes=ax1, aspect=aspect_ratio_dopp)
             plt.colorbar(cmap='coolwarm', ax=ax1, orientation='vertical')

             ax1.set_xlabel('Solar X (arsec)')
             ax1.set_ylabel('Solar Y (arsec)')

             ax2 = fig.add_subplot(132, projection=int_map)
             int_map.plot_settings['cmap'] = mpl.colormaps['Reds']
             norm = ImageNormalize(vmin = 0, vmax = 300, stretch=AsinhStretch())
             int_map.plot_settings['norm'] = norm
             int_map.plot(axes=ax2, aspect=aspect_ratio_int)
             plt.colorbar(cmap='Reds', ax=ax2, orientation='vertical')

             ax2.set_xlabel('Solar X (arsec)')
             ax2.set_ylabel('Solar Y (arsec)')

             ax3 = fig.add_subplot(133, projection=vnt_map)
             vnt_map.plot_settings['cmap'] = mpl.colormaps['viridis']
             vnt_map.plot_settings['norm'] = Normalize(0, 0.1)
             vnt_map.plot(axes=ax3, aspect=aspect_ratio_vnt)
             plt.colorbar(cmap='viridis', ax=ax3, orientation='vertical')

             ax3.set_xlabel('Solar X (arsec)')
             ax3.set_ylabel('Solar Y (arsec)')

             ax1.set_title("Doppler Map")
             ax2.set_title("Intensity Map")
             ax3.set_title("Line Width Map")
             fig.tight_layout()
             plt.savefig(output_loc+event+'/20241016.1403.f.png')

 