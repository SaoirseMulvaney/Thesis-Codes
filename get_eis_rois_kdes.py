
Original Author: David Long
Contact: david.long@dcu.ie

Updates Author: Saoirse Mulvaney
Contact: saoirse.mulvaney8@mail.dcu.ie

# %%
import os
import glob
import asdf
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sunpy.map import Map

if os.uname().sysname == 'Darwin':
    IRIS_data_loc = '/Users/dml/Data/IRIS/'
    output_loc = '/Users/dml/python_output/IRIS_output/'
else:
    IRIS_data_loc = '/home/ug/mulvans8/python_output/'
    output_loc = '/home/ug/mulvans8/python_output/'


iris_evts = ['20241016']
#iris_evts = ['20230329_111458']
#iris_evts = ['20231024_165953']

# Define the northern & southern hemispheres
upr_bl = [-1000, 0]
upr_tr = [1000,500]

lwr_bl = [-1000, -500]
lwr_tr = [1000, 0]

def get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map):

    upr_bottom_left = SkyCoord(upr_bl[0] * u.arcsec, upr_bl[1] * u.arcsec, frame=int_map.coordinate_frame)
    upr_top_right = SkyCoord(upr_tr[0] * u.arcsec, upr_tr[1] * u.arcsec, frame=int_map.coordinate_frame)

    int_upr_map = int_map.submap(upr_bottom_left, top_right = upr_top_right)
    dopp_upr_map = dopp_map.submap(upr_bottom_left, top_right = upr_top_right)
    vnt_upr_map = vnt_map.submap(upr_bottom_left, top_right = upr_top_right)

    int_upr = int_upr_map.data.flatten()
    dopp_upr = dopp_upr_map.data.flatten()
    vnt_upr = vnt_upr_map.data.flatten()

    lwr_bottom_left = SkyCoord(lwr_bl[0] * u.arcsec, lwr_bl[1] * u.arcsec, frame=int_map.coordinate_frame)
    lwr_top_right = SkyCoord(lwr_tr[0] * u.arcsec, lwr_tr[1] * u.arcsec, frame=int_map.coordinate_frame)

    int_lwr_map = int_map.submap(lwr_bottom_left, top_right = lwr_top_right)
    dopp_lwr_map = dopp_map.submap(lwr_bottom_left, top_right = lwr_top_right)
    vnt_lwr_map = vnt_map.submap(lwr_bottom_left, top_right = lwr_top_right)

    int_lwr = int_lwr_map.data.flatten()
    dopp_lwr = dopp_lwr_map.data.flatten()
    vnt_lwr = vnt_lwr_map.data.flatten()

    upr_data = pd.DataFrame({'int':int_upr, 'vdopp':dopp_upr, 'vnt':vnt_upr, 'type':np.tile('Northern Hemisphere',len(int_upr)), 'date':np.tile('20241016', len(int_upr)), 'window':np.tile(window, len(int_upr))})
    lwr_data = pd.DataFrame({'int':int_lwr, 'vdopp':dopp_lwr, 'vnt':vnt_lwr, 'type':np.tile('Southern Hemisphere',len(int_lwr)), 'date':np.tile('20241016', len(int_lwr)), 'window':np.tile(window, len(int_lwr))})

    df = pd.concat([upr_data, lwr_data], ignore_index=True)

    return df

def get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map,sep_map):

    upr_bottom_left = SkyCoord(upr_bl[0] * u.arcsec, upr_bl[1] * u.arcsec, frame=dv_map.coordinate_frame)
    upr_top_right = SkyCoord(upr_tr[0] * u.arcsec, upr_tr[1] * u.arcsec, frame=dv_map.coordinate_frame)

    dv_upr_map = dv_map.submap(upr_bottom_left, top_right = upr_top_right)
    sep_upr_map = sep_map.submap(upr_bottom_left, top_right = upr_top_right)

    dv_upr = dv_upr_map.data.flatten()
    sep_upr = sep_upr_map.data.flatten()

    lwr_bottom_left = SkyCoord(lwr_bl[0] * u.arcsec, lwr_bl[1] * u.arcsec, frame=dv_map.coordinate_frame)
    lwr_top_right = SkyCoord(lwr_tr[0] * u.arcsec, lwr_tr[1] * u.arcsec, frame=dv_map.coordinate_frame)

    dv_lwr_map = dv_map.submap(lwr_bottom_left, top_right = lwr_top_right)
    sep_lwr_map = sep_map.submap(lwr_bottom_left, top_right = lwr_top_right)

    dv_lwr = dv_lwr_map.data.flatten()
    sep_lwr = sep_lwr_map.data.flatten()

    upr_data = pd.DataFrame({'dv':dv_upr, 'sep':sep_upr, 'type':np.tile('Northern Hemisphere',len(dv_upr)), 'date':np.tile('20241016', len(dv_upr)), 'window':np.tile(window, len(dv_upr))})
    lwr_data = pd.DataFrame({'dv':dv_lwr, 'sep':sep_lwr, 'type':np.tile('Southern Hemisphere',len(dv_lwr)), 'date':np.tile('20241016', len(dv_lwr)), 'window':np.tile(window, len(dv_lwr))})

    df = pd.concat([upr_data, lwr_data], ignore_index=True)

    return df

def plot_kdes_siiv(df, wvl):

    fig = plt.figure(figsize=(15, 5))

    ax01 = fig.add_subplot(131, label='a)', title = ' ')
    h1=sns.kdeplot(df, x='int', hue='type', ax=ax01, fill = True, alpha = 0.4, legend='brief', bw_adjust=.5)
    h1.set(xlabel='Intensity')
    ax01.set_xlim(0,300)
    sns.move_legend(ax01, "lower center", bbox_to_anchor=(2, 1), ncol=2, title='16th October 2024', frameon=False)

    ax02 = fig.add_subplot(132, label='b)', title = ' ')
    h2=sns.kdeplot(df, x='vdopp', hue='type', ax=ax02, fill = True, alpha = 0.4, legend=False, bw_adjust=.5)
    h2.set(xlabel='Line-of-sight velocity (km/s)')
    h2.set(ylabel = 'Density')
    ax02.set_xlim(-20,20)

    ax03 = fig.add_subplot(133, label='c)', title = ' ')
    h3=sns.kdeplot(df, x='vnt', hue='type', ax=ax03, fill = True, alpha = 0.4, legend=False, bw_adjust=.5)
    h3.set(xlabel=r'Line width ($\AA$)')
    h3.set(ylabel = 'Density')
    ax03.set_xlim(0,0.1)
    fig.subplots_adjust(wspace=0.5)
    plt.savefig(output_loc+event+'/20241016.'+wvl+'.kde.png')

def plot_kdes_mgii(df, wvl):

    fig = plt.figure(figsize=(15, 5))

    ax01 = fig.add_subplot(121, label='a)', title = ' ')
    h1=sns.kdeplot(df, x='dv', hue='type', ax=ax01, fill = True, alpha = 0.4, legend='brief', bw_adjust=.5)
    h1.set(xlabel='Asymmetry')
    ax01.set_xlim(-10,10)
    sns.move_legend(ax01, "lower center", bbox_to_anchor=(1.3, 1), ncol=2, title='16th October 2024', frameon=False)

    ax02 = fig.add_subplot(122, label='c)', title = ' ')
    h2=sns.kdeplot(df, x='sep', hue='type', ax=ax02, fill = True, alpha = 0.4, legend=False, bw_adjust=.5)
    h2.set(xlabel=r'K-H Separation (km/s)')
    h2.set(ylabel = 'Density')
    ax02.set_xlim(0,60)
    fig.subplots_adjust(wspace=0.5)

    plt.savefig(output_loc+event+'/20241016.'+wvl+'.kde.png')



for event in iris_evts:

    os.makedirs(output_loc+event+'/', exist_ok='True')

# Si IV 1393
    window = 'Si IV 1393'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        int_map = asdf_file['int_map']
        dopp_map = asdf_file['dopp_map']
        vnt_map = asdf_file['vnt_map']

    df_si1393 = get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map)

    plot_kdes_siiv(df_si1393, '1393')

# Si IV 1403
    window = 'Si IV 1403'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        int_map = asdf_file['int_map']
        dopp_map = asdf_file['dopp_map']
        vnt_map = asdf_file['vnt_map']

    df_si1403 = get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map)

    plot_kdes_siiv(df_si1403, '1403')


# Mg II k
    window = 'Mg II k'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIk_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        dv_map = asdf_file['dv_k3_map']
        sep_map = asdf_file['k2_sep_map']

    dv_map_new = Map(dv_map.data.T, int_map.meta)
    sep_map_new = Map(sep_map.data.T, int_map.meta)

    df_mgk = get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map_new,sep_map_new)

    plot_kdes_mgii(df_mgk, 'k')

# Mg II h
    window = 'Mg II h'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIh_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        dv_map = asdf_file['dv_h3_map']
        sep_map = asdf_file['h2_sep_map']

    dv_map_new = Map(dv_map.data.T, int_map.meta)
    sep_map_new = Map(sep_map.data.T, int_map.meta)

    df_mgh = get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map_new,sep_map_new)

    plot_kdes_mgii(df_mgh, 'h')


