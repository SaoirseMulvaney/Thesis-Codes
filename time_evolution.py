
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


iris_evts = ['20131027', '20141020', '20151027', '20161024', '20171021', '20181025', '20191020', '20201019', '20211016', '20221009', '20231008', '20241016']
#iris_evts = ['20230329_111458']
#iris_evts = ['20231024_165953']

# Define the northern & southern hemispheres
upr_bl = [-1000, 0]
upr_tr = [1000,500]

lwr_bl = [-1000, -500]
lwr_tr = [1000, 0]

def get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map,event,window):

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

    upr_data = pd.DataFrame({'int':int_upr, 'vdopp':dopp_upr, 'vnt':vnt_upr, 'type':['Northern Hemisphere']*len(int_upr), 'date':[event]*len(int_upr), 'window':[window]*len(int_upr)})
    lwr_data = pd.DataFrame({'int':int_lwr, 'vdopp':dopp_lwr, 'vnt':vnt_lwr, 'type':['Southern Hemisphere']*len(int_lwr), 'date':[event]*len(int_lwr), 'window':[window]*len(int_lwr)})
    df = pd.concat([upr_data, lwr_data], ignore_index=True)

    return df

def get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map,sep_map,event,window):

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

    upr_data = pd.DataFrame({'dv':dv_upr, 'sep':sep_upr, 'type':['Northern Hemisphere']*len(dv_upr), 'date':[event]*len(dv_upr), 'window':[window]*len(dv_upr)})
    lwr_data = pd.DataFrame({'dv':dv_lwr, 'sep':sep_lwr, 'type':['Southern Hemisphere']*len(dv_lwr), 'date':[event]*len(dv_lwr), 'window':[window]*len(dv_lwr)})

    df = pd.concat([upr_data, lwr_data], ignore_index=True)

    return df

def plot_kdes_by_hemisphere(df, param, label, wvl):

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    sns.kdeplot(df, x='int', hue='date', fill = True, alpha = 0.2, bw_adjust=0.5, ax=axs[0])
    axs[0].set_xlabel('Intensity')
    axs[0].set_ylabel('Density')
    axs[0].set_xlim(0,300)

    sns.kdeplot(df, x='vdopp', hue='date', fill = True, alpha = 0.2, bw_adjust=0.5, ax=axs[1])
    axs[1].set_xlabel('Line-of-sight velocity (km/s)')
    axs[2].set_ylabel('Density')
    axs[1].set_xlim(-20,20)

    sns.kdeplot(df, x='vnt', hue='date', fill = True, alpha = 0.2, bw_adjust=0.5, ax=axs[2])
    axs[2].set_xlabel('Line width ($\AA$)')
    axs[2].set_ylabel('Density')
    axs[2].set_xlim(0,0.1)
    
    fig.suptitle(f' Si IV {wvl} KDE Evolution over The Solar Cyle')
    plt.tight_layout()
    plt.savefig(os.path.join(output_loc+ f"summary_kdes/SiIV{wvl}_all_years.png"))
    plt.close()

def plot_kdes_mgii_over_time(df, wvl):

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))

    sns.kdeplot(df, x='dv', hue='date', fill = True, alpha = 0.2, bw_adjust=0.5, ax=axs[0])
    axs[0].set_xlabel('Asymmetry')
    axs[0].set_ylabel('Density')
    axs[0].set_xlim(-10,10)

    sns.kdeplot(df, x='sep', hue='date', fill = True, alpha = 0.2, bw_adjust=0.5, ax=axs[1])
    axs[1].set_xlabel('K-H Separation (km/s)')
    axs[2].set_ylabel('Density')
    axs[1].set_xlim(0,60)
    
    fig.suptitle(f' Si IV {wvl} KDE Evolution over The Solar Cyle')
    plt.tight_layout()
    plt.savefig(output_loc+ f"time_evolution_{wvl}_{param}.png")

    fig.suptitle(f' Mg II {wvl} KDE Evolution over The Solar Cyle')
    plt.tight_layout()
    plt.savefig(os.path.join(output_loc+ f"summary_kdes/MgII{wvl}_all_years.png"))
    plt.close()

all_si1393 = []
all_si1403 = []
all_mgk = []
all_mgh = []

for event in iris_evts:

    os.makedirs(output_loc+event+'/', exist_ok='True')

# Si IV 1393
    window = 'Si IV 1393'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        int_map = asdf_file['int_map']
        dopp_map = asdf_file['dopp_map']
        vnt_map = asdf_file['vnt_map']

    df = get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map,event,window)
    all_si1393.append(df)

# Si IV 1403
    window = 'Si IV 1403'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_'+window.replace(' ', '_')+'_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        int_map = asdf_file['int_map']
        dopp_map = asdf_file['dopp_map']
        vnt_map = asdf_file['vnt_map']

    df = get_dataframe_siiv(upr_bl,upr_tr,lwr_bl,lwr_tr,int_map,dopp_map,vnt_map,event,window)
    all_si1403.append(df)

# Mg II k
    window = 'Mg II k'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIk_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        dv_map = asdf_file['dv_k3_map']
        sep_map = asdf_file['k2_sep_map']

    dv_map_new = Map(dv_map.data.T, int_map.meta)
    sep_map_new = Map(sep_map.data.T, int_map.meta)

    df = get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map_new,sep_map_new,event,window)
    all_mgk.append(df)

# Mg II h
    window = 'Mg II h'

    asdf_path = glob.glob(output_loc+event+'/IRIS_fitting_MgIIh_*.asdf')

    with asdf.open(asdf_path[0]) as asdf_file:
        dv_map = asdf_file['dv_h3_map']
        sep_map = asdf_file['h2_sep_map']

    dv_map_new = Map(dv_map.data.T, int_map.meta)
    sep_map_new = Map(sep_map.data.T, int_map.meta)

    df_mgh = get_dataframe_mgii(upr_bl,upr_tr,lwr_bl,lwr_tr,dv_map_new,sep_map_new,event,window)
    all_mgh.append(df)

os.makedirs(os.path.join(output_loc, "summary_kdes"), exist_ok='True')

if all_si1393:
    df_all_si1393 = pd.concat(all_si1393, ignore_index=True)
    plot_kdes_siiv_over_time(df_all_si1393, 'SiIV1393')

if all_si1403:
    df_all_si1403 = pd.concat(all_si1403, ignore_index=True)
    plot_kdes_siiv_over_time(df_all_si1403, 'SiIV1403')

if all_mgk:
    df_all_mgk = pd.concat(all_mgk, ignore_index=True)
    plot_kdes_mgii_over_time(df_all_mgk, 'MgIIk')

if all_mgh:
    df_all_mgh = pd.concat(all_mgh, ignore_index=True)
    plot_kdes_mgii_over_time(df_all_mgh, 'MgIIh')

