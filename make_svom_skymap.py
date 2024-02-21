#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 16:45:39 2022

@author: dt270490
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:37:49 2022

@author: dt270490
"""
import os
import numpy as np
import healpy as hp
from healpy.newvisufunc import projview, newprojplot
import astropy.units as u
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord

   
def get_grb_pos(grb_pos_file):
    
    grb_pos = pd.read_csv(grb_pos_file, delimiter=",")
    
    return grb_pos

def compute_theta_phi(grb_pos:pd.DataFrame, coord_sys:str):
    """
    

    Parameters
    ----------
    grb_pos;pd.DataFrame : TYPE
        DESCRIPTION.
    coord_sys : str
        DESCRIPTION.

    Returnsos.getcwd()
    -------
    None.

    """
    theta = []
    phi = []
    
    if coord_sys == 'gal':
        coord_theta = "b"
        coord_phi = "l"
    elif coord_sys == "equ":
        coord_theta = "dec"
        coord_phi = "ra"
        
    mask_long = grb_pos[coord_phi] >180 # need to separate
    if type(mask_long) == np.bool_ and mask_long:
        theta = np.deg2rad((grb_pos[coord_theta]* -1) + 90)
        phi = np.deg2rad(grb_pos[coord_phi]-360)
    elif type(mask_long) == np.bool_ and not mask_long:
        theta = np.deg2rad((grb_pos[coord_theta]* -1) + 90)
        phi = np.deg2rad(grb_pos[coord_phi])
    
    if not type(mask_long) == np.bool_:
        if mask_long.any():
            theta1 = np.deg2rad((grb_pos[mask_long][coord_theta]* -1) + 90)
            phi1 = np.deg2rad(grb_pos[mask_long][coord_phi]-360)
        if not mask_long.all():
            theta2 = np.deg2rad((grb_pos[~mask_long][coord_theta]* -1) + 90)
            phi2 = np.deg2rad(grb_pos[~mask_long][coord_phi])
        
        theta = np.concatenate((theta1.values,theta2.values))
        phi = np.concatenate((phi1.values,phi2.values))

    
    return theta, phi
    

class plots:
    """ A class to compute the SVOM skymap plots of the GRB and ToO catalog"""

    def __init__(self, coord_syst:str):
        """

        Parameters
        ----------
        coord_syst : str
            Whether you want a skymap in galactic ("gal") or equatorial ("equ")
            coordinates

        Returns
        -------
        None.

        """
        # plt.rcParams["font.family"] = "Yrsa"
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif Times New Roman')

        self.fontsize_label = 15
        self.fontsize_title = 30
        self.fontlabel = {
                          'color':  'black',
                          'weight': 'normal',
                          'size': self.fontsize_label,
                          }
        self.fontsize_legend = 30
        self.fontsize_tick = 15
        self.figsize = 14
        # the marker styles
        self.marker_lastgrb = "*"
        self.marker_grb = 'h'
        self.marker_too = 'v'
        self.marker_lasttoo = 'o'
        # The marker sizes
        self.marker_grb_size = 7
        self.marker_lastgrb_size = 15
        self.marker_too_size = 6
        self.marker_lasttoo_size = 7
        # The marker colors
        self.grb_color = "#f58021"
        self.lastgrb_color = "green"
        self.too_color = "#febe11"
        self.lasttoo_color = "green"
        
        self.label_color = "crimson"
        self.lw = 3
        self.majortick_width = 2
        self.minortick_width = 1
        self.majortick_size = 10
        self.minortick_size = 7
        
        self.coord_syst = coord_syst
        
    def plot_sky_position(self,grb_pos:pd.DataFrame,
                          too_pos:pd.DataFrame,
                          outdir:str,
                          savefig:bool):
        """
        Plot the GRB and ToO positions in an all sky map
    
        Parameters
        ----------
        grb_pos : pd.dataframe
            Gathers the burst_id and the equatorial coordinates of the SVOM GRBs
        too_pos : pd.dataframe
            Gathers the ToO_id and the equatorial coordinates of the SVOM ToOs.
        outdir: str
            The path where you want to save the figures.
        savefig: Boolean
            Whether you want to save the figure or not.
        Returns
        -------
        None.
    
        """        
        # plot settings
        sdf_map = hp.read_map('dustmap/EBV_SFD98_1_512.fits')
        if self.coord_syst =='gal':
            coord = 'G'
            x_label = "l [deg]"
            y_label = "b [deg]"
        elif self.coord_syst == 'equ':
            coord = 'C'
            x_label = r"$\alpha$ [deg]"
            y_label = r"$\delta$ [deg]"
        else:
            coord = 'C'
            x_label = r"$\alpha$ [deg]"
            y_label = r"$\delta$ [deg]"
        projview(sdf_map,
                 xsize=1000,
                 coord=["C",coord],
                 graticule=True,
                 graticule_labels=True,
                 projection_type="mollweide",
                 min=0,
                 max=2,
                 longitude_grid_spacing=30,
                 latitude_grid_spacing=15,
                 xlabel=x_label,
                 ylabel=y_label,
                 xtick_label_color=self.label_color,
                 ytick_label_color=self.label_color,
                 cbar=True,
                 cb_orientation='vertical',
                 unit="E(B-V)",
                 cmap='PuBu',
                 fontsize={'title':self.fontsize_label,
                           "xlabel":self.fontsize_label,
                           "ylabel":self.fontsize_label,
                           "xtick_label": self.fontsize_tick,
                           "ytick_label": self.fontsize_tick,
                           "cbar_label":self.fontsize_label,
                           "cbar_tick_label":self.fontsize_tick},
                 override_plot_properties={"figure_width": self.figsize,
                                           "figure_size_ratio": 0.5}
                 )
        
        #------------- The GRBs
        
        # Plot all the GRBs except the last one
        if len(grb_pos)>1:
            theta_grb,phi_grb = compute_theta_phi(grb_pos.iloc[:-1],
                                                  self.coord_syst)
            newprojplot(theta=theta_grb, phi=phi_grb,
                        marker=self.marker_grb,
                        color=self.grb_color,
                        markersize=self.marker_grb_size,
                        ls='',
                        mfc=self.grb_color);
        # Plot the last GRB
        theta_lastgrb,phi_lastgrb = compute_theta_phi(grb_pos.iloc[-1],
                                                      self.coord_syst)
        newprojplot(theta=theta_lastgrb, phi=phi_lastgrb,
                    marker=self.marker_lastgrb,
                    color=self.grb_color,
                    markersize=self.marker_lastgrb_size,
                    ls='',
                    mec = self.lastgrb_color);
        
        #------------- The ToOs
        # Plot all the ToO except the last one
        if len(too_pos)>1:
            theta_too,phi_too = compute_theta_phi(too_pos.iloc[:-1],
                                                  self.coord_syst)
            newprojplot(theta=theta_too, phi=phi_too,
                        marker=self.marker_too,
                        color=self.too_color,
                        markersize=self.marker_too_size,
                        ls='',
                        mfc="white");
        # Plot the last ToO
        theta_lasttoo,phi_lasttoo = compute_theta_phi(too_pos.iloc[-1],
                                                      self.coord_syst)
        newprojplot(theta=theta_lasttoo, phi=phi_lasttoo,
                    marker=self.marker_lasttoo,
                    color=self.too_color,
                    markersize=self.marker_lasttoo_size,
                    ls='',
                    mec = self.lasttoo_color);
        
        if savefig:
            
            plt.savefig(outdir+"SVOM_skymap_"+self.coord_syst+".png", dpi=200)
        
if __name__ == '__main__':
    
    # Get the SVOM GRB sky position
    grb_pos = get_grb_pos('data_example/svom_grb_position.csv')
    
    # add galactic coordinates
    coords = SkyCoord(ra=grb_pos.ra,
                      dec=grb_pos.dec,
                      unit=(u.deg, u.deg))
    grb_pos['l'] = coords.galactic.l.degree
    grb_pos['b'] = coords.galactic.b.degree
    # Get the SVOM ToO sky position
    too_pos = get_grb_pos('data_example/svom_too_position.csv')
    # add galactic coordinates
    coords = SkyCoord(ra=too_pos.ra,
                      dec=too_pos.dec,
                      unit=(u.deg, u.deg))
    too_pos['l'] = coords.galactic.l.degree
    too_pos['b'] = coords.galactic.b.degree
    
    #make the plot
    svom_skymap = plots("equ")
    svom_skymap.plot_sky_position(grb_pos,too_pos,os.getcwd()+"/",True)