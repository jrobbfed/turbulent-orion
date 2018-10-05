import numpy as np
import astropy.units as u
import pickle
from spectral_cube import SpectralCube
from glob import glob
import os
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
import subregion
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import astropy
from turbustat.statistics import VCS
from turbustat.statistics import PowerSpectrum 

distance = 414*u.pc

def main():

    high_cut = 0.04/u.pix
    high_cut_12co = 0.04/u.pix
    high_cut_13co = 0.04/u.pix
    high_cut_c18o = 0.025/u.pix

    refit = True
    pspec_str = "beamcorrect_apodizetukey"
    if refit:
        refit_str = "_beamcorrect_apodizetukey_pspechighcut{}".format(high_cut.value)
    else:
        refit_str = "_beamcorrect_apodizetukey"

    
### Plotting scatter plots of SCF and PSPEC slopes versus 
### surface density of YSOs (from Spitzer Orion) in Orion A regions.
# stat_list = glob("scf_results/scf_13co*davis*.p") +\
#         glob("scf_results/scf_13co*unger*.p")
# plot_stat_feedback(stat_list, plot_file="scf_13co_davis_nstars.pdf")
##Figure in paper.


    mols = ["12co", "13co", "c18o"]

    mols_dict = {
        "12co":r"$^{12}$CO",
        "13co":r"$^{13}$CO",
        "c18o":r"C$^{18}$O"}

    stat_names = ["scf", "pspec"]


    ytick_arr = np.array(
            [[[-0.10,-0.15],
              [-0.10,-0.20],
              [-0.20,-0.30]],
             [[-3.0,-4.0],
              [-3.0,-3.5],
              [-3.2,-3.6]]])
    
    # vcs = "vcs_results/vcs_12co_jrf_south.p"
    # plot_vcs_spectrum(vcs, ax=None, return_ax=False, plot_file="vcs_spectrum.pdf",
    #     mpl_style='presentation',
    #     plot_points=True, plot_fit=True,
    #     scatter_kwargs={'color':'red'}, plot_kwargs={'color':'black'},
    #     show_fitrange=True, refit=False, axvline_kwargs={'color':'black'})

             
    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):
            
    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))
            
    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'
 
    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=True, feedback_label=r"n$_{\rm YSO}$ (deg$^{-2}$)",
    #                 logx=True, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30}, within_data=True,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))
            
    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))

    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_nstars_logdensity{}.pdf".format(refit_str))
    # plt.close()

### Plotting scatter plots of SCF and PSPEC slopes versus 
### number of YSOs (from Spitzer Orion) in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))

    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=False, feedback_label=r"N$_{\rm YSO}$",
    #                 logx=True, within_data=True,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))
    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_nstars{}.pdf".format(refit_str))
    # plt.close()    

## Plotting scatter plots of SCF and PSPEC slopes versus 
## region area (from Spitzer Orion) in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))

    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=False, feedback_label=r"surface area (deg$^2$)",
    #                 logx=False, feedback_mode="area", within_data=False,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))
    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_surface_area{}.pdf".format(refit_str))
    # plt.close()  

## Plotting scatter plots of SCF and PSPEC slopes versus 
## outflow pdot surface density in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9.5,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.5)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'

    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co


    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="outflows_nro45m.csv",
    #                 return_cols="Pdot_flow", ra_col="RA_J2000", dec_col="DEC_J2000",
    #                 surface_density=True, feedback_label=r"$\dot P_{\rm out}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$ deg$^{-2}$)",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum", within_data=False,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))

    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))
    #         if irow < 2:
    #             ax.set_xlabel('')
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_outflow_pdot_surfacedensity{}.pdf".format(refit_str))
    # plt.close()

## Plotting scatter plots of SCF and PSPEC slopes versus 
## outflow pdot in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9.5,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.5)
    
    
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'

    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="outflows_nro45m.csv",
    #                 return_cols="Pdot_flow", ra_col="RA_J2000", dec_col="DEC_J2000",
    #                 surface_density=False, feedback_label=r"$\dot P_{\rm out}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$)",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum", within_data=False,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))

    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))
    #         if irow < 2:
    #             ax.set_xlabel('')
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_outflow_pdot{}.pdf".format(refit_str))
    # plt.close()

## Plotting scatter plots of SCF and PSPEC slopes versus 
## shell pdot in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #        figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
     
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #        range(len(mols)), axarr, mols):
    #    #iterate over each column.
    #    for icol, ax, stat_name in zip(
    #            range(len(stat_names)), axrow, stat_names):
     
    #        if stat_name == 'pspec':
    #            stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #            stat_name, stat_name, mol, pspec_str))
    #        else:
    #            stat_list = glob("{}_results/{}_{}*.p".format(
    #                stat_name, stat_name, mol))

    #        ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #        colors = np.array(['0.7'] * len(stat_list))
    #        colors[ii_davis] = '0.1'
     
    #        if mol == 'c18o':
    #            high_cut = high_cut_c18o 
    #        elif mol == '12co':
    #            high_cut = high_cut_12co
    #        elif mol == '13co':
    #            high_cut = high_cut_13co

    #        ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="shells_pdot_NtoS.csv",
    #                return_cols="Pdot_mid", ra_col="RA", dec_col="DEC",
    #                surface_density=False, feedback_label=r"$\dot P_{\rm shells}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$)",
    #                logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                scatter_kwargs={'color':colors, 's':30},
    #                feedback_mode="sum", within_data=False,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))
     
    #        ax.set_yticks(ytick_arr[icol,irow])
     
    #        mol_label = mols_dict[mol]
    #        ax.annotate(mol_label, (0.5, 0.9), 
    #                size=14, xycoords="axes fraction")
     
    #        if icol == 0:
    #            ax.yaxis.set_major_formatter(
    #                    ticker.StrMethodFormatter('{x:.2f}'))
    #        if icol == 1:
    #            # ax.set_ylabel('')
    #            ax.yaxis.set_major_formatter(
    #                    ticker.StrMethodFormatter('{x:.1f}'))
    #        if irow < 2:
    #            ax.set_xlabel('')
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("slope_shell_pdot{}.pdf".format(refit_str))
    # plt.close()

## Plot scatter plots of SCF and PSPEC slopes versus
## shell pdot surface density 
    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):
    #         if stat_name == 'pspec':
    #             stat_list = glob("{}_results/{}_{}*_{}*.p".format(
    #             stat_name, stat_name, mol, pspec_str))
    #         else:
    #             stat_list = glob("{}_results/{}_{}*.p".format(
    #                 stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'
    #         print("Running plot_stat_feedback")

    #         if mol == 'c18o':
    #             high_cut = high_cut_c18o 
    #         elif mol == '12co':
    #             high_cut = high_cut_12co
    #         elif mol == '13co':
    #             high_cut = high_cut_13co

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="shells_pdot_NtoS.csv",
    #                 return_cols="Pdot_mid", ra_col="RA", dec_col="DEC",
    #                 surface_density=True, feedback_label=r"$\dot P_{\rm shells}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$ deg$^{-2}$)",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum", within_data=False,
    #                 refit=stat_name is 'pspec',
    #                 refit_kwargs_list=[{"high_cut":high_cut}]*len(stat_list))
    #         print("finished plot_stat_feedback")
    #         ax.set_yticks(ytick_arr[icol,irow])

    #         mol_label = mols_dict[mol]
    #         ax.annotate(mol_label, (0.5, 0.9), 
    #                 size=14, xycoords="axes fraction")

    #         if icol == 0:
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.2f}'))
    #         if icol == 1:
    #             # ax.set_ylabel('')
    #             ax.yaxis.set_major_formatter(
    #                     ticker.StrMethodFormatter('{x:.1f}'))
    #         if irow < 2:
    #             ax.set_xlabel('')
    # # ax1.set_xlabel('')
    # # plt.tight_layout()
    # print("Saving figure.")
    # plt.savefig("slope_shell_pdot_density{}.pdf".format(refit_str))
    # plt.close()


## Plot covariance matrices.
## 12co
    # mols = ['12co', '13co', 'c18o']
    # for mol in mols:
    #     region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #             "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #             "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]
    #     pca_files = [glob("pca_results/pca_{}*{}*.p".format(mol, region))
    #             for region in region_names]

    #     f, axarr = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True,
    #             figsize=(10,8))
    #     f.subplots_adjust(hspace=0, wspace=0)

    #     pca_arr = np.reshape(pca_files, axarr.shape)
    #     region_arr = np.reshape(region_names, axarr.shape)

    #     #Iterate over each row.
    #     for irow, axrow in zip(
    #             range(axarr.shape[0]), axarr):
    #         #iterate over each column.
    #         for icol, ax in zip(
    #                 range(axarr.shape[1]), axrow):

    #             pca = pca_arr[irow, icol]

    #             ax = plot_cov(pca, ax=ax, imshow_kwargs={}, return_ax=True)

    #             ax.set_xlim([4,14.9])
    #             ax.set_ylim([4,14.9])
                

    #     #         mol_label = mols_dict[mol]
    #             ax.annotate(region_arr[irow,icol], (0.06, 0.9), 
    #                     size=14, xycoords="axes fraction",
    #                     color='white')
    #             ax.set_yticks([5,10])
    #             ax.set_xticks([5,10])

    #             #Remove tick mark.
    #             ax.tick_params(axis='both', which='both', length=0)
                
    #             # ax.grid(True, linestyle=':', color='gray')
    #             # ax.xaxis.set_ticks_position('both')
    #             # ax.yaxis.set_ticks_position('both')
    #             if icol > 0:
    #                 ax.set_ylabel('')
    #                 # ax.yaxis.set_major_formatter(
    #                 #         ticker.StrMethodFormatter('{x:.1f}'))
    #                 ax.tick_params(axis='y', which='both', length=4, direction='in')
    #             else:
    #                 ax.tick_params(axis='y', which='both', length=4)
    #             if irow < 2:
    #                 ax.set_xlabel('')
    #                 ax.tick_params(axis='x', which='both', length=4, direction='in')
    #             else:
    #                 ax.tick_params(axis='x', which='both', length=4)
    #     # # ax1.set_xlabel('')
    #     # # plt.tight_layout()
    #     plt.savefig("cov_{}.pdf".format(mol))
    #     plt.close()


#### Plot color-color plots of PCA distances between regions, ordered by 
#### the YSO surface density.
    
    # region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #         "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #         "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]

    # feedback_list = []

    # for region_name in region_names:

    #     feedback_table = region_stars(region_name=region_name,
    #             table="spitzer_orion.fit",
    #             return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
    #             ra_col="_RAJ2000", dec_col="_DEJ2000",
    #             subcube_dir="subcubes/", within_data=True)

    #     feedback_total = len(feedback_table.as_array())

    #     area = region_area(unit="deg", dist=414*u.pc,
    #         region_name=region_name, subcube_dir="subcubes/")

    #     feedback_total = feedback_total / area.value

    #     feedback_list.append(feedback_total)
    
    # ii_sort_nstars = np.argsort(feedback_list)
    # region_names = np.array(region_names)[ii_sort_nstars]
    # nstars = np.array(feedback_list)[ii_sort_nstars]

    # mols = ["12co", "13co", "c18o"]

    # mols_dict = {
    #     "12co":r"$^{12}$CO",

    # f, axarr = plt.subplots(nrows=1, ncols=3, sharey=True,
    #         figsize=(12,4))
    # f.subplots_adjust(wspace=0.)

    # for mol, ax in zip(mols,axarr):        

    #     pca_files = [glob("pca_results/pca_{}*{}*.p".format(mol, region))[0]
    #             for region in region_names]

    #     # print(pca_files)
       
    #     ax = plot_colorcolor(pca_files, ax=ax, return_ax=True,
    #             dist_func=dist_pca, dist_func_kwargs={},
    #             region_names=region_names, imshow_kwargs={})

    #     ax.set_xticklabels(region_names, size=10)
    #     ax.set_yticklabels(region_names, size=10)

    #     ax.set_title(mols_dict[mol], size=15)

    #     # ax1.set_xlabel('')
    # # plt.tight_layout()
    # plt.savefig("pca_colorcolor_sortedbynstar.pdf".format(mol), bbox_inches='tight')
    # plt.close()

### Plot SCF spectra of subregions together.
###
    # mols = ['12co', '13co', 'c18o']
    # yticks = [[0.4,0.6,0.8], [0.3,0.5,0.7], [0.1,0.2,0.3]]
    # # xticks = [[5,10,30], [5,10,30], [5,10,30]]

    # for mol,ytick in zip(mols,yticks):

    #     region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #             "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #             "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]

    #     bkgrd_colors = ["0.7"]*len(region_names)

    #     scf_files = [glob("scf_results/scf_{}*{}*.p".format(mol, region))[0]
    #             for region in region_names]

    #     f, axarr = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True,
    #             figsize=(10,8))
    #     f.subplots_adjust(hspace=0, wspace=0)

    #     scf_arr = np.reshape(scf_files, axarr.shape)
    #     region_arr = np.reshape(region_names, axarr.shape)

    #     #Iterate over each row.
    #     for irow, axrow in zip(
    #             range(axarr.shape[0]), axarr):
    #         #iterate over each column.
    #         for icol, ax in zip(
    #                 range(axarr.shape[1]), axrow):

    #             ii = scf_files.index(scf_arr[irow, icol])
    #             ax = plot_scf_spectrum(scf_files[:ii]+scf_files[ii+1:],
    #                     ax=ax, return_ax=True,
    #                     colors=bkgrd_colors[:ii] + bkgrd_colors[ii+1:],
    #                     labels=region_names[:ii] + region_names[ii+1:],
    #                     xlabel="Lag [pixels]", ylabel="SCF",
    #                     plot_points=True, plot_fit=True,
    #                     errorbar_kwargs={'ms':5, "elinewidth":2})

    #             #Plot one region's SCF spectrum in black.
    #             print(scf_arr[irow,icol])
    #             ax = plot_scf_spectrum([scf_files[ii]], ax=ax, return_ax=True,
    #                     colors=['black'], labels=[region_names[ii]],
    #                     xlabel="Lag [pixels]", ylabel="SCF",
    #                     plot_points=True, plot_fit=True, errorbar_kwargs={"zorder":3},
    #                     plot_kwargs={"zorder":3},
    #                     show_fitrange=False,
    #                     show_beamwidth=False)
    #             if mol == '12co': 
    #                 region_label_loc = (0.08,0.08) 
    #             if mol == '13co': 
    #                 region_label_loc = (0.08,0.08) 
    #             if mol == 'c18o': 
    #                 region_label_loc = (0.3,0.9) 
    #             ax.annotate(region_arr[irow,icol], region_label_loc, 
    #                     size=11, xycoords="axes fraction",
    #                     color='black')
    #             ax.set_xscale('log')
    #             ax.set_yscale('log')
                
    #             ax.set_xticks([5,10,30])
    #             ax.set_yticks(ytick)
                
    #             ax.get_xaxis().set_major_formatter(
    #                 ticker.ScalarFormatter())
                                                         
    #             ax.get_yaxis().set_major_formatter(
    #                 ticker.ScalarFormatter())

    #             ax.tick_params(axis='y', which='minor', labelleft=False)
              
    #             if icol > 0:
    #                 ax.set_ylabel('')
    #                 # ax.tick_params(axis='y', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='y', which='both', length=4)
    #             if irow < 2:
    #                 ax.set_xlabel('')
    #                 ax.tick_params(axis='x', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='x', which='both', length=4)
    #     plt.savefig("scf_spectrum_{}.pdf".format(mol))
    #     plt.close()

## Plot power spectra together.
##
    # mols = ['12co', '13co', 'c18o']
    # for mol in mols:

    #     region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #             "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #             "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]

    #     bkgrd_color = "0.7"

    #     pspec_files = [glob("pspec_results/pspec_{}*{}*.p".format(mol, region))[0]
    #             for region in region_names]

    #     f, axarr = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True,
    #             figsize=(10,8))
    #     f.subplots_adjust(hspace=0, wspace=0)

    #     pspec_arr = np.reshape(pspec_files, axarr.shape)
    #     region_arr = np.reshape(region_names, axarr.shape)

    #     #Iterate over each row.
    #     for irow, axrow in zip(
    #             range(axarr.shape[0]), axarr):
    #         #iterate over each column.
    #         for icol, ax in zip(
    #                 range(axarr.shape[1]), axrow):

    #             ii = pspec_files.index(pspec_arr[irow, icol])
    #             # if icol > 0: 
    #             #     break
    #             for i, pspec in enumerate(pspec_files):

    #                 if i == ii:
    #                     line_color = 'black'
    #                     fill_between = True
    #                     fb_color='tab:red'
    #                     scatter=True
    #                     sc_color='tab:red'
    #                     line=True
    #                     lw=3
    #                     zorder=3
    #                     if irow == 0:
    #                         twin_axis=True
    #                     else:
    #                         twin_axis=False
    #                 else:
    #                     line_color = 'gray'
    #                     fill_between = False
    #                     fb_color='tab:red'
    #                     scatter = False
    #                     sc_color='tab:red'
    #                     line = True
    #                     lw=1
    #                     zorder=1
    #                     twin_axis=False


    #                 ax = plot_pspec(pspec,
    #                         ax=ax, return_ax=True,
    #                         twin_axis=twin_axis,# twin_axis_labels=twin_axis_labels,
    #                         twin_axis_ticks=np.array([0.01,0.1,1]),
    #                         fill_between=fill_between,
    #                         scatter=scatter,
    #                         line=line,
    #                         fill_between_kwargs={
    #                             'zorder':zorder, 'alpha':0.3, 'color':fb_color},
    #                         scatter_kwargs={
    #                             'zorder':zorder, 's':7, 'color':sc_color, 'alpha':0.7},
    #                         plot_kwargs={
    #                             'color':line_color, 'zorder':zorder, 'linestyle':'-', 'lw':lw})
            
                    

    #             if mol == '12co': 
    #                 region_label_loc = (0.08,0.08) 
    #                 ax.set_ylim(10.1,23)
    #             if mol == '13co': 
    #                 region_label_loc = (0.08,0.08) 
    #                 ax.set_ylim(8,22)
    #             if mol == 'c18o': 
    #                 region_label_loc = (0.3,0.9)
    #                 ax.set_ylim(8,21)

                    
    #             ax.annotate(region_arr[irow,icol], region_label_loc, 
    #                     size=11, xycoords="axes fraction",
    #                     color='black')
                
    #             # ax.set_xscale('log')
    #             # ax.set_yscale('log')
                
    #             # ax.set_xticks([5,10,30])
    #             # ax.set_yticks(ytick)
                
    #             # ax.get_xaxis().set_major_formatter(
    #             #     ticker.ScalarFormatter())
                                                         
    #             # ax.get_yaxis().set_major_formatter(
    #             #     ticker.ScalarFormatter())

    #             ax.tick_params(axis='y', which='minor', labelleft=False)
              
    #             if icol > 0:
    #                 ax.set_ylabel('')
    #                 # ax.tick_params(axis='y', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='y', which='both', length=4)
    #             if irow < 2:
    #                 ax.set_xlabel('')
    #                 ax.tick_params(axis='x', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='x', which='both', length=4)
    #     plt.savefig("pspec_spectrum_{}.pdf".format(mol))
    #     plt.close()

## Plot power spectra together after beam correction and apodization.
## Refit the power spectra by restricting the fitting range.
    # mols = ['12co', '13co', 'c18o']
    # for mol in mols:

    #     region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #             "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #             "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]

    #     bkgrd_color = "0.7"

    #     pspec_files = [glob("pspec_results/pspec_{}*{}*beamcorrect_apodizetukey*.p".format(mol, region))[0]
    #             for region in region_names]

    #     f, axarr = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True,
    #             figsize=(10,8))
    #     f.subplots_adjust(hspace=0, wspace=0)

    #     pspec_arr = np.reshape(pspec_files, axarr.shape)
    #     region_arr = np.reshape(region_names, axarr.shape)

    #     #Iterate over each row.
    #     for irow, axrow in zip(
    #             range(axarr.shape[0]), axarr):
    #         #iterate over each column.
    #         for icol, ax in zip(
    #                 range(axarr.shape[1]), axrow):

    #             ii = pspec_files.index(pspec_arr[irow, icol])
    #             # if icol > 0: 
    #             #     break
    #             for i, pspec in enumerate(pspec_files):

    #                 if i == ii:
    #                     line_color = 'black'
    #                     fill_between = True
    #                     fb_color='tab:red'
    #                     scatter=True
    #                     sc_color='tab:red'
    #                     line=True
    #                     lw=3
    #                     zorder=3
    #                     if irow == 0:
    #                         twin_axis=True
    #                     else:
    #                         twin_axis=False
    #                     show_fitrange=True
    #                     show_beamwidth=True
    #                 else:
    #                     line_color = 'gray'
    #                     fill_between = False
    #                     fb_color='tab:red'
    #                     scatter = False
    #                     sc_color='tab:red'
    #                     line = True
    #                     lw=1
    #                     zorder=1
    #                     twin_axis=False
    #                     show_fitrange=False
    #                     show_beamwidth=False
    #                 if mol == 'c18o':
    #                     high_cut = high_cut_c18o 
    #                 elif mol == '12co':
    #                     high_cut = high_cut_12co
    #                 elif mol == '13co':
    #                     high_cut = high_cut_13co
    #                 print("Running plot_pspec on {}".format(mol))
    #                 ax = plot_pspec(pspec,
    #                         ax=ax, return_ax=True,
    #                         twin_axis=twin_axis,# twin_axis_labels=twin_axis_labels,
    #                         twin_axis_ticks=np.array([0.01,0.1,1]),
    #                         fill_between=fill_between,
    #                         scatter=scatter,
    #                         line=line,
    #                         fill_between_kwargs={
    #                             'zorder':zorder, 'alpha':0.3, 'color':fb_color},
    #                         scatter_kwargs={
    #                             'zorder':zorder, 's':7, 'color':sc_color, 'alpha':0.7},
    #                         plot_kwargs={
    #                             'color':line_color, 'zorder':zorder, 'linestyle':'-', 'lw':lw},
    #                         refit=True,
    #                         refit_kwargs={'high_cut':high_cut},
    #                         show_fitrange=show_fitrange,
    #                         show_beamwidth=show_beamwidth)
            
                    

    #             if mol == '12co': 
    #                 region_label_loc = (0.08,0.08) 
    #                 ax.set_ylim(10.1,23)
    #             if mol == '13co': 
    #                 region_label_loc = (0.08,0.08) 
    #                 ax.set_ylim(8,22)
    #             if mol == 'c18o': 
    #                 region_label_loc = (0.3,0.9)
    #                 ax.set_ylim(8,21)

                    
    #             ax.annotate(region_arr[irow,icol], region_label_loc, 
    #                     size=11, xycoords="axes fraction",
    #                     color='black')
                
    #             # ax.set_xscale('log')
    #             # ax.set_yscale('log')
                
    #             # ax.set_xticks([5,10,30])
    #             # ax.set_yticks(ytick)
                
    #             # ax.get_xaxis().set_major_formatter(
    #             #     ticker.ScalarFormatter())
                                                         
    #             # ax.get_yaxis().set_major_formatter(
    #             #     ticker.ScalarFormatter())

    #             ax.tick_params(axis='y', which='minor', labelleft=False)
              
    #             if icol > 0:
    #                 ax.set_ylabel('')
    #                 # ax.tick_params(axis='y', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='y', which='both', length=4)
    #             if irow < 2:
    #                 ax.set_xlabel('')
    #                 ax.tick_params(axis='x', which='both', length=4, direction='in')
    #             else:
    #                 pass
    #                 # ax.tick_params(axis='x', which='both', length=4)
    #     plt.savefig("pspec_spectrum_{}_beamcorrect_apodizetukey_highcut{}.pdf".format(mol,high_cut.value))
    #     plt.close()


## Plot power spectra of red/central/blue line wings in each jrf region
## after beam correction and apodization.
## 
## 


    region_names = ["jrf_north", "jrf_central", "jrf_south"]#, "jrf_l1641n"]
    mols = ['12co', '13co', 'c18o']
    colors = ['blue', 'black', 'red']
    for mol in mols:
        for region in region_names:
            

            pspec_files = glob(
                    "pspec_results/pspec_{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region))
            vlow_list = [float(p.split('_')[5].split('to')[0].replace('p','.')) for p in pspec_files]
            colors = np.array(['blue', 'black', 'red'])[np.argsort(vlow_list)]
            central_pspec = PowerSpectrum.load_results(np.array(pspec_files)[colors == 'black'][0])
            central_norm = central_pspec.ps1D[0]
            bmaj_pix = central_pspec.header['BMAJ'] / central_pspec.header['CDELT2']

            f, ax = plt.subplots(nrows=1, ncols=1)

            for pspec_file, color in zip(pspec_files, colors):
               ax = plot_pspec(pspec_file,
                            ax=ax, return_ax=True,
                            norm_factor=central_norm,
                            fill_between=True,
                            scatter=False,
                            line=False,
                            connect_points=True,
                            fill_between_kwargs={
                                'zorder':1, 'color':color, 'alpha':0.3},
                            # scatter_kwargs={
                            #     'zorder':3, 's':7, 'color':color, 'alpha':0.7},
                            connect_plot_kwargs={'color':color, 'linestyle':'-'},
                            refit=False,
                            show_fitrange=False,
                            show_beamwidth=False,
                            twin_axis=True,
                            twin_axis_ticks=np.array([0.01,0.1,1]))
            # plt.xlim(right=np.log10(1/bmaj_pix))

            ax.set_xlim(right=np.log10(1/bmaj_pix))
            ax.set_ylim(top=2, bottom=-8)
            print(ax.get_xlim())
      
      
    region = "jrf_l1641n"
    mols = ['12co', '13co', 'c18o']
    
    for mol in mols:
        red_l1641n_file = glob(
                    "pspec_results/pspec_{}*{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region, "10p81"))[0]
        red_ngc1999_file = glob(
                    "pspec_results/pspec_{}*{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region, "10p31"))[0]
        central_l1641n_file = glob(
                    "pspec_results/pspec_{}*{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region, "6p56"))[0]
        central_ngc1999_file = glob(
                    "pspec_results/pspec_{}*{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region, "8p56"))[0]
        blue_file = glob(
                    "pspec_results/pspec_{}*{}*{}*kms*beamcorrect_apodizetukey*.p".format(mol, region, "4p31"))[0]
        
        central_l1641n_pspec = PowerSpectrum.load_results(central_l1641n_file)
        central_norm = central_l1641n_pspec.ps1D[0]
        bmaj_pix = central_l1641n_pspec.header['BMAJ'] / central_l1641n_pspec.header['CDELT2']


        pspec_files = [red_l1641n_file, central_l1641n_file, blue_file]
        colors = ['red', 'black', 'blue']

        f, ax = plt.subplots(nrows=1, ncols=1)

        for pspec_file, color in zip(pspec_files, colors):

            ax = plot_pspec(pspec_file,
                ax=ax, return_ax=True,
                norm_factor=central_norm,
                fill_between=True,
                scatter=False,
                line=False,
                connect_points=True,
                fill_between_kwargs={
                    'zorder':1, 'color':color, 'alpha':0.3},
                # scatter_kwargs={
                #     'zorder':3, 's':7, 'color':color, 'alpha':0.7},h
                connect_plot_kwargs={'color':color, 'linestyle':'-'},
                refit=False,
                show_fitrange=False,
                show_beamwidth=False,
                twin_axis=True,
                twin_axis_ticks=np.array([0.01,0.1,1])) 


        ax.set_xlim(right=np.log10(1/bmaj_pix))
        ax.set_ylim(top=2, bottom=-8)
        print(ax.get_xlim())
        plt.tight_layout()
        plt.savefig("pspec_linewings_{}_{}_beamcorrect_apodizetukey.pdf".format(mol, region))
        plt.close()


        central_ngc1999_pspec = PowerSpectrum.load_results(central_ngc1999_file)
        central_norm = central_ngc1999_pspec.ps1D[0]

        pspec_files = [red_ngc1999_file, central_ngc1999_file, blue_file]
        colors = ['red', 'black', 'blue']

        f, ax = plt.subplots(nrows=1, ncols=1)

        for pspec_file, color in zip(pspec_files, colors):

            ax = plot_pspec(pspec_file,
                ax=ax, return_ax=True,
                norm_factor=central_norm,
                fill_between=True,
                scatter=False,
                line=False,
                connect_points=True,
                fill_between_kwargs={
                    'zorder':1, 'color':color, 'alpha':0.3},
                # scatter_kwargs={
                #     'zorder':3, 's':7, 'color':color, 'alpha':0.7},h
                connect_plot_kwargs={'color':color, 'linestyle':'-'},
                refit=False,
                show_fitrange=False,
                show_beamwidth=False,
                twin_axis=True,
                twin_axis_ticks=np.array([0.01,0.1,1])) 


        ax.set_xlim(right=np.log10(1/bmaj_pix))
        ax.set_ylim(top=2, bottom=-8)
        print(ax.get_xlim())
        plt.tight_layout()
        plt.savefig("pspec_linewings_{}_{}_ngc1999_beamcorrect_apodizetukey.pdf".format(mol, region))
        plt.close()


        

 ### Plot overview.
 ###
    # plt.clf()
    # plt.close()
    # cube = SpectralCube.read("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits")
    # wcs = WCS(cube.hdu).celestial
    # Tpeak = cube.max(axis=0).data[:]
    # print(Tpeak.shape)
    # fig = plt.figure(figsize=(6,5))
    # ax1 = plt.subplot(121, projection=wcs)
    # ax2 = plt.subplot(122, projection=wcs, sharey=ax1)
    
    # fig.subplots_adjust(wspace=0)

    # ax1 = plot_overview(Tpeak, wcs, ax=ax1, return_ax=True,
    #         imshow_kwargs={'origin':'lower', 'cmap':'viridis','interpolation':"None"})

    # ra1, dec1 = ax1.coords[0], ax1.coords[1]
    # ra1.set_ticks(number=3)
    # dec1.set_ticks_position('l')
   

    # ax2 = plot_overview(Tpeak, wcs, ax=ax2, return_ax=True, plot_regions=False,
    #         plot_ysos=True,
    #         yso_scatter_kwargs={'color':'tab:red', 'marker':'.', 's':2, 'edgecolors':'none', 'zorder':3},
    #         imshow_kwargs={'origin':'lower', 'cmap':'viridis', 'interpolation':"None"},
    #         plot_outflows=True,
    #         outflow_plot_kwargs={'markersize':7, 'linestyle':'None', 'color':'tab:olive', 'zorder':2},
    #         plot_shells=True,
    #         shell_patch_kwargs={'edgecolor':'white', 'facecolor':'none', 'zorder':1})
            
    # ax2.set_ylabel('')
    # ra2, dec2 = ax2.coords[0], ax2.coords[1]
    # dec2.set_ticklabel_visible(False)
    # dec2.set_ticks_position('r')
    # ra2.set_ticks(number=3)

    # ax1.set_xlim(-0.5, Tpeak.shape[1] - 0.5)
    # ax1.set_ylim(-0.5, Tpeak.shape[0] - 0.5)
    # ax2.set_xlim(-0.5, Tpeak.shape[1] - 0.5)
    # ax2.set_ylim(-0.5, Tpeak.shape[0] - 0.5)

    # plt.savefig("overview.pdf", bbox_inches='tight')
    # plt.clf()
    # plt.close()
    # ax1.cla()
    # ax2.cla()
    
    
    
    # ax1 = plot_overview(Tpeak, wcs, ax=ax1, return_ax=True)


    # ax1.imshow(Tpeak, origin='lower')
    # ax2.imshow(Tpeak, origin='lower')

    # fig.subplots_adjust(hspace=0, wspace=0)

    # plt.savefig("overview.pdf")
    


# Run VCS
#
    # print("Run VCS.")
    # cubes = glob("subcubes/*.fits")
    # names = [cube.split('/')[-1].split('.')[0] for cube in cubes]
    # pickle_files = ["vcs_results/vcs_"+name+".p" for name in names]
    
    # for i in range(len(cubes)):
    #     if not glob(pickle_files[i]+"*"):
    #         print("Running VCS on {}. Pickling to {}.".format(cubes[i], pickle_files[i]))
    #         run_vcs(cubes[i], distance=414*u.pc, keep_data=False,
    #                 run_kwargs={"xunit":u.m/u.s},
    #                 pickle_file=pickle_files[i])


## Plot VCS velocity power spectra together.
    # mols = ['12co', '13co', 'c18o']
    # for mol in mols:

    #     region_names = ["jrf_north", "jrf_central", "jrf_south", "jrf_l1641n",
    #             "davis_ngc1977", "davis_omc23", "ungerechts_omc1", "davis_omc4",
    #             "davis_omc5", "davis_hh34", "davis_l1641n", "davis_v380"]

    #     bkgrd_color = "0.7"

    #     vcs_files = [glob("vcs_results/vcs_{}*{}*.p".format(mol, region))[0]
    #             for region in region_names]
    #     print(vcs_files)
    #     f, axarr = plt.subplots(nrows=3, ncols=4, sharex=True, sharey=True,
    #             figsize=(10,8))
    #     f.subplots_adjust(hspace=0, wspace=0)

    #     vcs_arr = np.reshape(vcs_files, axarr.shape)
    #     region_arr = np.reshape(region_names, axarr.shape)

    #     #Iterate over each row.
    #     for irow, axrow in zip(
    #             range(axarr.shape[0]), axarr):
    #         #iterate over each column.
    #         for icol, ax in zip(
    #                 range(axarr.shape[1]), axrow):

    #             ii = vcs_files.index(vcs_arr[irow, icol])
    #             # if icol > 0: 
    #             #     break
    #             for i, vcs in enumerate(vcs_files):

    #                 if i == ii:
    #                     line_color = 'black'
    #                     plot_points=True
    #                     sc_color='tab:red'
    #                     plot_fit=True
    #                     lw=1
    #                     zorder=3
    #                     show_fitrange=True
    #                 else:
    #                     line_color = 'gray'
    #                     plot_points = False
    #                     sc_color='tab:red'
    #                     plot_fit = True
    #                     lw=1
    #                     zorder=1
    #                     show_fitrange=False
    #                 print("Running plot_vcs on {}".format(mol))
    #                 ax = plot_vcs_spectrum(vcs,
    #                         ax=ax, return_ax=True,
    #                         plot_points=plot_points,
    #                         plot_fit=plot_fit,
    #                         scatter_kwargs={
    #                             'zorder':zorder, 's':20, 'color':sc_color, 'alpha':1},
    #                         plot_kwargs={
    #                             'color':line_color, 'zorder':zorder, 'linestyle':'-', 'lw':lw},
    #                         axvline_kwargs={'color':line_color},
    #                         refit=False,
    #                         show_fitrange=show_fitrange)
            
                    

    #             if mol == '12co': 
    #                 region_label_loc = (0.27,0.9) 
    #                 # ax.set_ylim(10.1,23)
    #                 ax.set_xlim(3e-5,4e-3)
    #             if mol == '13co': 
    #                 region_label_loc = (0.27,0.9) 
    #                 # ax.set_ylim(8,22)
    #                 ax.set_xlim(3e-5,7e-3)
    #             if mol == 'c18o': 
    #                 region_label_loc = (0.29,0.9)
    #                 # ax.set_ylim(8,21)
    #                 ax.set_xlim(3e-5,7e-3)

                    
    #             ax.annotate(region_arr[irow,icol], region_label_loc, 
    #                     size=11, xycoords="axes fraction",
    #                     color='black')
                
    #             ax.set_xscale('log')
    #             ax.set_yscale('log')
                
        #         # ax.set_xticks([5,10,30])
        #         # ax.set_yticks(ytick)
                
        #         # ax.get_xaxis().set_major_formatter(
        #         #     ticker.ScalarFormatter())
                                                         
        #         # ax.get_yaxis().set_major_formatter(
        #         #     ticker.ScalarFormatter())

        #         ax.tick_params(axis='y', which='minor', labelleft=False)
              
        #         if icol > 0:
        #             ax.set_ylabel('')
        #             # ax.tick_params(axis='y', which='both', length=4, direction='in')
        #         else:
        #             pass
        #             # ax.tick_params(axis='y', which='both', length=4)
        #         if irow < 2:
        #             ax.set_xlabel('')
        #             ax.tick_params(axis='x', which='both', length=4, direction='in')
        #         else:
        #             pass
        #             # ax.tick_params(axis='x', which='both', length=4)
        # plt.savefig("vcs_spectrum_{}.pdf".format(mol))
        # plt.close()


    


######## OLD CODE

    ### Plot scatter plots of SCF and PSPEC slopes versus Pdot_shell in JRF regions.
    # stats = ['pspec','scf']
    # column = "Pdot_shell_mid"
    # column_low = "Pdot_shell_low"
    # column_hi = "Pdot_shell_hi"
    # for stat in stats:
    #     for mol in ['12co', '13co']:
    #         region_names = np.array(["Central", "L1641N", "North", "South_no_outliers"])
    #         stat_files = glob("{}_results/{}_{}*jrf*.p".format(stat, stat, mol))
    #         plot_scatter(stat_files, stat=stat,
    #             column=column, column_low=column_low, column_hi=column_hi,
    #             table="feedback_table.ecsv", 
    #             region_list=region_names,
    #             feedback_label="Shell Momentum Injection Rate",
    #             plot_file = "{}_{}_jrf_slope_{}.pdf".format(stat, mol, column))

    # stats = ['pspec','scf']
    # column = "Pdot_shell_mid_plus_out"
    # column_low = "Pdot_shell_low_plus_out"
    # column_hi = "Pdot_shell_hi_plus_out"
    # for stat in stats:
    #     for mol in ['12co', '13co']:
    #         region_names = np.array(["Central", "L1641N", "North", "South_no_outliers"])
    #         stat_files = glob("{}_results/{}_{}*jrf*.p".format(stat, stat, mol))
    #         plot_scatter(stat_files, stat=stat,
    #             column=column, column_low=column_low, column_hi=column_hi,
    #             table="feedback_table.ecsv", 
    #             region_list=region_names,
    #             feedback_label="Shell + Outflow Momentum Injection Rate",
    #             plot_file = "{}_{}_jrf_slope_{}.pdf".format(stat, mol, column),
    #             logy=True)

    # ### Plot scatter plots of SCF and PSPEC slopes versus Pdot_shell in JRF regions.
    # stats = ['pspec','scf']
    # column = "Pdot_shell_mid"
    # column_low = None 
    # column_hi = None
    # for stat in stats:
    #     for mol in ['12co', '13co']:
    #         region_names = np.array(["Central", "L1641N", "North", "South_no_outliers"])
    #         stat_files = glob("{}_results/{}_{}*jrf*.p".format(stat, stat, mol))
    #         plot_scatter(stat_files, stat=stat,
    #             column=column, column_low=column_low, column_hi=column_hi,
    #             table="feedback_table.ecsv", 
    #             region_list=region_names,
    #             feedback_label="Shell Momentum Injection Rate",
    #             plot_file = "{}_{}_jrf_slope_{}_noerrors.pdf".format(stat, mol, column))

    # stats = ['pspec','scf']
    # column = "Pdot_shell_mid_plus_out"
    # column_low = None 
    # column_hi = None
    # for stat in stats:
    #     for mol in ['12co', '13co']:
    #         region_names = np.array(["Central", "L1641N", "North", "South_no_outliers"])
    #         stat_files = glob("{}_results/{}_{}*jrf*.p".format(stat, stat, mol))
    #         plot_scatter(stat_files, stat=stat,
    #             column=column, column_low=column_low, column_hi=column_hi,
    #             table="feedback_table.ecsv", 
    #             region_list=region_names,
    #             feedback_label="Shell + Outflow Momentum Injection Rate",
    #             plot_file = "{}_{}_jrf_slope_{}_noerrors.pdf".format(stat, mol, column),
    #             logy=True)
    ### Plot color-color plots of only the JRF regions (removing outliers from the South, ordered by the
    ### median momentum injection rate of shells.
    #for mol in ['12co', '13co', 'c18o']:
        
    #     region_names = np.array(["North", "Central", "South_no_outliers", "L1641N"])
    #     pca_files = glob("pca_results/pca_{}*jrf*.p".format(mol))
    #     print(pca_files)
    #     l = [r.lower().split('_')[0] in p for p in pca_files for r in region_names]
    #     print(l)
    #     ii = np.where(l)[0] % len(pca_files)
    #     print(ii)
    #     #pca_files is in same order as region_names
    #     pca_files = np.array([pca_files[i] for i in ii])

    #     pdot_shell_mid = region_feedback(column="Pdot_shell_mid", region_list=region_names)
    #     pdot_shell_low = region_feedback(column="Pdot_shell_low", region_list=region_names)
    #     pdot_shell_hi = region_feedback(column="Pdot_shell_hi", region_list=region_names)
    #     print(pdot_shell_mid)
    #     ii_sort = np.argsort(pdot_shell_mid.data)
    #     print(pca_files[ii_sort])
    #     plot_colorcolor(pca_files[ii_sort], dist_func=dist_pca,
    #             region_names=region_names[ii_sort],
    #             plot_file="colorcolor_{}_pca_jrf_shellsort.pdf".format(mol))


#
    # pca_files = glob("pca_results/pca_12co*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pca_files]
    # plot_colorcolor(pca_files, dist_func=dist_pca, region_names=region_names,
    #         plot_file="colorcolor_12co_pca.pdf")

    # pca_files = glob("pca_results/pca_13co*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pca_files]
    # plot_colorcolor(pca_files, dist_func=dist_pca, region_names=region_names,
    #         plot_file="colorcolor_13co_pca.pdf")

    # pca_files = glob("pca_results/pca_c18o*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pca_files]
    # plot_colorcolor(pca_files, dist_func=dist_pca, region_names=region_names,
    #         plot_file="colorcolor_c18o_pca.pdf")

    # scf_files = glob("scf_results/scf_12co*[!pix].p.pkl")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in scf_files]
    # plot_colorcolor(scf_files, dist_func=dist_scf, region_names=region_names,
    #         plot_file="colorcolor_12co_scf.pdf")
    
    # scf_files = glob("scf_results/scf_13co*.pkl")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in scf_files]
    # plot_colorcolor(scf_files, dist_func=dist_scf, region_names=region_names,
    #         plot_file="colorcolor_13co_scf.pdf")

    # scf_files = glob("scf_results/scf_c18o*.pkl")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in scf_files]
    # plot_colorcolor(scf_files, dist_func=dist_scf, region_names=region_names,
    #         plot_file="colorcolor_c18o_scf.pdf")

    # pspec_files = glob("pspec_results/pspec_12co*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pspec_files]
    # plot_colorcolor(pspec_files, dist_func=dist_pspec, region_names=region_names,
    #         plot_file="colorcolor_12co_pspec.pdf")
    
    # pspec_files = glob("pspec_results/pspec_13co*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pspec_files]
    # plot_colorcolor(pspec_files, dist_func=dist_pspec, region_names=region_names,
    #         plot_file="colorcolor_13co_pspec.pdf")

    # pspec_files = glob("pspec_results/pspec_c18o*.p")
    # region_names = [
    #         p.split('.')[0].split('_')[-2] +
    #         "_" + p.split('.')[0].split('_')[-1] for p in pspec_files]
    # plot_colorcolor(pspec_files, dist_func=dist_pspec, region_names=region_names,
    #         plot_file="colorcolor_c18o_pspec.pdf")
    # pass
    # cubes = glob("subcubes/*c18o*.fits")

### Run SCF.
### Lags from 0 to 60" (coresponding to roughly the same physical scales as probed in 
### Boyden et al. 2015, or 0 to 0.12 pc at 414 pc distance).
### Each pixel is 2", and the smallest beam is ~6",
### so I'll go for lags separated by 3 pixels.

    # roll_lags = np.arange(-30, 31, 3)
    # names = [cube.split('/')[-1].split('.')[0] for cube in cubes]
    # pickle_files = ["scf_results/scf_"+name+".p" for name in names]
    # for i in range(len(cubes)):
    #     if not glob(pickle_files[i]+"*"):
    #         print("Running SCF on {}. Pickling to {}.".format(cubes[i], pickle_files[i]))
    #         run_scf(cubes[i], distance=None,
    #                 xunit=None, roll_lags=roll_lags,
    #                 pickle_file=pickle_files[i])

### Test SCF with different max lags.
###
    # cube = "subcubes/12co_ungerechts_omc1.fits"
    # name = cube.split('/')[-1].split('.')[0]
    # max_lags = [15,30,45,60]
    # scf_files=[]
    # for m in max_lags:
    #     pickle_file = "scf_results/scf_"+name+"_"+str(m)+"pix"
    #     run_scf(cube, distance=None,
    #             xunit=None, roll_lags=np.arange(-m, m+1, 3),
    #             pickle_file=pickle_file)
   
    # plot_scf_spectrum(glob("scf_results/scf_12co*pix.p.pkl"),
    #         "scf_ungerechts_omc1_maxlagtest.pdf",
    #         colors=["red", "blue", "green", "brown"],
    #         labels=["15 pix", "30 pix", "45 pix", "60 pix"],
    #         scale_lags=1)


### Run PCA. 
###
    # cubes = glob("subcubes/*.fits")
    # names = [cube.split('/')[-1].split('.')[0] for cube in cubes]
    # pickle_files = ["pca_results/pca_"+name+".p" for name in names]
   
    # for i in range(len(cubes)):
    #     if not glob(pickle_files[i]+"*"):
    #         print("Running PCA on {}. Pickling to {}.".format(cubes[i], pickle_files[i]))
    #         run_pca(cubes[i], distance=414*u.pc,
    #                 min_eigval=1e-4, spatial_output_unit=u.pc,
    #                 spectral_output_unit=u.m/u.s, brunt_beamcorrect=True,
    #                 pickle_file=pickle_files[i])

     # run_pca("12co_ungerechts_omc1.fits", distance=414*u.pc,
     #         pickle_file="pca_12co_ungerechts_omc1.p",
     #         min_eigval=1e-4, spatial_output_unit=u.pc,
     #         spectral_output_unit=u.m/u.s,

### Run Power Spectrum
###

    # cubes = glob("subcubes/*.fits")
    # names = [cube.split('/')[-1].split('.')[0] for cube in cubes]
    # pickle_files = ["pspec_results/pspec_"+name+"_beamcorrect_apodizetukey"+".p" for name in names]
    
    # for i in range(len(cubes)):
    #     if not glob(pickle_files[i]+"*"):
    #         print("Running Power Spectrum on {}. Pickling to {}.".format(cubes[i], pickle_files[i]))
    #         run_pspec(cubes[i], distance=414*u.pc, keep_data=False,
    #                 run_kwargs={'xunit':u.pix**-1, 'beam_correct':True,
    #                     'apodize_kernel':'tukey'},
    #                 pickle_file=pickle_files[i])

### Run power spectra for velocity ranges specifed in velocity_ranges.txt. These
#are the red/central/blue line wings for different parts of Orion A, to compare to the
#power spectra in Swift and Welch 2008.
    # vel_range_file = "velocity_ranges.txt"
    # vel_range_tab = ascii.read(vel_range_file)

    # distance=414*u.pc
    # run_kwargs={'xunit':u.pix**-1, 'beam_correct':True, 'apodize_kernel':"tukey"}

    # for region, blue_range, central_range, red_range in vel_range_tab:
    #     for vel_range in [blue_range, central_range, red_range]:
    #         vel_range_str = vel_range.strip('[]').replace(
    #                                      '.','p').replace(
    #                                      ',','to') + 'kms'
    #         vel_range = [float(s)*u.km/u.s
    #                 for s in vel_range.strip('[]').split(',')]
    #         print(vel_range_str)
    #         print(vel_range)
            
    #         for mol in mols:
    #             cube_file = "subcubes/{}_{}.fits".format(mol, region)
    #             print(cube_file)

    #             mom0 = calc_mom0(cube_file, vel_range=vel_range)
    #             pspec = PowerSpectrum(mom0, distance=distance)
    #             pspec.run(**run_kwargs)
    #             pspec.save_results(
    #                     "pspec_results/pspec_{}_{}_{}_beamcorrect_apodizetukey".format(
    #                         mol, region, vel_range_str, keep_data=False))


### Plot Power Spectra
### 
    # pspec_files = glob("pspec_results/*.p")
    # for f in pspec_files:
    #     plot_file = f.split('.')[0] + ".pdf"
    #     plot_pspec(f, plot_file)
        
            
# ## Plot covariance matrices.
# ##
#     pca_files = glob("pca_results/*.p")
#     for pca_file in pca_files:
#         print(pca_file)
#         plot_file = pca_file.split('.')[0] + "_cov.pdf"
#         try:
#             plot_cov(pca_file, plot_file)
#         except:
#             pass

# #### Plot SCFs with fits.
#     scf_files = glob("scf_results/*.pkl")
#     region_names = list(set([scf_file.split('o_')[-1].split('.')[0] for scf_file in scf_files]))
#     print(region_names) 
#     for region in region_names:
#         #print("Plotting from {}.".format(scf_file))
#         plot_file = "scf_results/scf_{}_spectrumfit.pdf".format(region)
#         print(plot_file)
#         scf_files = glob("scf_results/*{}*.pkl".format(region))
#         scf_files.sort()
#         print(scf_files)
#         #Scale lags from pixels to degrees at 2"/pixel
#         plot_scf_spectrum(scf_files, plot_file, colors=['red','blue'],
#                 labels=[r"$^{12}$CO", r"$^{13}$CO"],
#                 scale_lags = 1)
        
# ### Plot all Davis region SCF fits on the same plot.
#     mol = "12co"
    
#     scf_files = glob("scf_results/*{}*davis*.pkl".format(mol))
#     plot_file = "scf_results/scf_{}_spectrumfit.pdf".format(mol)
#     regions = [f.split("davis_")[1].split(".")[0] for f in scf_files]
    
#     cmap = matplotlib.cm.get_cmap('Dark2')
#     colors = cmap.colors[:len(regions)] 
    
#     plot_scf_spectrum(scf_files, plot_file, colors=colors,
#             labels=regions)

#     mol = "13co"

#     scf_files = glob("scf_results/*{}*davis*.pkl".format(mol))
#     plot_file = "scf_results/scf_{}_spectrumfit.pdf".format(mol)
#     regions = [f.split("davis_")[1].split(".")[0] for f in scf_files]
    
#     cmap = matplotlib.cm.get_cmap('Dark2')
#     colors = cmap.colors[:len(regions)] 
    
#     plot_scf_spectrum(scf_files, plot_file, colors=colors,
#             labels=regions)

#     mol = "c18o"

#     scf_files = glob("scf_results/*{}*davis*.pkl".format(mol))
#     plot_file = "scf_results/scf_{}_spectrumfit.pdf".format(mol)
#     regions = [f.split("davis_")[1].split(".")[0] for f in scf_files]
    
#     cmap = matplotlib.cm.get_cmap('Dark2')
#     colors = cmap.colors[:len(regions)] 
    
#     plot_scf_spectrum(scf_files, plot_file, colors=colors,
#labels=regions)

## Plot all JRF region SCF fits on the same plot.

    # mol = "12co"
    
    # scf_files = glob("scf_results/*{}*jrf*.p".format(mol))
    # scf_files.sort()
    # plot_file = "scf_results/scf_jrf_{}_spectrumfit.pdf".format(mol)
    # regions = [f.split("jrf_")[1].split(".")[0] for f in scf_files]
    
    # cmap = matplotlib.cm.get_cmap('Dark2')
    # colors = cmap.colors[:len(regions)] 
    
    # plot_scf_spectrum(scf_files, plot_file, colors=colors,
    #         labels=regions)

    # mol = "13co"

    # scf_files = glob("scf_results/*{}*jrf*.p".format(mol))
    # scf_files.sort()
    # plot_file = "scf_results/scf_jrf_{}_spectrumfit.pdf".format(mol)
    # regions = [f.split("jrf_")[1].split(".")[0] for f in scf_files]
    
    # cmap = matplotlib.cm.get_cmap('Dark2')
    # colors = cmap.colors[:len(regions)] 
    
    # plot_scf_spectrum(scf_files, plot_file, colors=colors,
    #         labels=regions)

    # mol = "c18o"

    # scf_files = glob("scf_results/*{}*jrf*.p".format(mol))
    # scf_files.sort()
    # plot_file = "scf_results/scf_jrf_{}_spectrumfit.pdf".format(mol)
    # regions = [f.split("jrf_")[1].split(".")[0] for f in scf_files]
    
    # cmap = matplotlib.cm.get_cmap('Dark2')
    # colors = cmap.colors[:len(regions)] 
    
    # plot_scf_spectrum(scf_files, plot_file, colors=colors,
    #     labels=regions)


def region_area(unit="deg", dist=414*u.pc,
        region_name="davis_omc23", subcube_dir="subcubes/"):
    """
    returns astropy quantity
    """
    from astropy.wcs.utils import proj_plane_pixel_area
    cube = SpectralCube.read(
            glob("{}*{}*".format(subcube_dir, region_name))[0])
    pixel_area = proj_plane_pixel_area(cube.wcs)*u.deg**2.
    npixels = np.count_nonzero(~np.isnan(cube[0]))
    area = pixel_area * npixels 
    return area.to(u.Unit(unit)*u.Unit(unit)) 

def region_stars(return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
        ra_col="_RAJ2000", dec_col="_DEJ2000",
        region_name="davis_omc23",
        table="spitzer_orion.fit", subcube_dir="subcubes/",
        within_data=False):
    """
    Return the rows from table which fall within the region specified
    by the subcube in subcube_dir with the name given by region_name.
    Returns an astropy.table with columns specified in return_cols.
    """
    from astropy.table import Table
    from glob import glob
    from spectral_cube import SpectralCube
    table = Table.read(table)
    cube = SpectralCube.read(
            glob("{}*{}*".format(subcube_dir, region_name))[0])
    ii = np.where((table[ra_col].data*u.deg > cube.longitude_extrema[0]) &
            (table[ra_col].data*u.deg < cube.longitude_extrema[1]) &
            (table[dec_col].data*u.deg > cube.latitude_extrema[0]) &
            (table[dec_col].data*u.deg < cube.latitude_extrema[1]))
    
    if within_data:
        ra, dec = table[ra_col][ii], table[dec_col][ii]

        coord = SkyCoord(ra, dec, frame='fk5', unit='deg')
        xy = astropy.wcs.utils.skycoord_to_pixel(coord, WCS(cube.hdu).celestial)
        xy_round = np.round(xy).astype(np.int)

        
        data = np.array(cube[0])

        xy_round[1][xy_round[1] == data.shape[0]] = data.shape[0] - 1
        xy_round[0][xy_round[0] == data.shape[1]] = data.shape[1] - 1

        xy_chan0 = data[xy_round[1], xy_round[0]]
        chan0_not_nan = ~np.isnan(xy_chan0)  
        ra,dec = ra[chan0_not_nan], dec[chan0_not_nan]
    
        return table[ii][chan0_not_nan][return_cols]

    else:
        return table[ii][return_cols]

# def region_shellimpact(column='Pdot_mid',
#         table="shells_pdot_NtoS.csv", ra_col="RA", dec_col="DEC",
#         region_name="davis_omc23", subcube_dir="subcubes/"):
#     """
#     Returns the total impact from shells in the specified region. 
#     """
#     from astropy.io import ascii
#     table = ascii.read(table)

#     region_shells = region_stars(return_cols=[column], ra_col=ra_col, dec_col=dec_col,
#             region_name=region_name, subcube_dir=subcube_dir)

#     total = region_shells[column].sum()
    
#     return total


# stat_list, stat="scf", column="Pdot_shell_mid", column_low="Pdot_shell_low", column_hi="Pdot_shell_hi",
#         table="feedback_table.ecsv", 
#         region_list=["North", "Central", "South", "L1641N"],
#         feedback_label="Shell Momentum Injection Rate",
#         plot_file = "scf_slope_pdotshell.pdf", mpl_style='presentation', logy=False


def plot_overview(data, wcs=None, plotfile="plot_overview.pdf",
        ax=None, return_ax=True,
        plot_data=True, imshow_kwargs={'origin':'lower', 'cmap':'viridis'}, nan_color='black',
        plot_regions=True, region_radec=[subregion.regions_ra, subregion.regions_dec],
        region_names=subregion.regions_name, 
        region_rectangle_kwargs={'edgecolor':'tab:red','facecolor':'none'},
        small_region_color='white',
        plot_ysos=False, yso_file="spitzer_orion.fit",
        yso_scatter_kwargs={},
        plot_shells=False, shell_file="../../shells/shell_candidates/AllShells_NtoS.reg",
        shell_patch_kwargs={'edgecolor':'white', 'facecolor':'none'},
        plot_outflows=False, outflow_file='outflows_nro45m.csv',
        outflow_plot_kwargs={'markersize':20, 'linestyle':'None', 'color':'white'},
        mpl_style="presentation",
        xlabel="RA [J2000]", ylabel="DEC [J2000]"):
    from matplotlib.patches import Rectangle
    from astropy.table import Table
    matplotlib.style.use(mpl_style)
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111, projection=wcs)
   
    if plot_data:
        current_cmap = matplotlib.cm.get_cmap(imshow_kwargs['cmap'])
        current_cmap.set_bad(color=nan_color)
        ax.imshow(data, **imshow_kwargs)

    if plot_regions:
        for ra, dec, name in zip(region_radec[0], region_radec[1], region_names):
            fontsize = 9
            if "davis" in name or "unger" in name:
                r = Rectangle(
                        (ra[0].value, dec[0].value),
                        ra[1].value - ra[0].value, dec[1].value - dec[0].value,
                        transform=ax.get_transform('fk5'),
                        edgecolor=small_region_color, facecolor='none'
                        )

                s = name.split('_')[1].upper()
                
                if "ngc1977" in name:
                    xy = (ra[0].value+0.023, 0.5*(dec[1] + dec[0]).value - 0.1)
                    va,ha = 'center', 'right'

                elif "omc1" in name:
                    xy = (ra[1].value-0.04, 0.5*(dec[1] + dec[0]).value)
                    va,ha = 'center', 'left'

                else:
                    xy = (0.5*(ra[0] + ra[1]).value, 0.5*(dec[1] + dec[0]).value)
                    va,ha = 'center', 'center'

                ax.annotate(
                        s,
                        xy,
                        xycoords=ax.get_transform('fk5'),
                        color=small_region_color,
                        verticalalignment=va,
                        horizontalalignment=ha,
                        fontsize=fontsize)
                
            else: 
                r = Rectangle(
                        (ra[0].value, dec[0].value),
                        ra[1].value - ra[0].value, dec[1].value - dec[0].value,
                        transform=ax.get_transform('fk5'),
                        **region_rectangle_kwargs)
                if '1641' in name:
                    s = name.split('_')[1].upper()
                    xy = (ra[0].value - 0.06, 0.5*(dec[1] + dec[0]).value - 0.05) 
                else:
                    s = name.split('_')[1].capitalize()
                    xy = (ra[0].value - 0.06, 0.5*(dec[1] + dec[0]).value) 
                ax.annotate(
                        s,
                        xy,
                        xycoords=ax.get_transform('fk5'),
                        color=region_rectangle_kwargs['edgecolor'],
                        verticalalignment='center',
                        fontsize=fontsize)
            ax.add_patch(r)

    if plot_ysos:
        t = Table.read("spitzer_orion.fit")
        yso_ra, yso_dec = t['_RAJ2000'], t['_DEJ2000']

        ra_lim = [min(wcs.calc_footprint()[:,0]), max(wcs.calc_footprint()[:,0])]
        dec_lim = [min(wcs.calc_footprint()[:,1]), max(wcs.calc_footprint()[:,1])]

        in_map = (yso_ra > ra_lim[0]) & (yso_ra < ra_lim[1]) &\
                 (yso_dec > dec_lim[0]) & (yso_dec < dec_lim[1])
        
        yso_ra, yso_dec = yso_ra[in_map], yso_dec[in_map]

        yso_coord = SkyCoord(yso_ra, yso_dec, frame='fk5', unit='deg')
        yso_xy = astropy.wcs.utils.skycoord_to_pixel(yso_coord, wcs)
        yso_xy_round = np.round(yso_xy).astype(np.int)


        yso_xy_round[1][yso_xy_round[1] == data.shape[0]] = data.shape[0] - 1
        yso_xy_round[0][yso_xy_round[0] == data.shape[1]] = data.shape[1] - 1
        print(yso_xy_round.shape)
        yso_xy_Tmb = np.array(data)[yso_xy_round[1], yso_xy_round[0]]
        data_not_nan = ~np.isnan(yso_xy_Tmb)  
        yso_ra,yso_dec = yso_ra[data_not_nan], yso_dec[data_not_nan]

        ax.scatter(yso_ra, yso_dec,
                transform=ax.get_transform('fk5'), **yso_scatter_kwargs)

    if plot_outflows:
        t = Table.read("outflows_nro45m.csv")
        ra_outflows, dec_outflows, pa_outflows = t["RA_J2000"], t["DEC_J2000"], t["PA"] 
        for ra, dec, pa in zip(ra_outflows, dec_outflows, pa_outflows):
            ax.plot(ra, dec, transform=ax.get_transform('fk5'),
                    marker=(2, 0, 180+pa), **outflow_plot_kwargs)


    if plot_shells:
        from astropy.visualization.wcsaxes import SphericalCircle
        import pyregion
        region_list = pyregion.open(shell_file)
        for region in region_list:
            
            circle = SphericalCircle(
                    (region.coord_list[0] * u.deg, region.coord_list[1] * u.deg),
                     region.coord_list[2] * u.deg,
                     transform=ax.get_transform('fk5'),
                     **shell_patch_kwargs)

            ax.add_patch(circle)
      

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if return_ax:
        return ax

    else:
        plt.tight_layout()
        plt.savefig(plotfile)


def plot_stat_feedback(stat_list, ax=None, catalog="spitzer_orion.fit",
        return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
        ra_col="_RAJ2000", dec_col="_DEJ2000",
        region_names=["davis_omc23", "davis_ngc1977"],
        regions_from_stat_list=True,
        subcube_dir = "subcubes/", feedback_mode='count', within_data=False,
        errorbar_kwargs={}, scatter_kwargs={},
        surface_density=True, feedback_label=r"n$_{\rm YSO}$ (deg$^{-2}$)",
        logx=False, logy=False, mpl_style="presentation",
        return_ax=True,return_eb=False, plot_file="plot_stat_feedback.pdf",
        refit=False, refit_kwargs_list=[{}]):
    """
    Make scatter plot of subregion stat metric
    vs. feedback metric:
    1. number (or density of) stars from a catalog of ysos.
    2. total (or density) momentum injection rate from shells.
    3. total (or density) momentum injection rate from outflows.

    stat_list and region_names should be in the same order,
    or can automatically find region_names from stat_list if
    regions_from_stat_list = True
    """
    # import subregion 
    matplotlib.style.use(mpl_style)

    if ax is None:
        
        fig = plt.figure(figsize=(6,6))
        ax = plt.subplot(111)
    
    if "scf" in stat_list[0]:
        stat_measure, stat_measure_err = slope_scf(stat_list)
        stat_label = "SCF Slope"
    elif "pspec" in stat_list[0]:
        stat_measure, stat_measure_err = slope_pspec(stat_list, refit=refit,
                refit_kwargs_list=refit_kwargs_list)
        stat_label = "SPS Slope"
    else:
        raise(Exception("stat must be one of ['scf', 'pspec']"))
    
    if regions_from_stat_list:
        region_names_choices = subregion.regions_name
        region_names = []
        for stat_file in stat_list:
            for region_name in region_names_choices:
                if region_name in stat_file:
                   region_names.append(region_name) 
                
    feedback_list = []
    print("Looping over regions.\n")
    for region_name in region_names:
        print("Now doing region {}.\n".format(region_name))
        

        feedback_table = region_stars(region_name=region_name, table=catalog,
                return_cols=return_cols, ra_col=ra_col,
                dec_col=dec_col, subcube_dir=subcube_dir, within_data=within_data)

        if feedback_mode == 'count':
            feedback_total = len(feedback_table.as_array())
            
        elif feedback_mode == 'sum':
            print(feedback_table)
            feedback_total = feedback_table.sum()
        
        elif feedback_mode == 'area':
            area = region_area(unit="deg", dist=414*u.pc,
                region_name=region_name, subcube_dir="subcubes/")
            feedback_total = area.value

        if surface_density:
            print("Calculating surface density.\n")
            area = region_area(unit="deg", dist=414*u.pc,
                region_name=region_name, subcube_dir="subcubes/")

            feedback_total = feedback_total / area.value

        feedback_list.append(feedback_total)

    
    print("Plotting.\n")
    eb = ax.errorbar(feedback_list, stat_measure, yerr=stat_measure_err, fmt="None",
            **errorbar_kwargs)
    sc = ax.scatter(feedback_list, stat_measure, **scatter_kwargs)

    if logy:
        ax.set_yscale('log')

    if logx:
        ax.set_xscale('log')
        
   
    xla = ax.set_xlabel("{}".format(feedback_label))
    yla = ax.set_ylabel(stat_label)

    
    if return_ax:
        if return_eb:
            return ax,eb
        else:
            return ax
        
    else: 
        plt.tight_layout()
        plt.savefig(plot_file)


# def region_feedback(column="Pdot_shell_mid", table="feedback_table.ecsv",
#         region_list=["North", "Central", "South", "South_no_outliers", "L1641N"],
#         return_region_names=False):
#     """
#     Returns the feedback impact from feedback_table for the regions in region_list.
#     """
#     from astropy.io import ascii
#     table = ascii.read(table)
#     table = table[np.where(np.isin(table["Region_Name"], region_list))]
#     if return_region_names:
#         return table["Region_Name", column]
#     else:
#         return table[column]


### Functions for running turbustat programs.
def run_pca(cube, distance=414*u.pc, pickle_file=None, **kwargs):
    from turbustat.statistics import PCA
    cube = SpectralCube.read(cube)
    pca = PCA(cube, distance=distance)
    pca.run(**kwargs)
    
    if pickle_file:
        pickle.dump(pca, open(pickle_file, "wb"))    
    return pca

def run_scf(cube, size=None, roll_lags=np.arange(-15,15+1,3),
        distance=414*u.pc, xunit=u.pc, pickle_file=None, keep_data=True, **kwargs):
    from turbustat.statistics import SCF
    cube = SpectralCube.read(cube)
    scf = SCF(cube, size=size, roll_lags=roll_lags, distance=distance)
    scf.run(boundary='cut', xunit=xunit, **kwargs)
    if pickle_file:   
        scf.save_results(pickle_file, keep_data=keep_data)

def run_pspec(cube, distance=414*u.pc,
        run_kwargs={'xunit':u.pix**-1, 'beam_correct':False, 'apodize_kernel':"tukey"},
        pickle_file=None, keep_data=True):
    from turbustat.statistics import PowerSpectrum
    cube = SpectralCube.read(cube)
    mom0_hdu = cube.moment0().hdu
    pspec = PowerSpectrum(mom0_hdu, distance=distance)

    pspec.run(**run_kwargs)
    
    try:
        pspec.save_results(pickle_file, keep_data=keep_data)
    except AttributeError:
        pass

    try:
        os.rename("{}.pkl".format(pickle_file), pickle_file)
    except FileNotFoundError:
        pass

    return pspec

def run_vcs(cube, distance=414*u.pc,
        run_kwargs={},
        pickle_file=None, keep_data=True):
    from turbustat.statistics import VCS
    cube = SpectralCube.read(cube)
    vcs = VCS(cube)
    vcs.run(**run_kwargs)
     
    try:
        vcs.save_results(pickle_file, keep_data=keep_data)
    except AttributeError:
        pass

    try:
        os.rename("{}.pkl".format(pickle_file), pickle_file)
    except FileNotFoundError:
        pass

    return vcs   
    



### Functions for calculating distance metrics of stats between subregions.
def dist_pca(pca_list, remove_first_eigval=True, normalize=True):
    """
    Calculate the distance metric each pair of a set of PCA results.
    Returns an array that can be plotted like Boyden et al. 2016 Figure 16.
    """

    ## If mean subtraction was not done, do not use the first eigenvalue.
    print(type(pca_list[0]))
    try:
        print("Unpickling pca files")
        pca_list = [pickle.load(open(f, 'rb')) for f in pca_list]
    except:
        pass
    
    if remove_first_eigval:
        ind = 1
    else:
        ind = 0

    if normalize:
        eigval_norm_list = [pca.eigvals[ind:]/np.sum(pca.eigvals[ind:]) for pca in pca_list]
    else:
        eigval_norm_list = [pca.eigvals[ind:] for pca in pca_list]
    # add a phantom axis to make a row (a) and column (b) vector with each
    # subregion's eigenvalue array in each element.
    
    a = np.expand_dims(np.array(eigval_norm_list),1)
    b = np.expand_dims(np.array(eigval_norm_list),0)
    
    dist_array = np.linalg.norm(a-b, axis=2)
    
    return dist_array

def dist_scf(scf_list, weight=True, weight_correct=True):
    """
    weight_correct = False corresponds to the 
    weighting scheme in Turbustat:
    http://turbustat.readthedocs.io/en/latest/_modules/turbustat/statistics/scf/scf.html
    in the distance_metric method of the SCF_Distance class.

    weight_correct = True calculates the weights of each 
    pixel in the scf surface using their center-to-center 
    distance from the central pixel.
    """
    try:
        scf_list = [pickle.load(open(f, 'rb')) for f in scf_list]
    except:
        pass
    
    scf0 = scf_list[0]
    dx = np.arange(scf0.size) - scf0.size // 2
    dy = dx
    #Pixel coordinates of the scf surface.
    a,b = np.meshgrid(dx, dy)

    if weight:
        if weight_correct:
            dist_weight = 1 / np.sqrt(a ** 2 + b ** 2)
            # Centre pixel set to 1
            dist_weight[np.where((a == 0) & (b == 0))] = 1.

        else:
            a[np.where(a == 0)] = 1.
            b[np.where(b == 0)] = 1.
            dist_weight = 1 / np.sqrt(a ** 2 + b ** 2)
    else:
        dist_weight = np.ones((scf0.size, scf0.size))

    scf_surface_list = np.array([scf.scf_surface for scf in scf_list])
    
    #Add a phantom axis to make a row and column vector with each subregions
    #scf_surface array in each element
    scf_surface_row = np.expand_dims(np.array(scf_surface_list),0)
    scf_surface_col = np.expand_dims(np.array(scf_surface_list),1)
    
    scf_surface_diff_array = (
            scf_surface_row - scf_surface_col) ** 2. * dist_weight
    
    dist_array = np.sqrt(
            np.sum(scf_surface_diff_array, axis=(2,3)) / np.sum(dist_weight))
    
    return dist_array

def slope_scf(scf_list):
    try:
        scf_list = [pickle.load(open(f, 'rb')) for f in scf_list]
    except:
        pass
    
    slope = np.array([scf.fit.params[1] for scf in scf_list])
    slope_err = np.array([scf.fit.bse[1] for scf in scf_list])

    return slope, slope_err

def dist_pspec(pspec_list):
    try:
        pspec_list = [pickle.load(open(f, 'rb')) for f in pspec_list]
    except:
        pass
    
    slope = np.array([pspec.slope for pspec in pspec_list])
    slope_err = np.array([pspec.slope_err for pspec in pspec_list])

    slope_row = np.expand_dims(np.array(slope),0)
    slope_col = np.expand_dims(np.array(slope),1)

    slope_err_row = np.expand_dims(np.array(slope_err),0)
    slope_err_col = np.expand_dims(np.array(slope_err),1)
    
    dist_array = abs(
            (slope_row - slope_col) / np.sqrt(
                slope_err_row ** 2 + slope_err_col ** 2 ))
    
    return dist_array

def slope_pspec(pspec_list, refit=False, refit_kwargs_list=[{}]):
    
    try:
        pspec_list = [pickle.load(open(f, 'rb')) for f in pspec_list]
    except:
        pass

    if refit:       
        for pspec,refit_kwargs in zip(pspec_list, refit_kwargs_list):
            pspec.fit_pspec(**refit_kwargs)
    

    slope = np.array([pspec.slope for pspec in pspec_list])
    slope_err = np.array([pspec.slope_err for pspec in pspec_list])

    return slope, slope_err

### Functions plotting color-color plots of statistical distance metrics
### between subregions of data.
def plot_colorcolor(stat_list, ax=None, return_ax=True,
        dist_func=dist_pca, dist_func_kwargs={},
        region_names=None, mpl_style="presentation",
        imshow_kwargs={},
        savefig_kwargs={"dpi":300}, plot_file="color_color.pdf"):
    """
    stat_list: can be a list of str filenames pointing to pickled results
    of the stat method corresponding to dist_func,
        or a list of turbustat.statistics objects with the results.
    dist_func_args: optional, dict of keyword arguments to dist_func.
    """
    matplotlib.style.use(mpl_style)
    dist_array = dist_func(stat_list, **dist_func_kwargs)

    if return_ax is None:
        fig, ax = plt.subplots(1, figsize=(8,8))
    
    im = ax.imshow(dist_array, origin='lower', **imshow_kwargs)
    # xla = ax.set_xlabel("subregion") 
    # yla = ax.set_ylabel("subregion")

    xticks = ax.set_xticks(np.arange(len(stat_list)))
    yticks = ax.set_yticks(np.arange(len(stat_list)))

    xtick_labels = ax.set_xticklabels(region_names, rotation='vertical')
    ytick_labels = ax.set_yticklabels(region_names)

    # cbar = fig.colorbar(im)
    
    if return_ax:
        return ax
    
    else:
        plt.tight_layout()
        plt.savefig(plot_file, **savefig_kwargs)

###
def plot_scatter(stat_list, stat="scf", column="Pdot_shell_mid", column_low="Pdot_shell_low", column_hi="Pdot_shell_hi",
        table="feedback_table.ecsv", 
        region_list=["North", "Central", "South", "L1641N"],
        feedback_label="Shell Momentum Injection Rate",
        plot_file = "scf_slope_pdotshell.pdf", mpl_style='presentation', logy=False):
    """
    Plot a scatter plot of a stat measurement against a feedback measurement.
    """
    matplotlib.style.use(mpl_style)
    print(region_list)
    mid = region_feedback(column, region_list=region_list).data
    if column_low:
        low = region_feedback(column_low, region_list=region_list).data
    else:
        low = mid
    if column_hi:
        hi = region_feedback(column_hi, region_list=region_list).data
    else:
        hi = mid
    print(low, mid, hi)
    fig = plt.figure(figsize=(6,6))
    ax = plt.subplot(111)
    
    if stat == "scf":
        stat_measure, stat_measure_err = slope_scf(stat_list)
        stat_label = "SCF Slope"
    elif stat == "pspec":
        stat_measure, stat_measure_err = slope_pspec(stat_list)
        stat_label = "Power Spectrum Slope"
    else:
        raise(Exception("stat must be one of ['scf', 'pspec']"))
    print(stat_measure, mid, stat_measure_err,[mid-low, hi-mid]) 
    eb = ax.errorbar(stat_measure, mid, xerr=stat_measure_err, yerr=[mid-low, hi-mid], fmt=".")
    if logy:
        plt.semilogy()
    
    yla = ax.set_ylabel("{}".format(feedback_label))
    xla = ax.set_xlabel(stat_label)

    plt.tight_layout()
    plt.savefig(plot_file)
    
    
 ### Functions for plotting stats.
def plot_cov(pca, ax=None, plot_file="cov.pdf", mpl_style='presentation', imshow_kwargs={},
        return_ax=True, axis_label=r"v$_{\rm LSR}$ (km s$^{-1}$)"):
    matplotlib.style.use(mpl_style)

    try:
        pca = pickle.load(open(pca, 'rb'))
    except:
        pass

    cov = pca.cov_matrix

    if ax is None:
        
        fig = plt.figure()
        ax = plt.subplot(111)
        
    extent = 2*[pca.header['CRVAL3']/1000.,
            pca.header['CRVAL3']/1000. + pca.header['CDELT3']*pca.header['NAXIS3']/1000.]
    im = ax.imshow(cov, origin='lower', extent=extent, interpolation='none', **imshow_kwargs)
    ax.set_xlabel(axis_label)
    ax.set_ylabel(axis_label)
    ax.set_xlim([4,16])
    ax.set_ylim([4,16])
    ax.set_aspect(1)
 
    if return_ax:
        return ax
        
    else: 
        plt.tight_layout()
        plt.savefig(plot_file)


def plot_scf_spectrum(scf_list, ax=None, return_ax=False, plot_file="scf_spectrum.pdf",
        colors=[None], labels=[None], mpl_style='presentation',
        scale_lags=1, xlabel="Lag [pixels]", ylabel="SCF",
        plot_points=True, plot_fit=True,
        errorbar_kwargs={}, plot_kwargs={},
        show_fitrange=False, show_beamwidth=False):

    """
    pix_scale scales the lags by a factor to, say,
    convert from lags in pixel to lags in degrees or pc. 
    """
    matplotlib.style.use(mpl_style)

    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)
    
    

    for scf,label,color in zip(list(scf_list),list(labels),list(colors)):
        try:
            scf = pickle.load(open(scf, 'rb'))
        except:
            raise
        
        lags = scf.lags.value*scale_lags

        if plot_points:
            ax.errorbar(lags, scf.scf_spectrum, fmt='.', yerr=scf.scf_spectrum_stddev,
                    color=color, **errorbar_kwargs)
        if plot_fit:
            ax.plot(lags,
                     scf.fitted_model(scf.lags),
                     '--', color=color,
                     label=r"{}: $\alpha = {} \pm {}$".format(
                         label, round(scf.fit.params[1],3), round(scf.fit.bse[1],3)),
                     **plot_kwargs)

        if show_fitrange:
            ax.axvline(scf.xlow.value, ymax=0.2)
            ax.axvline(scf.xhigh.value, ymax=0.2)
            pass
        if show_beamwidth:
            bmaj_pix = scf.header['BMAJ'] / scf.header['CDELT2']
            ax.axvline(bmaj_pix, ymax=0.2, ls='dotted')
            pass

    if ax is None:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.legend()a

    if return_ax:
        return ax
    else:
        plt.tight_layout()
        plt.savefig(plot_file)
    
def plot_pspec(pspec, ax=None, return_ax=True,
        norm_factor=1.,
        plot_file="pspec.pdf",
        mpl_style="presentation",
        ax2_units="pc",
        fill_between=True, twin_axis=True,
        twin_axis_ticks=np.array([0.01,0.1,1,10]),
        twin_axis_labels=True,
        fill_between_kwargs={'alpha':0.3, 'color':'grey'},
        scatter=True,
        scatter_kwargs={'s':10, 'color':'grey'},
        line=True,
        plot_kwargs={'color':'black'},
        connect_points=False,
        connect_plot_kwargs={'color':'grey', 'linestyle':'-'},
        refit=False,
        refit_kwargs={},
        show_fitrange=False, show_beamwidth=False):
    matplotlib.style.use(mpl_style)
    try:
        pspec = pickle.load(open(pspec, 'rb'))
    except:
        raise

    if refit:
        pspec.fit_pspec(**refit_kwargs)
    
    logerr = np.log10(np.e) * pspec.ps1D_stddev / pspec.ps1D
    
    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)

    if fill_between:
        ax.fill_between(np.log10(pspec.freqs.value), np.log10(pspec.ps1D/norm_factor) - logerr,
                        np.log10(pspec.ps1D/norm_factor) + logerr, **fill_between_kwargs)
    if scatter:
        ax.scatter(np.log10(pspec.freqs.value), np.log10(pspec.ps1D/norm_factor), **scatter_kwargs)

    if connect_points:
        ax.plot(np.log10(pspec.freqs.value), np.log10(pspec.ps1D/norm_factor), **connect_plot_kwargs)

    if line:
        ax.plot(np.log10(pspec.freqs.value),
                (pspec.fit.params[0] + pspec.fit.params[1]*np.log10(pspec.freqs.value))/norm_factor,
                label=r"slope: {} $\pm$ {}".format(round(pspec.fit.params[1],2), round(pspec.fit.bse[1],2)),
                **plot_kwargs)

    

    # ax.legend()
    if twin_axis:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())

        ax2_tick_locations = np.array(twin_axis_ticks)
        
        if ax2_units == 'arcmin':
            logk = arcmin_to_logk(ax2_tick_locations)
            ax2.set_xlabel(r"$\theta$ [arcmin]")

        elif ax2_units == 'pc':
            logk = pc_to_logk(ax2_tick_locations)
            
            ax2.set_xlabel(r"$r$ [pc]")

        ax2.set_xticks(logk)
        ax2.set_xticklabels(["{0:g}".format(a) for a in ax2_tick_locations])
        if not twin_axis_labels:
            ax2.set_visible('False')

    if show_fitrange:
        ax.axvline(np.log10(pspec.low_cut.value), ymin=0.9, ymax=1)
        ax.axvline(np.log10(pspec.high_cut.value), ymin=0.9, ymax=1)
        print(np.log10(pspec.high_cut.value))
        print(np.log10(pspec.low_cut.value))

    if show_beamwidth:
        bmaj_pix = pspec.header['BMAJ'] / pspec.header['CDELT2']
        ax.axvline(np.log10(1/bmaj_pix), ymin=0.9, ymax=1, ls='dotted')
        print(np.log10(1/bmaj_pix))

    ax.set_ylabel(r"log P$_2$(k)")
    ax.set_xlabel("log k [1 / pix]")

    if return_ax:
        return ax

    else:
        plt.tight_layout()
        plt.savefig(plot_file)
   

def plot_vcs_spectrum(vcs, ax=None, return_ax=False, plot_file="vcs_spectrum.pdf",
        mpl_style='presentation',
        xlabel=r"log $k_v$ [s/m]", ylabel=r"log $P_1(k_v)$",
        plot_points=True, plot_fit=True, xunit=u.s/u.m,
        scatter_kwargs={}, plot_kwargs={}, axvline_kwargs={},
        show_fitrange=False, refit=False, refit_kwargs={}):

    """
    pix_scale scales the lags by a factor to, say,
    convert from lags in pixel to lags in degrees or pc. 
    """
    matplotlib.style.use(mpl_style)

    if ax is None:
        fig = plt.figure()
        ax = plt.subplot(111)


    try:
        vcs = pickle.load(open(vcs, 'rb'))
    except:
        raise
    
    shape = vcs.freqs.size
    #Freqs and ps1D have redundant values, take first half only.
    rfreqs = vcs.freqs[1:shape // 2]
    ps1D = vcs.ps1D[1:shape // 2]
    freq = vcs._spectral_freq_unit_conversion(rfreqs, xunit)

    if refit:
        vcs.fit_pspec(**refit_kwargs)
        lags = vcs.lags.value*scale_lags

    if plot_points:
        ax.scatter(freq.value, ps1D, marker='.',
                **scatter_kwargs)

    if plot_fit:
        ax.plot(freq.value,
                 10.**vcs.fit.fit.fittedvalues,
                 '-',
                 **plot_kwargs)

    if show_fitrange:
        low_cut = vcs._spectral_freq_unit_conversion(vcs.low_cut, xunit)
        high_cut = vcs._spectral_freq_unit_conversion(vcs.high_cut, xunit)

        ax.axvline(low_cut.value, ymax=0.2, **axvline_kwargs)
        ax.axvline(high_cut.value, ymax=0.2, **axvline_kwargs)

    if not return_ax:
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    # ax.legend()a

    if return_ax:
        return ax
    else:
        plt.tight_layout()
        plt.savefig(plot_file)



def calc_mom0(cube=None, vel_range=None, writeto=None):
    if not isinstance(cube, SpectralCube):
        cube = SpectralCube.read(cube)
    
    slab = cube.spectral_slab(vel_range[0], vel_range[1])
    mom0 = slab.moment0()

    if writeto:
        mom0.write(writeto, overwrite=True)
    return mom0
        

def arcmin_to_logk(arcmin, arcsec_per_pix=2.):
    logk = np.log10(1/(arcmin*60/arcsec_per_pix))
    return logk 
def pc_to_logk(pc, arcsec_per_pix=2., dist=414*u.pc):
    logk = np.log10(1/(206265. * pc / dist.value / arcsec_per_pix))
    return logk

if __name__ == "__main__":
    main()

