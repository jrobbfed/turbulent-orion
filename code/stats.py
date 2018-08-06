import numpy as np
import astropy.units as u
import pickle
from spectral_cube import SpectralCube
from glob import glob
import os
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
distance = 414*u.pc

def main():

### Plotting scatter plots of SCF and PSPEC slopes versus 
### surface density of YSOs (from Spitzer Orion) in Orion A regions.
# stat_list = glob("scf_results/scf_13co*davis*.p") +\
#         glob("scf_results/scf_13co*unger*.p")
# plot_stat_feedback(stat_list, plot_file="scf_13co_davis_nstars.pdf")
##Figure in paper.


    # mols = ["12co", "13co", "c18o"]

    # mols_dict = {
    #     "12co":r"$^{12}$CO",
    #     "13co":r"$^{13}$CO",
    #     "c18o":r"C$^{18}$O"}

    # stat_names = ["scf", "pspec"]


    # ytick_arr = np.array(
    #         [[[-0.10,-0.15],
    #           [-0.10,-0.20],
    #           [-0.20,-0.30]],
    #          [[-3.0,-4.0],
    #           [-3.0,-3.5],
    #           [-3.2,-3.6]]])
             
    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))
            
    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'
            
    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=True, feedback_label=r"n$_{\rm YSO}$ (deg$^{-2}$)",
    #                 logx=True, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30})
            
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
    # plt.savefig("slope_nstars_logdensity.pdf")
    # plt.close()

    
### Plotting scatter plots of SCF and PSPEC slopes versus 
### number of YSOs (from Spitzer Orion) in Orion A regions.

    #     f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=False, feedback_label=r"N$_{\rm YSO}$",
    #                 logx=True)
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
    # plt.savefig("slope_nstars.pdf")
    # plt.close()    

### Plotting scatter plots of SCF and PSPEC slopes versus 
### region area (from Spitzer Orion) in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True,
    #                 surface_density=False, feedback_label=r"surface area (deg$^2$)",
    #                 logx=False, feedback_mode="area")
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
    # plt.savefig("slope_surface_area.pdf")
    # plt.close()  

### Plotting scatter plots of SCF and PSPEC slopes versus 
### shell pdot surface density in Orion A regions.

    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9.5,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.5)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="outflows_nro45m.csv",
    #                 return_cols="Pdot_flow", ra_col="RA_J2000", dec_col="DEC_J2000",
    #                 surface_density=True, feedback_label=r"$\dot P_{\rm out}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$ deg$^{-2}$)",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum")

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
    # plt.savefig("slope_outflow_pdot_surfacedensity.pdf")
    # plt.close()


    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9.5,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.5)
    
    
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="outflows_nro45m.csv",
    #                 return_cols="Pdot_flow", ra_col="RA_J2000", dec_col="DEC_J2000",
    #                 surface_density=False, feedback_label=r"$\dot P_{\rm out}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$)",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum")

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
    # plt.savefig("slope_outflow_pdot.pdf")
    # plt.close()

# ### Plotting scatter plots of SCF and PSPEC slopes versus 
# ### shell pdot in Orion A regions.

    #f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #        figsize=(9,9.5))
    #f.subplots_adjust(hspace=0, wspace=0.4)
     
    ##Iterate over each row.
    #for irow, axrow, mol in zip(
    #        range(len(mols)), axarr, mols):
    #    #iterate over each column.
    #    for icol, ax, stat_name in zip(
    #            range(len(stat_names)), axrow, stat_names):
     
    #        stat_list = glob("{}_results/{}_{}*.p".format(
    #            stat_name, stat_name, mol))
     
    #        ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #        colors = np.array(['0.7'] * len(stat_list))
    #        colors[ii_davis] = '0.1'
     
    #        ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="shells_pdot_NtoS.csv",
    #                return_cols="Pdot_mid", ra_col="RA", dec_col="DEC",
    #                surface_density=False, feedback_label=r"$\dot P_{\rm shells}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$)l",
    #                logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                scatter_kwargs={'color':colors, 's':30},
    #                feedback_mode="sum")
     
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
    ## ax1.set_xlabel('')
    ## plt.tight_layout()
    #plt.savefig("slope_shell_pdot.pdf")
    #plt.close()

## Plot scatter plots of SCF and PSPEC slopes versus
## outflow pdot surface density (from outflows_nro45m.csv)
    # f, axarr = plt.subplots(nrows=3, ncols=2, sharex=True, sharey=False,
    #         figsize=(9,9.5))
    # f.subplots_adjust(hspace=0, wspace=0.4)
   
    # #Iterate over each row.
    # for irow, axrow, mol in zip(
    #         range(len(mols)), axarr, mols):
    #     #iterate over each column.
    #     for icol, ax, stat_name in zip(
    #             range(len(stat_names)), axrow, stat_names):

    #         stat_list = glob("{}_results/{}_{}*.p".format(
    #             stat_name, stat_name, mol))

    #         ii_davis = ["davis" in stat_file for stat_file in stat_list]
    #         colors = np.array(['0.7'] * len(stat_list))
    #         colors[ii_davis] = '0.1'

    #         ax = plot_stat_feedback(stat_list, ax=ax, return_ax=True, catalog="shells_pdot_NtoS.csv",
    #                 return_cols="Pdot_mid", ra_col="RA", dec_col="DEC",
    #                 surface_density=False, feedback_label=r"$\dot P_{\rm shells}$ (M$_\odot$ km s$^{-1}$ yr$^{-1}$)l",
    #                 logx=False, errorbar_kwargs={'ecolor':colors, 'linewidth':2},
    #                 scatter_kwargs={'color':colors, 's':30},
    #                 feedback_mode="sum")

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
    # plt.savefig("slope_outflow_pdot.pdf")
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
    # pickle_files = ["pspec_results/pspec_"+name+".p" for name in names]
    
    # for i in range(len(cubes)):
    #     if not glob(pickle_files[i]+"*"):
    #         print("Running Power Spectrum on {}. Pickling to {}.".format(cubes[i], pickle_files[i]))
    #         run_pspec(cubes[i], distance=414*u.pc,
    #                 xunit=u.pix**-1,
    #                 pickle_file=pickle_files[i])

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
        table="spitzer_orion.fit", subcube_dir="subcubes/"):
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

def plot_stat_feedback(stat_list, ax=None, catalog="spitzer_orion.fit",
        return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
        ra_col="_RAJ2000", dec_col="_DEJ2000",
        region_names=["davis_omc23", "davis_ngc1977"],
        regions_from_stat_list=True,
        subcube_dir = "subcubes/", feedback_mode='count',
        errorbar_kwargs={}, scatter_kwargs={},
        surface_density=True, feedback_label=r"n$_{\rm YSO}$ (deg$^{-2}$)",
        logx=False, logy=False, mpl_style="presentation",
        return_ax=True,return_eb=False, plot_file="plot_stat_feedback.pdf"):
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
    import subregion 
    matplotlib.style.use(mpl_style)

    if ax is None:
        
        fig = plt.figure(figsize=(6,6))
        ax = plt.subplot(111)
    
    if "scf" in stat_list[0]:
        stat_measure, stat_measure_err = slope_scf(stat_list)
        stat_label = "SCF Slope"
    elif "pspec" in stat_list[0]:
        stat_measure, stat_measure_err = slope_pspec(stat_list)
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
    for region_name in region_names:

        feedback_table = region_stars(region_name=region_name, table=catalog,
                return_cols=return_cols, ra_col=ra_col,
                dec_col=dec_col, subcube_dir=subcube_dir)

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
            area = region_area(unit="deg", dist=414*u.pc,
                region_name=region_name, subcube_dir="subcubes/")

            feedback_total = feedback_total / area.value

        feedback_list.append(feedback_total)

    
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


def plot_colorcolor_stars(star_catalog="spitzer_orion.fit",
        star_return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
        star_ra_col="_RAJ2000", star_dec_col="_DEJ2000",
        region_name="davis_omc23", subcube_dir = "subcubes/"):
    """
    Make color-color plots of stat distance metrics between subregions,
    ranking the subregions by stellar density from a catalog.
    """
    
    
    
    star_table = region_stars(region_name=region_name, table=star_catalog,
            return_cols=star_return_cols, ra_col=star_ra_col,
            dec_col=star_dec_col, subcube_dir=subcube_dir)



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

def run_pspec(cube, distance=414*u.pc, xunit=u.pix**-1, pickle_file=None, **kwargs):
    from turbustat.statistics import PowerSpectrum
    cube = SpectralCube.read(cube)
    mom0_hdu = cube.moment0().hdu
    pspec = PowerSpectrum(mom0_hdu, distance=distance)
    pspec.run(xunit=xunit, **kwargs)
    if pickle_file:
        pickle.dump(pspec, open(pickle_file, "wb"))
    return pspec


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

def slope_pspec(pspec_list):
    try:
        pspec_list = [pickle.load(open(f, 'rb')) for f in pspec_list]
    except:
        pass

    slope = np.array([pspec.slope for pspec in pspec_list])
    slope_err = np.array([pspec.slope_err for pspec in pspec_list])

    return slope, slope_err

### Functions plotting color-color plots of statistical distance metrics
### between subregions of data.
def plot_colorcolor(stat_list,
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
    fig, ax = plt.subplots(1, figsize=(8,8))
    
    im = ax.imshow(dist_array, origin='lower', **imshow_kwargs)
    # xla = ax.set_xlabel("subregion") 
    # yla = ax.set_ylabel("subregion")

    xticks = ax.set_xticks(np.arange(len(stat_list)))
    yticks = ax.set_yticks(np.arange(len(stat_list)))

    xtick_labels = ax.set_xticklabels(region_names, rotation='vertical')
    ytick_labels = ax.set_yticklabels(region_names)

    cbar = fig.colorbar(im)
    
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
    
    
    
def plot_stat_feedback(stat_list, ax=None, catalog="spitzer_orion.fit",
        return_cols=["_RAJ2000", "_DEJ2000", "Cl"],
        ra_col="_RAJ2000", dec_col="_DEJ2000",
        region_names=["davis_omc23", "davis_ngc1977"],
        regions_from_stat_list=True,
        subcube_dir = "subcubes/", feedback_mode='count',
        errorbar_kwargs={}, scatter_kwargs={},
        surface_density=True, feedback_label=r"n$_{\rm YSO}$ (deg$^{-2}$)",
        logx=False, logy=False, mpl_style="presentation",
        return_ax=True,return_eb=False, plot_file="plot_stat_feedback.pdf"):


### Functions for plotting stats.
def plot_cov(pca, ax=None, plot_file="cov.pdf", mpl_style='presentation', imshow_kwargs={},
        return_ax=True, axis_label=r"v$_{\rm LSR}$ (km s$^{-1}$)"):
    matplotlib.style.use(mpl_style)
    if type(pca) == str:
        pca = pickle.load(open(pca, 'rb'))

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


def plot_scf_spectrum(scf_list, plot_file="scf_spectrum.pdf",
        colors=[None], labels=[None], mpl_style='presentation',
        scale_lags=1, xlabel="Lag [pixels]", ylabel="SCF",
        plot_points=True, plot_fit=True,
        **kwargs):

    """
    pix_scale scales the lags by a factor to, say,
    convert from lags in pixel to lags in degrees or pc. 
    """
    matplotlib.style.use(mpl_style)
    fig = plt.figure()
    ax = plt.subplot(111)
    
    for scf,label,color in zip(scf_list,labels,colors):
        if type(scf) == str:
            scf = pickle.load(open(scf, 'rb'))
        
        lags = scf.lags.value*scale_lags

        if plot_points:
            plt.errorbar(lags, scf.scf_spectrum, fmt='.', yerr=scf.scf_spectrum_stddev,
                    color=color)
        if plot_fit:
            plt.plot([lags[0], lags[-1]],
                     [scf.fitted_model(scf.lags[0]), scf.fitted_model(scf.lags[-1])],
                     '--', color=color,
                     label=r"{}: $\alpha = {} \pm {}$".format(
                         label, round(scf.fit.params[1],3), round(scf.fit.bse[1],3)))
    plt.semilogx()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_file)
    
def plot_pspec(pspec, plot_file="pspec.pdf", mpl_style="presentation", ax2_units="pc"):
    matplotlib.style.use(mpl_style)
    if type(pspec) == str:
        pspec = pickle.load(open(pspec, 'rb'))
    
    logerr = np.log10(np.e) * pspec.ps1D_stddev / pspec.ps1D
    


    fig = plt.figure()
    ax = plt.subplot(111)
    ax.fill_between(np.log10(pspec.freqs.value), np.log10(pspec.ps1D) - logerr,
                    np.log10(pspec.ps1D) + logerr, alpha=0.3, color='grey')
    
    ax.scatter(np.log10(pspec.freqs.value), np.log10(pspec.ps1D),10, color='grey')
    print(pspec.freqs.value, pspec.fit.fittedvalues)
    ax.plot([np.log10(pspec.freqs.value)[0],np.log10(pspec.freqs.value)[-1]],
            [pspec.fit.fittedvalues[0], pspec.fit.fittedvalues[-1]], 'black',
         label=r"slope: {} $\pm$ {}".format(round(pspec.fit.params[1],2), round(pspec.fit.bse[1],2)))
    ax.legend()
   
    ax2 = plt.twiny()
    ax2.set_xlim(ax.get_xlim())

    ax2_tick_locations = np.array([0.01, 0.1, 1, 10])
    
    if ax2_units == 'arcmin':
        logk = arcmin_to_logk(ax2_tick_locations)
        ax2.set_xlabel(r"$\theta$ [arcmin]")

    elif ax2_units == 'pc':
        logk = pc_to_logk(ax2_tick_locations)
        
        ax2.set_xlabel(r"$r$ [pc]")

    ax2.set_xticks(logk)
    ax2.set_xticklabels(["{0:g}".format(a) for a in ax2_tick_locations])
    
    
    ax.set_ylabel(r"log P$_2$(k)")
    ax.set_xlabel("log k [1 / pix]")

    plt.tight_layout()
    plt.savefig(plot_file)
    
def arcmin_to_logk(arcmin, arcsec_per_pix=2.):
    logk = np.log10(1/(arcmin*60/arcsec_per_pix))
    return logk 
def pc_to_logk(pc, arcsec_per_pix=2., dist=414*u.pc):
    logk = np.log10(1/(206265. * pc / dist.value / arcsec_per_pix))
    return logk

if __name__ == "__main__":
    main()