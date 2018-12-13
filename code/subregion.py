#This is code to divide spectral cubes into subregions.
import astropy.units as u

### Subregions of Interest
## From the PCA study of OMC1 by Ungerechts et al. 1997
ungerechts_omc1_ra = [83.85105467*u.deg, 83.77605467*u.deg]
ungerechts_omc1_dec = [-5.492928*u.deg, -5.292928*u.deg]

### Expanded omc1 region (2x width) from above in order to fit SCF lags of up to 200 pixels.
omc1_ra = [83.88855467*u.deg, 83.73855467*u.deg] 
omc1_dec = [-5.492928*u.deg, -5.292928*u.deg]
## From Berne et al. 2014 IRAM survey of Orion A
# berne_dec = [-5.8*u.deg, -4.95*u.deg]
# berne_ra = [5.62*15*u.deg, 5.554*15*u.deg]
## From Feddersen et al. 2018, Expanding Shells in Orion A
jrf_north_ra = [84.39490277777777*u.deg, 83.27602222222221*u.deg]
jrf_north_dec = [-5.323802777777778*u.deg, -4.805402777777778*u.deg]
jrf_central_ra = [84.39551250000001*u.deg, 83.27301250000001*u.deg]
jrf_central_dec = [-5.567415138888889*u.deg, -5.322412638888889*u.deg]
jrf_south_ra = [84.40058055555555*u.deg, 83.27853611111111*u.deg]
jrf_south_dec = [-6.284161111111111*u.deg, -5.567411111111111*u.deg]
jrf_l1641n_ra = [84.75317916666667*u.deg, 83.59196250000001*u.deg]
jrf_l1641n_dec = [-7.195631944444444*u.deg, -6.284156944444445*u.deg]
## From Davis et al. 2009, A census of molecular hydrogen outflows
## and their sources along the Orion A molecular ridge.
davis_omc23_ra = [84.01166666666667*u.deg, 83.67833333333333*u.deg]
davis_omc23_dec = [-5.276972222222223*u.deg, -4.943638888888889*u.deg]
davis_ngc1977_ra = [83.94583333333334*u.deg, 83.61666666666666*u.deg]
davis_ngc1977_dec = [-4.952083333333333*u.deg, -4.622916666666666*u.deg]
davis_omc4_ra = [83.92001805555556*u.deg, 83.58305694444444*u.deg]
davis_omc4_dec = [-5.7821855555555555*u.deg, -5.549914444444445*u.deg]
davis_omc5_ra = [84.00006666666668*u.deg, 83.67006666666667*u.deg]
davis_omc5_dec = [-6.251913888888889*u.deg, -5.834913888888889*u.deg]
davis_l1641n_ra = [84.33233750000001*u.deg, 84.0010875*u.deg]
davis_l1641n_dec = [-6.704083333333333*u.deg, -6.317316666666667*u.deg]
davis_hh34_ra = [84.00102500000001*u.deg, 83.66769166666667*u.deg]
davis_hh34_dec = [-6.668663888888888*u.deg, -6.283247222222222*u.deg]
davis_v380_ra = [84.270833333*u.deg, 84.00000*u.deg]
davis_v380_dec = [-6.91666667*u.deg, -6.616666667*u.deg]



regions_ra = [
    ungerechts_omc1_ra, 
    # omc1_ra,
    # berne_ra, 
    jrf_north_ra, 
    jrf_central_ra, 
    jrf_south_ra, 
    jrf_l1641n_ra, 
    davis_omc23_ra, 
    davis_ngc1977_ra, 
    davis_omc4_ra, 
    davis_omc5_ra, 
    davis_l1641n_ra, 
    davis_hh34_ra, 
    davis_v380_ra
    ] 

regions_dec = [
    ungerechts_omc1_dec, 
    # omc1_dec,
    # berne_dec, 
    jrf_north_dec, 
    jrf_central_dec, 
    jrf_south_dec, 
    jrf_l1641n_dec, 
    davis_omc23_dec, 
    davis_ngc1977_dec, 
    davis_omc4_dec, 
    davis_omc5_dec, 
    davis_l1641n_dec, 
    davis_hh34_dec, 
    davis_v380_dec
    ]
regions_name = [
   'ungerechts_omc1', 
   # 'o_omc1',
   # 'berne', 
   'jrf_north', 
   'jrf_central', 
   'jrf_south', 
   'jrf_l1641n', 
   'davis_omc23', 
   'davis_ngc1977', 
   'davis_omc4', 
   'davis_omc5', 
   'davis_l1641n', 
   'davis_hh34', 
   'davis_v380'
        ]
def main():
    print("Executing main function in subregion.py")

    plot_subregions()

    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_omc4_ra,davis_omc4_dec, out_file="13co_davis_omc4.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_omc5_ra,davis_omc5_dec, out_file="13co_davis_omc5.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_l1641n_ra,davis_l1641n_dec, out_file="13co_davis_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_hh34_ra,davis_hh34_dec, out_file="13co_davis_hh34.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_omc23_ra,davis_omc23_dec, out_file="13co_davis_omc23.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_ngc1977_ra,davis_ngc1977_dec, out_file="13co_davis_ngc1977.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits", davis_v380_ra,davis_v380_dec, out_file="13co_davis_v380.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_omc4_ra,davis_omc4_dec, out_file="12co_davis_omc4.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_omc5_ra,davis_omc5_dec, out_file="12co_davis_omc5.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_l1641n_ra,davis_l1641n_dec, out_file="12co_davis_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_hh34_ra,davis_hh34_dec, out_file="12co_davis_hh34.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_omc23_ra,davis_omc23_dec, out_file="12co_davis_omc23.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_ngc1977_ra,davis_ngc1977_dec, out_file="12co_davis_ngc1977.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", davis_v380_ra,davis_v380_dec, out_file="12co_davis_v380.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_omc4_ra,davis_omc4_dec, out_file="c18o_davis_omc4.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_omc5_ra,davis_omc5_dec, out_file="c18o_davis_omc5.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_l1641n_ra,davis_l1641n_dec, out_file="c18o_davis_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_hh34_ra,davis_hh34_dec, out_file="c18o_davis_hh34.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_omc23_ra,davis_omc23_dec, out_file="c18o_davis_omc23.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_ngc1977_ra,davis_ngc1977_dec, out_file="c18o_davis_ngc1977.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits", davis_v380_ra,davis_v380_dec, out_file="c18o_davis_v380.fits")


    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits",
    #         jrf_north_ra,jrf_north_dec, out_file="13co_jrf_north.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits",
    #         jrf_south_ra,jrf_south_dec, out_file="13co_jrf_south.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits",
    #         jrf_l1641n_ra,jrf_l1641n_dec, out_file="13co_jrf_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits",
    #         jrf_central_ra,jrf_central_dec, out_file="13co_jrf_central.fits") 
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits",
    #         jrf_north_ra,jrf_north_dec, out_file="12co_jrf_north.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits",
    #         jrf_south_ra,jrf_south_dec, out_file="12co_jrf_south.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits",
    #         jrf_l1641n_ra,jrf_l1641n_dec, out_file="12co_jrf_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits",
    #         jrf_central_ra,jrf_central_dec, out_file="12co_jrf_central.fits") 
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits",
    #         jrf_north_ra,jrf_north_dec, out_file="c18o_jrf_north.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits",
    #         jrf_south_ra,jrf_south_dec, out_file="c18o_jrf_south.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits",
    #         jrf_l1641n_ra,jrf_l1641n_dec, out_file="c18o_jrf_l1641n.fits")
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits",
    #         jrf_central_ra,jrf_central_dec, out_file="c18o_jrf_central.fits") 

    # subregion("../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits",
    #         ungerechts_omc1_ra,ungerechts_omc1_dec, out_file="12co_ungerechts_omc1.fits") 
    # subregion("../../carma_orion/mask_imfit_13co_pix_2_Tmb.fits",
    #         ungerechts_omc1_ra,ungerechts_omc1_dec, out_file="13co_ungerechts_omc1.fits") 
    # subregion("../../carma_orion/mask_imfit_c18o_pix_2_Tmb.fits",
    #         ungerechts_omc1_ra,ungerechts_omc1_dec, out_file="c18o_ungerechts_omc1.fits") 

def subregion(in_file, ra_range, dec_range, vel_range=None, out_file="subregion.fits", write=True):
    """
    in_file: str
    Input spectral cube file.

    ra_range: tuple or list of angular astropy Quantities 
    dec_range: tuple or list of angular astropy Quantities 
    vel_range: tuple or list of velocity astropy Quantities 

    """
    from spectral_cube import SpectralCube
    cube = SpectralCube.read(in_file)
    print(cube)
    subcube = cube.subcube(ra_range[0], ra_range[1], dec_range[0], dec_range[1])
    if vel_range:
        subcube = subcube.spectral_slab(vel_range)
    if write:
        subcube.write(out_file)
    else:
        return subcube

def plot_subregions(cube="../../carma_orion/mask_imfit_12co_pix_2_Tmb.fits", plotfile="regions.pdf",
        mode='peak', regions_ra=regions_ra, regions_dec=regions_dec, regions_name=regions_name):

    import matplotlib.pyplot as plt
    from astropy.wcs import WCS
    from astropy.io import fits
    from spectral_cube import SpectralCube
    from matplotlib.patches import Rectangle
    
    cube = SpectralCube.read(cube)
    if mode == 'peak':
        data = cube.max(axis=0).data

    wcs = cube.wcs
    ax = plt.subplot(projection=wcs.celestial)
    ax.imshow(data, origin='lower')

    for ra, dec, name in zip(regions_ra, regions_dec, regions_name):
        print((ra[0] - ra[1]).value, (dec[1] - dec[0]).value)
        r = Rectangle((ra[1].value, dec[0].value),
                (ra[0] - ra[1]).value, (dec[1] - dec[0]).value,
                facecolor='none', edgecolor='red', label=name,
                transform=ax.get_transform('world'))
        ax.add_patch(r)
    # ax.legend()
    plt.savefig(plotfile)
    
if __name__ == '__main__':
    main()