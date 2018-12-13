import turbustat
from turbustat.statistics import PowerSpectrum
from spectral_cube import SpectralCube
import pickle
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from turbustat.tests.generate_test_images import make_extended
from turbustat.io.sim_tools import create_fits_hdu
from astropy.io import fits

def main():
    # pixel_scale = 2*u.arcsec
    # bmin = 8.09*u.arcsec
    # bmaj = 10.01*u.arcsec
    # bpa = -12.9*u.deg
    # restfreq=115.271203*u.GHz
    # bunit = u.K

    # nbmaj_highcut = 3
    # test_pspec(plotname="pspec_rnoise_beamsmooth_beamcorrect.pdf",
    #         run_kwargs={
    #             'verbose':True, 'beam_correct':True, 'apodize_kernel':None,
    #             'xunit':u.pix**-1,
    #             'low_cut':0.01 / u.pix,
    #             'high_cut': 1. /(nbmaj_highcut * bmaj/pixel_scale * u.pix)},
    #         pixel_scale=pixel_scale, bmin=bmin, bmaj=bmaj, bpa=bpa,
    #         restfreq=restfreq, bunit=bunit)

    stat = "pspec_beamcorrect_apodizetukey"
    mols = ["12co", "13co", "c18o"]
    region = "jrf_l1641n"

    # pspec_noise(run_kwargs={'verbose':True})
    # pspec_noise(cube="../subcubes/12co_davis_hh34.fits", vel_low = -2*u.km/u.s, vel_hi=2*u.km/u.s,
    #      run_kwargs={'verbose':True}, pspec_file="pspec_noise_12co_hh34")
    pspec_noise(cube="../subcubes/13co_jrf_l1641n.fits", vel_low = 13*u.km/u.s, vel_hi=16*u.km/u.s,
         run_kwargs={'verbose':True}, pspec_file="pspec_noise_13co_l1641n")
    

    

    # plot_pspec("pspec_noise_12co_hh34.pkl", save_fig="pspec_noise_12co_hh34.png")
    # plot_pspec("pspec_noise_13co_l1641n.pkl", save_fig="pspec_noise_13co_l1641n.png")
    # for mol in mols:
    #     basename = "{}_{}".format(mol, region) 
    #     file = "../subcubes/{}.fits".format(basename)
        
    #     turbustat_version = turbustat.__version__
    #     print(turbustat_version)
    #     pickle_file = "{}_{}_turbustat-{}.p".format(stat, basename, turbustat_version)
    #     plot_file = "{}_{}_turbustat-{}.pdf".format(stat, basename, turbustat_version)
      

    #     # pspec = run_pspec(file, distance=414*u.pc, xunit=u.pix**-1,
    #     #         pickle_file=pickle_file, beam_correct=True, apodize_kernel='tukey')
    #     pl = plot_pspec(pickle_file, save_fig=plot_file,
    #         plot_fit_kwargs={"fit_color":"black", 'show':False}, refit=True,
    #         refit_kwargs={"low_cut":low_cut, 'high_cut':high_cut,
    #             'weighted_fit':False})
    #     plt.gcf().clear()


def pspec_noise(cube="../subcubes/c18o_jrf_l1641n.fits",
        vel_low=12*u.km/u.s, vel_hi=16*u.km/u.s, run_kwargs={},
        pspec_file="pspec_noise_c18o_jrf_l1641n"):
    try:
        cube = SpectralCube.read(cube)
    except:
        pass

    noise_hdu = cube.spectral_slab(vel_low, vel_hi).moment0().hdu
    noise_pspec = PowerSpectrum(noise_hdu)
    noise_pspec.run(**run_kwargs)

    noise_pspec.save_results(pspec_file)
    
    return noise_pspec




def test_pspec(plotname="pspec_rnoise_beamsmooth_apodizetukey.pdf",size=256, powerlaw=3.,
        run_kwargs={'verbose':False, 'apodize_kernel':'tukey'},
        plot_kwargs={'fit_color':'black'},
        beam_smooth=True, pixel_scale = 2*u.arcsec, bmin=8.09*u.arcsec,
        bmaj=10.01*u.arcsec, bpa=-12.9*u.deg, restfreq=1.4*u.GHz, bunit=u.K
        ):
    from spectral_cube import Projection
    from radio_beam import Beam
    
    rnoise_img = make_extended(size, powerlaw)
    # Create a FITS HDU
    rnoise_hdu = create_fits_hdu(rnoise_img, 2*u.arcsec, 2*u.arcsec,
            rnoise_img.shape, 1.4 * u.GHz, u.K)

    pspec = PowerSpectrum(rnoise_hdu)

    if beam_smooth:
        pencil_beam = Beam(0*u.deg)
        rnoise_proj = Projection.from_hdu(rnoise_hdu).with_beam(pencil_beam)
        new_beam = Beam(bmaj, bmin, bpa)
        rnoise_conv = rnoise_proj.convolve_to(new_beam)

        # hdr = fits.Header(header)
        # rnoise_hdu = fits.PrimaryHDU(rnoise_img, header=hdr)
        pspec = PowerSpectrum(rnoise_conv)

    pspec.run(**run_kwargs) 
    pspec.plot_fit(save_name=plotname, **plot_kwargs)

    return pspec

def plot_pspec(pspec_file, save_fig="pspec.pdf",
        plot_fit_kwargs={"fit_color":"black", 'show':False},
        refit=False, refit_kwargs={"low_cut":0.01/u.pix, "high_cut":0.1/u.pix,
            'weighted_fit':False,}):
    try:
        pspec = PowerSpectrum.load_results(pspec_file)
    except:
        pass
    
    if refit:
       pspec.fit_pspec(**refit_kwargs)
       pspec.save_results(pspec_file, keep_data=True)
    pl = pspec.plot_fit(save_name=save_fig, **plot_fit_kwargs)
    plt.gcf().clear()


def run_pspec(cube, distance=414*u.pc, xunit=u.pix**-1, pickle_file=None,
        keep_data=False, **kwargs):
    cube = SpectralCube.read(cube)
    mom0_hdu = cube.moment0().hdu
    pspec = PowerSpectrum(mom0_hdu, distance=distance)
    pspec.run(xunit=xunit, **kwargs)
    if pickle_file:
        pspec.save_results(pickle_file, keep_data=keep_data) 
    return pspec

if __name__ == "__main__":
    main()