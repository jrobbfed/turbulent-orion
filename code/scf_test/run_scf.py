from glob import glob
from spectral_cube import SpectralCube
from turbustat.statistics import SCF
import numpy as np
import astropy.units as u


files = np.array(glob("../subcubes/13co_*.fits"))
files = files[['unger' not in file for file in files]] 

for file in files:
    save_name = "scf_" + file.split('/')[-1].split('.')[0] + "_5to200pix"
    if not glob(save_name+"*"):
        scf = SCF(SpectralCube.read(file), roll_lags=np.arange(-200, 201, 5
    )*u.pix)
        print("Running SCF on {}".format(save_name))
        scf.run(boundary='cut')
        scf.fit_plaw(50*u.pix, 100*u.pix)
        scf.save_results(save_name)


files = np.array(glob("../subcubes/c18o_*.fits"))
files = files[['unger' not in file for file in files]] 

for file in files:
    save_name = "scf_" + file.split('/')[-1].split('.')[0] + "_5to200pix"
    if not glob(save_name+"*"):
        scf = SCF(SpectralCube.read(file), roll_lags=np.arange(-200, 201, 5
    )*u.pix)
        print("Running SCF on {}".format(save_name))
        scf.run(boundary='cut')
        scf.fit_plaw(50*u.pix, 100*u.pix)
        scf.save_results(save_name)