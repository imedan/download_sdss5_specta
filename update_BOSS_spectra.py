from astropy.io import fits
from sdss5_spectra import SDSSV_Spectra, update_spALL_file, compare_gaia_to_spAll
import numpy as np
from sdss_access import Access
from sdss_access.path import Path
import os
import shutil
from os.path import exists


# first update the spALL file
path = Path(release='sdss5', preserve_envvars=True)
access = Access(release='sdss5')
update_spALL_file(path, access)

# compare gaia list to whats in spALL
file = 'gaia_v05_match.txt'
catalogids = compare_gaia_to_spAll(file)
print('There are %d stars of interest in the spAll file!' % len(catalogids))

# download the spectra
chunk_size = 1000
for c in range(0, len(catalogids), chunk_size):
    spec = SDSSV_Spectra(catalogid=catalogids[c: c + chunk_size])
    # start rsync
    paths = spec.rsync()

path ="/Users/imedan/sas/sdsswork/bhm/boss/spectro/redux/master/spectra/full"
#we shall store all the file names in this list
filelist = []

for root, dirs, files in os.walk(path):
    for file in files:
        #append the file name to the list
        filelist.append(os.path.join(root,file))

save_dir = '/Users/imedan/Dropbox/sdssV_spectra/'

for file in filelist:
    if not exists(save_dir + file.split('/')[-1]):
        shutil.copy(file, save_dir)
