from astropy.io import fits
from sdss5_spectra import SDSSV_Spectra
import numpy as np


file = fits.open('spAll-master.fits')[1].data
high_pm = ((file['GAIA_PMRA'] **2 + file['GAIA_PMDEC'] **2 > 40 **2) &
           (file['GAIA_PMRA'] != -100.) &
           (file['GAIA_PMDEC'] != -100.) &
           (file['SURVEY'] == 'MWM'))
catalogids = list(file['CATALOGID'][high_pm])
# init the class and when doing so update the spAll file
spec = SDSSV_Spectra(catalogid=catalogids,
                     update_spAll=True)
# get the urls for the spectra (not best way to do this)
urls = spec.get_urls()
print(len(urls))
# download first two spectra with rsync
spec = SDSSV_Spectra(catalogid=[catalogids[0],
                                catalogids[1]])
paths = spec.rsync()
#show paths where spectra are now on your computer
print(paths)
