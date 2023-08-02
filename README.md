# Download SDSS-V Spectra

This repository provides a python class to download SDSS-V spectra based on information in the spAll file. Currently, the class sdss5_spectra.SDSSV_Spectra only works with BOSS spectra and can download data via rsync.

# Example

There is a short example in the file ``example_get_high_pm_spectra.py``, but in summary a spectrum for some catalogid can be downloaded by:

```
from sdss5_spectra import SDSSV_Spectra

catalogid = 27021597911984196
spec = SDSSV_Spectra(catalogid=[catalogid],
                    update_spAll=True)
paths = spec.rsync()
```
