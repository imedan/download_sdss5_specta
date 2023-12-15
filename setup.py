from setuptools import setup


setup(
    name="sdss5_spectra",
    version="0.1.0",
    author="Ilija Medan",
    author_email="ilija.medan@vanderbilt.edu",
    description="Download SDSS-V Spectra",
    url="https://github.com/imedan/download_sdss5_specta",
    license="BSD-3-Clause license",
    py_modules=['sdss5_spectra.sdss5_spectra'],
    install_requires=['sdss_access','numpy','astropy',]
)
