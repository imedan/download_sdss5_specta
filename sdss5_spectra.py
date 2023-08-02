from astropy.io import fits
from sdss_access.path import Path
import numpy as np
from sdss_access import Access
from os.path import exists
from sdss_access.path.path import AccessError


def update_spALL_file(path, access):
    """
    update the local spAll file
    """
    access.remote()
    access.add('spAll', run2d='master')
    access.set_stream()
    access.commit()


def compare_gaia_to_spAll(file):
    """
    Check which catalogids in gaia
    comparison file have been observed so far

    Parameters
    ---------
    file: str
        path to gaia file with catalogids

    Returns
    -------
    catalogids: list
        List of catalogids currently observed
    """
    path = Path(release='sdss5', preserve_envvars=True)
    spAll_path = path.full('spAll', run2d='master')
    if exists(spAll_path):
        spAll = fits.open(spAll_path)[1].data
    else:
        spAll = fits.open('spAll-master.fits')[1].data
    gaia = np.genfromtxt(file, dtype=int, skip_header=1,
                         usecols=(0))
    obs = np.isin(gaia,
                  list(spAll['CATALOGID'][spAll['SURVEY'] == 'MWM']))
    catalogids = list(gaia[obs])
    return catalogids


class SDSSV_Spectra(object):
    """
    class for downloading SDSSV spectra

    Parameters
    ----------
    catalogid: int or list
        Either one or multiple catalogids (in the form of a list)
        of stars' sepctra you want to download

    instrument: str
        Either 'BOSS' or 'APOGEE'. Depending on instrument,
        will either query spALL (for BOSS) or allStar? (for APOGEE)

    update_spAll: boolean
        Whether or not to update the downloaded spAll file
        before checking for spectra.

    Attributes
    ----------
    file: np.array
        File containing summary data for instrument

    path: sdss_access.path.Path object
        The sdss_access path object

    access: sdss_access.Access object
        The object to rsync data

    ind_where: list
        List of the indexes in the file that have the
        requested catalogids.
    """
    def __init__(self, catalogid, instrument='BOSS',
                 update_spAll=False):
        self.catalogid = catalogid
        self.path = Path(release='sdss5', preserve_envvars=True)
        self.access = Access(release='sdss5')
        if update_spAll:
            update_spALL_file(self.path, self.access)
        # use sas version of spAll if it exists
        if instrument ==  'BOSS':
            spAll_path = self.path.full('spAll', run2d='master')
            if exists(spAll_path):
                self.file = fits.open(spAll_path)[1].data
            else:
                self.file = fits.open('spAll-master.fits')[1].data
        # get the indexes in the filed where catalogids are
        self.ind_where = []
        # no longer str it seems like?
        if isinstance(self.catalogid, int) or isinstance(self.catalogid, np.int64):
            cat_list = [self.catalogid]
        else:
            cat_list = list(self.catalogid)
        self.ind_where = list(np.where(np.isin(self.file['CATALOGID'],
                                               cat_list))[0])

    def get_urls(self):
        """
        get urls for spectra desired

        Returns
        -------
        urls: list
            List of urls for downloading spectra
        """
        urls = []
        for i in self.ind_where:
            url = self.path.url('specFull',
                                run2d=self.file['RUN2D'][i],
                                mjd=self.file['MJD'][i],
                                fieldid='%06d' % self.file['FIELD'][i],
                                isplate='',
                                catalogid=self.file['CATALOGID'][i])
            urls.append(url)
        return urls

    def get_local_paths(self):
        """
        get the location of the spectra on your local
        machine

        Returns
        -------
        paths: dict
            Dictonary of the paths to the spectra where the dict
            keys are the individual catalogids
        """
        paths = {}
        for i in self.ind_where:
            path = self.path.full('specFull',
                                  run2d=self.file['RUN2D'][i],
                                  mjd=self.file['MJD'][i],
                                  fieldid='%06d' % self.file['FIELD'][i],
                                  isplate='',
                                  catalogid=self.file['CATALOGID'][i])
            paths[self.file['CATALOGID'][i]] = path
        return paths

    def rsync(self, return_paths=True,
              check_duplicate_downloads=True,
              chunk_size=None):
        """
        download data via rsync

        Parameters
        ----------
        return_paths: boolean
            Whether or not to return paths of
            downloaded spectra

        check_duplicate_downloads: boolean
            Whether or not to check is spectra has
            already been downloaded

        Returns
        -------
        paths: dict
            Returns paths if return_paths = True
        """
        # make connection
        self.access.remote()
        # add spectra to download
        if check_duplicate_downloads:
            paths = self.get_local_paths()
        if chunk_size is None:
            spec_added = 0
            for i in self.ind_where:
                # check if copy of spectra already exists
                if check_duplicate_downloads:
                    if not exists(paths[self.file['CATALOGID'][i]]):
                        self.access.add('specFull',
                                        run2d=self.file['RUN2D'][i],
                                        mjd=self.file['MJD'][i],
                                        fieldid='%06d' % self.file['FIELD'][i],
                                        isplate='',
                                        catalogid=self.file['CATALOGID'][i])
                        spec_added += 1
                else:
                    self.access.add('specFull',
                                    run2d=self.file['RUN2D'][i],
                                    mjd=self.file['MJD'][i],
                                    fieldid='%06d' % self.file['FIELD'][i],
                                    isplate='',
                                    catalogid=self.file['CATALOGID'][i])
                    spec_added += 1
            # download all spectra
            if spec_added > 0:
                self.access.set_stream()
                self.access.commit()
            else:
                print('All Spectra Already Downloaded')
        else:
            for c in range(0, len(self.ind_where), chunk_size):
                spec_added = 0
                for i in self.ind_where[c: c + chunk_size]:
                    # check if copy of spectra already exists
                    if check_duplicate_downloads:
                        if not exists(paths[self.file['CATALOGID'][i]]):
                            self.access.add('specFull',
                                            run2d=self.file['RUN2D'][i],
                                            mjd=self.file['MJD'][i],
                                            fieldid='%06d' % self.file['FIELD'][i],
                                            isplate='',
                                            catalogid=self.file['CATALOGID'][i])
                            spec_added += 1
                    else:
                        self.access.add('specFull',
                                        run2d=self.file['RUN2D'][i],
                                        mjd=self.file['MJD'][i],
                                        fieldid='%06d' % self.file['FIELD'][i],
                                        isplate='',
                                        catalogid=self.file['CATALOGID'][i])
                        spec_added += 1
                # download all spectra
                if spec_added > 0:
                    self.access.set_stream()
                    self.access.commit()
                else:
                    print('All Spectra Already Downloaded')
                self.access.reset()
        if return_paths:
            return self.get_local_paths()
