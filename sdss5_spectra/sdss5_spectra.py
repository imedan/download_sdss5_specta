from astropy.io import fits
from sdss_access.path import Path
import numpy as np
from sdss_access import Access, RsyncAccess
from os.path import exists
from sdss_access.path.path import AccessError


def update_spALL_file(path=None, access=None, run2d='master',
                      release='sdsswork'):
    """
    update the local spAll file
    """
    if path is None:
        path = Path(release=release, preserve_envvars=True)
        access = RsyncAccess(release=release)
    access.remote()
    access.add('spAll', run2d=run2d)
    access.set_stream()
    access.commit()


def update_mwmAllStar_file(path=None, access=None, v_astra='0.5.0',
                           release='ipl3'):
    """
    update the local spAll file
    """
    if path is None:
        path = Path(release=release, preserve_envvars=True)
        access = RsyncAccess(release=release)
    access.remote()
    access.add('mwmAllStar', v_astra=v_astra)
    access.set_stream()
    access.commit()



def compare_gaia_to_spAll(file, run2d='master',
                          release='sdsswork'):
    """
    Check which catalogids in gaia
    comparison file have been observed so far

    Parameters
    ---------
    file: str
        path to gaia file with catalogids

    run2d: str
        The run version of the pipeline to us

    release: str
        The release of sdss to use

    Returns
    -------
    catalogids: list
        List of catalogids currently observed
    """
    path = Path(release=release, preserve_envvars=True)
    spAll_path = path.full('spAll', run2d=run2d)
    if exists(spAll_path):
        spAll = fits.open(spAll_path)[1].data
    else:
        spAll = fits.open('spAll-{run2d}.fits'.format(run2d=run2d))[1].data
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
        of stars' sepctra you want to download for BOSS spectra.
        If downloading APOGEE spectra, this should be sdss_id.

    instrument: str
        Either 'BOSS' or 'APOGEE'. Depending on instrument,
        will either query spALL (for BOSS) or allStar? (for APOGEE)

    update_spAll: boolean
        Whether or not to update the downloaded spAll file
        before checking for spectra.

    run2d: str
        The run version of the pipeline to us

    release: str
        The release of sdss to use

    v_astra: str
        Version of astra to use. Needed to get APOGEE spectra.

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
                 update_spAll=False, run2d='master',
                 release='sdsswork', v_astra='0.5.0'):
        self.catalogid = catalogid
        self.path = Path(release=release, preserve_envvars=True)
        self.access = RsyncAccess(release=release)
        self.instrument = instrument
        self.v_astra = v_astra
        if update_spAll:
            if self.instrument == 'BOSS':
                update_spALL_file(path=self.path, access=self.access,
                                  run2d=run2d, release=release)
            else:
                update_mwmAllStar_file(path=self.path, access=self.access,
                                       v_astra=v_astra, release=release)
        # use sas version of spAll if it exists
        if self.instrument ==  'BOSS':
            spAll_path = self.path.full('spAll', run2d=run2d)
            if exists(spAll_path):
                self.file = fits.open(spAll_path)[1].data
            else:
                self.file = fits.open('spAll-{run2d}.fits'.format(run2d=run2d))[1].data
        else:
            mwmAllStar_path = self.path.full('mwmAllStar', v_astra=v_astra)
            if exists(mwmAllStar_path):
                self.file = fits.open(mwmAllStar_path)[2].data
            else:
                self.file = fits.open('mwmAllStar-{v_astra}.fits'.format(v_astra=v_astra))[2].data
        # get the indexes in the filed where catalogids are
        self.ind_where = []
        # no longer str it seems like?
        if isinstance(self.catalogid, int) or isinstance(self.catalogid, np.int64):
            cat_list = [self.catalogid]
        else:
            cat_list = list(self.catalogid)
        if self.instrument ==  'BOSS':
            self.file_id_col = 'CATALOGID'
        else:
            self.file_id_col = 'sdss_id'
        self.ind_where = list(np.where(np.isin(self.file[self.file_id_col],
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
            if self.instrument ==  'BOSS':
                url = self.path.url('specFull',
                                    run2d=self.file['RUN2D'][i],
                                    mjd=self.file['MJD'][i],
                                    fieldid='%06d' % self.file['FIELD'][i],
                                    isplate='',
                                    catalogid=self.file['CATALOGID'][i])
            else:
                url = self.path.url('mwmStar',
                                     v_astra=self.v_astra,
                                     sdss_id=self.file['sdss_id'][i],
                                     component='')
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
            if self.instrument ==  'BOSS':
                path = self.path.full('specFull',
                                      run2d=self.file['RUN2D'][i],
                                      mjd=self.file['MJD'][i],
                                      fieldid='%06d' % self.file['FIELD'][i],
                                      isplate='',
                                      catalogid=self.file['CATALOGID'][i])
            else:
                path = self.path.full('mwmStar',
                                      v_astra=self.v_astra,
                                      sdss_id=self.file['sdss_id'][i],
                                      component='')
            paths[self.file[self.file_id_col][i]] = path
        return paths

    def add_spectra_access(self, i):
        """
        add the spectrum of interest to the rsync
        """
        if self.instrument ==  'BOSS':
            self.access.add('specFull',
                            run2d=self.file['RUN2D'][i],
                            mjd=self.file['MJD'][i],
                            fieldid='%06d' % self.file['FIELD'][i],
                            isplate='',
                            catalogid=self.file['CATALOGID'][i])
        else:
            self.access.add('mwmStar',
                            v_astra=self.v_astra,
                            sdss_id=self.file['sdss_id'][i],
                            component='')

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
                    if not exists(paths[self.file[self.file_id_col][i]]):
                        self.add_spectra_access(i)
                        spec_added += 1
                else:
                    self.add_spectra_access(i)
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
                        if not exists(paths[self.file[self.file_id_col][i]]):
                            self.add_spectra_access(i)
                            spec_added += 1
                    else:
                        self.add_spectra_access(i)
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
