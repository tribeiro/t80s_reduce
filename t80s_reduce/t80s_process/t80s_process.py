import os
import numpy as np
import logging
import yaml
import copy
import time
from astropy.io import fits, ascii
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
from padova.isocdata import IsochroneSet

import datetime
from distutils.dir_util import mkpath
from shutil import copyfile, move
# import imreg_dft as ird
import image_registration

from t80s_reduce.core.constants import *
from t80s_reduce.util.imcombine import imcombine
from t80s_reduce.util.imarith import imarith
from t80s_reduce.util.sextractor import SExtractor, SExtractorException
from t80s_reduce.util.scamp import Scamp
from t80s_reduce.util.swarp import Swarp
from t80s_reduce.util.astrometrynet import AstrometryNet, AstrometryNetException
from t80s_reduce.t80s_process.t80s_preprocess import T80SPreProc

log = logging.getLogger(__name__)


class T80SProcess:
    def __init__(self, config=None):

        log.debug('Reading configuration file: %s' % config)
        if config is not None:
            with open(config, 'r') as fp:
                self.config = yaml.load(fp)
        else:
            # Todo: Read sample configuration file
            raise NotImplementedError("Read from sample configuration file not implemented yet.")
        self._config_file = config

    def __del__(self):
        log.debug('Rewriting configuration file %s' % self._config_file)
        filename,ext = os.path.splitext(self._config_file)
        backup_name = "%s-%s%s" % (filename,
                                   time.strftime("%Y%m%d-%H%M%S"),
                                   ext)
        copyfile(self._config_file,
                 backup_name)
        with open(self._config_file, 'w') as fp:
            try:
                yaml.dump(self.config, fp, default_flow_style=False)
            except Exception,e:
                copyfile(backup_name,
                         self._config_file)
                # print self.config
                log.exception(e)
                log.warning('Could not save configuration file. Restoring last backup.')



    @property
    def masterbias(self):
        if 'master' in self.config['calibrations']['bias']:
            master = os.path.join(self.config['path'],
                                  self.config['calibrations']['bias']['night'],
                                  'bias',
                                  self.config['calibrations']['bias']['master'])
            if os.path.exists(master):
                return master
            else:
                raise IOError('Master bias file %s does not exists in defined path. Check your configuration file!'
                              % master)
        else:
            log.debug('Master bias not defined in configuration structure.')
            return None

    @masterbias.setter
    def masterbias(self, value):
        self.config['calibrations']['bias']['master'] = value
        # if not os.path.exists(self.masterbias):
        #     raise IOError('Master bias file does not exists! Create the file before inserting it in the database.')

    @property
    def biaspath(self):
        return os.path.join(self.config['path'],
                            self.config['calibrations']['bias']['night'],
                            'bias')

    def get_bias_list(self):
        return [os.path.join(self.config['path'],
                             self.config['calibrations']['bias']['night'],
                             'bias',
                             'raw_files',
                             raw) for raw in self.config['calibrations']['bias']['raw files']]

    def get_flat_list(self, get_file_type='raw', write_file_type=None, overwrite=False, getfilter=None):
        config = self.config
        img_list = []
        # get list of flats
        filterlist = config['calibrations']['sky-flat']['filters'] if getfilter is None else getfilter

        for filter in filterlist:
            path = os.path.join(config['path'],
                                config['calibrations']['sky-flat'][filter]['night'],
                                'flat',
                                filter)
            if (write_file_type is not None) and overwrite and \
                    (write_file_type in config['calibrations']['sky-flat'][filter]):
                log.warning('Running in overwrite mode. Cleaning existing %s list.' % write_file_type)
                config['calibrations']['sky-flat'][filter][write_file_type] = []
            elif write_file_type is not None and write_file_type in config['calibrations']['sky-flat'][filter]:
                log.warning('Flat field frames on %s already processed but running in non-overwrite mode. '
                            'Skipping.' % filter)
                continue
            elif write_file_type is not None:
                config['calibrations']['sky-flat'][filter][write_file_type] = []

            for flat in config['calibrations']['sky-flat'][filter][get_file_type]:
                if write_file_type is not None:
                    img_list.append((os.path.join(path, get_file_type.replace(' ', '_'), flat),
                                     os.path.join(path, write_file_type.replace(' ', '_'), flat),
                                     ('flat', filter)))
                else:
                    img_list.append((os.path.join(path, get_file_type.replace(' ', '_'), flat),
                                     None,
                                     ('flat', filter)))

        return img_list

    def get_target_list(self, get_file_type='raw files', write_file_type=None, overwrite=False, getobject=None,
                        getfilter=None):
        config = self.config
        img_list = []

        filter_list = getfilter if getfilter is not None else FILTERS
        object_list = getobject if getobject is not None else config['objects']

        for object in object_list:
            if object not in config['objects']:
                log.warning('Object %s not in configuration file! Skipping...' % object)
                continue

            for filter in filter_list:
                if filter not in config['objects'][object]:
                    continue

                rpath = os.path.join(config['path'],
                                     config['objects'][object]['night'],
                                     object.replace(' ', '_'),
                                     filter)
                wpath = rpath if 'wpath' not in self.config else os.path.join(config['wpath'],
                                                                              config['objects'][object]['night'],
                                                                              object.replace(' ', '_'),
                                                                              filter)

                if (write_file_type is not None) and (overwrite) and \
                        (write_file_type in config['objects'][object][filter]):
                    log.warning('Running in overwrite mode. Cleaning existing %s list.' % write_file_type)
                    config['objects'][object][filter][write_file_type] = []
                elif (write_file_type is not None) and (write_file_type in config['objects'][object][filter]):
                    log.warning('%s frames in %s already processed but running in non-overwrite mode. '
                                'Checking files to process.' % (object, filter))
                    for raw in config['objects'][object][filter][get_file_type]:
                        get = os.path.join(rpath, get_file_type.replace(' ', '_'), raw)
                        write = os.path.join(wpath, write_file_type.replace(' ', '_'), raw)
                        if not os.path.exists(write):
                            log.debug('Image %s not processed. Adding to processing list.' % raw)
                            img_list.append((get,
                                             write,
                                             ('target', object, filter)))
                        else:
                            log.debug('Image %s processed, skipping.' % raw)
                            # 'Skipping.' % (object, filter))
                    continue
                else:
                    config['objects'][object][filter][write_file_type] = []

                for raw in config['objects'][object][filter][get_file_type]:
                    rfile = os.path.join(rpath, get_file_type.replace(' ', '_'), raw)
                    wfile = os.path.join(wpath, write_file_type.replace(' ', '_'), raw) if write_file_type is not None \
                        else ''
                    img_list.append((rfile,
                                     wfile,
                                     ('target', object, filter)))
        return img_list

    def imcombine(self, key, overwrite=False):
        if key == 'master-bias':
            biasname = 'masterbias.fits'
            if (self.masterbias is not None) and (not overwrite):
                raise IOError('Master bias already in database. Run in overwrite mode do proceed.')
            mbias = self.masterbias if self.masterbias is not None else os.path.join(self.biaspath, biasname)
            try:
                imcombine(self.get_bias_list(), output=mbias, overwrite=overwrite)
            except Exception, e:
                log.exception(e)
                raise
            else:
                self.masterbias = biasname
                self.config['calibrations']['bias']['master'] = biasname
        elif key == 'master-flat':
            for filter in self.config['calibrations']['sky-flat']['filters']:
                if 'master' in self.config['calibrations']['sky-flat'][filter] and not overwrite:
                    log.warning('Master flat already created. Run in overwrite mode to continue.')
                    continue

                log.debug('Creating master flat in %s' % filter)
                imglist = self.get_flat_list(get_file_type='norm', getfilter=[filter])

                if 'flat-par' in self.config:

                    log.debug('Selecting frames with average count between %s and %s' % (self.config['flat-par']['min'],
                                                                                         self.config['flat-par']['max']
                                                                                         ))
                    good_flats = []
                    for img in imglist:
                        hdr = fits.getheader(img[0])
                        dmean = float(hdr['NORMFAC'])
                        log.debug('%s: %10.2f' % (os.path.basename(img[0]),
                                                  dmean))
                        if ( self.config['flat-par']['min'] < dmean < self.config['flat-par']['max'] ):
                            good_flats.append(img)
                    if len(good_flats) < 3:
                        log.critical('No good flats found in this sequence!')
                        continue
                    elif len(good_flats) < NFLATS:
                        log.warning('Number of good flats is smaller than expected! Found %i of %i. Minumum is %i' % (
                            len(good_flats),
                            len(imglist),
                            NFLATS))
                    imglist = good_flats


                path = os.path.join(self.config['path'],
                                    self.config['calibrations']['sky-flat'][filter]['night'],
                                    'flat',
                                    filter)
                mflat = os.path.join(path,
                                     'masterflat.fits')
                try:
                    flatlist = [entry[0] for entry in imglist]
                    imcombine(flatlist, output=mflat, overwrite=overwrite)
                except Exception, e:
                    log.exception(e)
                    raise
                else:
                    self.config['calibrations']['sky-flat'][filter]['master'] = 'masterflat.fits'

    def biascorr(self, overwrite=False):

        # Building file list
        img_list = self.get_flat_list(get_file_type='raw files', write_file_type='biascorr', overwrite=overwrite)
        for item in self.get_target_list(get_file_type='raw files', write_file_type='biascorr', overwrite=overwrite):
            img_list.append(item)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                imarith(str(img[0]), '-', str(self.masterbias), str(img[1]), overwrite=overwrite)
            except Exception, e:
                log.exception(e)
                continue
            else:
                if img[2][0] == 'flat':
                    self.config['calibrations']['sky-flat'][img[2][1]]['biascorr'].append(os.path.basename(img[1]))
                elif img[2][0] == 'target':
                    self.config['objects'][img[2][1]][img[2][2]]['biascorr'].append(os.path.basename(img[1]))

    def overcorr(self, overwrite=False):
        # Building file list
        img_list = self.get_flat_list(get_file_type='biascorr', write_file_type='overcorr', overwrite=overwrite)
        for item in self.get_target_list(get_file_type='biascorr', write_file_type='overcorr', overwrite=overwrite):
            img_list.append(item)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                overcorr = T80SPreProc()
                log.info('Reading in %s' % img[0])
                overcorr.read('%s' % img[0])

                log.info('Loading configuration from %s' % self.config['overscan_config'])
                overcorr.loadConfiguration(self.config['overscan_config'])

                overcorr.overscan()
                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                elif os.path.exists(img[1]) and overwrite:
                    log.debug('Removing existing file %s' % img[1])
                    os.remove(img[1])
                log.debug('Writing %s' % img[1])
                overcorr.ccd.write(img[1])
            except Exception, e:
                log.exception(e)
                continue
            else:
                if img[2][0] == 'flat':
                    self.config['calibrations']['sky-flat'][img[2][1]]['overcorr'].append(os.path.basename(img[1]))
                elif img[2][0] == 'target':
                    self.config['objects'][img[2][1]][img[2][2]]['overcorr'].append(os.path.basename(img[1]))

    def linearize(self, overwrite=False):

        if 'linearity_coefficients' not in self.config:
            raise IOError('Linearity coefficients file not in configuration file.')

        if 'overscan_config' not in self.config:
            raise IOError('Overscan configuration not in configuration file.')

        log.debug('Reading linearity coefficients from %s.' % self.config['linearity_coefficients'])
        linearity_coefficients = np.load(os.path.expanduser(self.config['linearity_coefficients']))

        # Building file list
        # Todo: Check if overcorr is in the list and use it insted of biascorr
        img_list = self.get_flat_list(get_file_type='overcorr', write_file_type='lincorr', overwrite=overwrite)
        for item in self.get_target_list(get_file_type='overcorr', write_file_type='lincorr', overwrite=overwrite):
            img_list.append(item)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                lincorr = T80SPreProc()
                log.info('Reading in %s' % img[0])
                lincorr.read('%s' % img[0])

                log.info('Loading configuration from %s' % self.config['overscan_config'])
                lincorr.loadConfiguration(self.config['overscan_config'])

                lincorr.linearize(linearity_coefficients)
                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                elif os.path.exists(img[1]) and overwrite:
                    log.debug('Removing existing file %s' % img[1])
                    os.remove(img[1])

                log.debug('Writing %s' % img[1])
                lincorr.ccd.write(img[1])

                # completed[i] = imarith(img[0], '-', masterbias, img[1], overwrite=args.overwrite) == 0
            except Exception, e:
                log.exception(e)
                continue
            else:
                if img[2][0] == 'flat':
                    self.config['calibrations']['sky-flat'][img[2][1]]['lincorr'].append(os.path.basename(img[1]))
                elif img[2][0] == 'target':
                    self.config['objects'][img[2][1]][img[2][2]]['lincorr'].append(os.path.basename(img[1]))

    def trim(self, overwrite=False):

        if 'overscan_config' not in self.config:
            raise IOError('Overscan configuration not in configuration file.')

        # Building file list
        # Todo: Trim different file list
        img_list = self.get_flat_list(get_file_type='lincorr', write_file_type='trim', overwrite=overwrite)
        for item in self.get_target_list(get_file_type='lincorr', write_file_type='trim', overwrite=overwrite):
            img_list.append(item)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                trimcorr = T80SPreProc()
                log.info('Reading in %s' % img[0])
                trimcorr.read('%s' % img[0])

                log.info('Loading configuration from %s' % self.config['overscan_config'])
                trimcorr.loadConfiguration(self.config['overscan_config'])

                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                elif os.path.exists(img[1]) and overwrite:
                    log.debug('Removing existing file %s' % img[1])
                    os.remove(img[1])

                log.debug('Writing %s' % img[1])
                ccdout = trimcorr.trim()
                ccdout.write(img[1])

                # completed[i] = imarith(img[0], '-', masterbias, img[1], overwrite=args.overwrite) == 0
            except Exception, e:
                log.exception(e)
                continue
            else:
                if img[2][0] == 'flat':
                    self.config['calibrations']['sky-flat'][img[2][1]]['trim'].append(os.path.basename(img[1]))
                elif img[2][0] == 'target':
                    self.config['objects'][img[2][1]][img[2][2]]['trim'].append(os.path.basename(img[1]))

    def normalize_flatfield(self, overwrite=False, nhdu=0):
        img_list = self.get_flat_list(get_file_type='trim', write_file_type='norm', overwrite=overwrite)

        # Todo: Paralellize this!
        for img in img_list:
            try:
                log.info('Reading in %s' % img[0])
                hdulist = fits.open(img[0])
                normfactor = np.mean(hdulist[nhdu].data)
                log.debug('Normalizing')
                newdata = np.zeros_like(hdulist[nhdu].data, dtype=np.float32)
                newdata += hdulist[nhdu].data
                newdata /= normfactor
                output_hdu = fits.PrimaryHDU(header=hdulist[nhdu].header,
                                             data=newdata)
                header_comments = ["NORMALIZE: %s" % datetime.datetime.now(),
                                   "NORMALIZE: factor = %f" % normfactor]
                output_hdu.header['NORMFAC'] = normfactor

                for comment in header_comments:
                    output_hdu.header["COMMENT"] = comment

                output_hdulist = fits.HDUList([output_hdu])

                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))

                log.info('Saving output to %s' % img[1])
                output_hdulist.writeto(img[1], clobber=overwrite)

            except Exception, e:
                log.exception(e)
                continue
            else:
                self.config['calibrations']['sky-flat'][img[2][1]]['norm'].append(os.path.basename(img[1]))

    def normalize_flatfield_channel(self, overwrite=False, nhdu=0):
        if 'trimmed_overscan_config' not in self.config:
            raise IOError('Overscan configuration for trimmed data not in configuration file.')

        img_list = self.get_flat_list(get_file_type='trim', write_file_type='norm', overwrite=overwrite)

        # Todo: Paralellize this!
        for img in img_list:
            try:
                normcorr = T80SPreProc()
                log.info('Reading in %s' % img[0])
                normcorr.read('%s' % img[0])

                log.info('Loading configuration from %s' % self.config['trimmed_overscan_config'])
                normcorr.loadConfiguration(self.config['trimmed_overscan_config'])

                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                elif os.path.exists(img[1]) and overwrite:
                    log.debug('Removing existing file %s' % img[1])
                    os.remove(img[1])

                normcorr.norm_section()
                log.debug('Writing %s' % img[1])
                normcorr.ccd.write(img[1])

                # log.info('Reading in %s' % img[0])
                # hdulist = fits.open(img[0])
                # normfactor = np.mean(hdulist[nhdu].data)
                # log.debug('Normalizing')
                # newdata = np.zeros_like(hdulist[nhdu].data, dtype=np.float32)
                # newdata += hdulist[nhdu].data
                # newdata /= normfactor
                # output_hdu = fits.PrimaryHDU(header=hdulist[nhdu].header,
                #                              data=newdata)
                # header_comments = ["NORMALIZE: %s" % datetime.datetime.now(),
                #                    "NORMALIZE: factor = %f" % normfactor]
                # output_hdu.header['NORMFAC'] = normfactor
                #
                # for comment in header_comments:
                #     output_hdu.header["COMMENT"] = comment
                #
                # output_hdulist = fits.HDUList([output_hdu])
                #
                # if not os.path.exists(os.path.dirname(img[1])):
                #     log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                #     os.mkdir(os.path.dirname(img[1]))
                #
                # log.info('Saving output to %s' % img[1])
                # output_hdulist.writeto(img[1], clobber=overwrite)

            except Exception, e:
                log.exception(e)
                continue
            else:
                self.config['calibrations']['sky-flat'][img[2][1]]['norm'].append(os.path.basename(img[1]))

    def flatcorr(self, overwrite=False):

        # Building file list
        img_list = self.get_target_list(get_file_type='trim', write_file_type='flatcorr', overwrite=overwrite)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                path = os.path.join(self.config['path'],
                                    self.config['calibrations']['sky-flat'][img[2][2]]['night'],
                                    'flat',
                                    img[2][2])
                masterflat = os.path.join(path,
                                          'masterflat.fits')
                imarith(img[0], '/', masterflat, img[1], overwrite=overwrite)
            except Exception, e:
                log.exception(e)
                continue
            else:
                self.config['objects'][img[2][1]][img[2][2]]['flatcorr'].append(os.path.basename(img[1]))

    def background_subtraction(self, overwrite=False):
        '''

        :param overwrite:
        :return:
        '''
        # Building file list
        img_list = self.get_target_list(get_file_type='astrometry', write_file_type='backgrd', overwrite=overwrite)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                backgrdcorr = T80SPreProc()
                log.info('Reading in %s' % img[0])
                backgrdcorr.read('%s' % img[0])

                log.info('Loading configuration from %s' % self.config['trimmed_overscan_config'])
                backgrdcorr.loadConfiguration(self.config['trimmed_overscan_config'])

                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    os.mkdir(os.path.dirname(img[1]))
                elif os.path.exists(img[1]) and overwrite:
                    log.debug('Removing existing file %s' % img[1])
                    os.remove(img[1])

                # if there is a segmentation map, use it as a mask
                mask = None
                if os.path.exists(img[0].replace('.fits','.segm.fits')):
                    log.debug('Found segmentation map. Using as source mask...')
                    mask = fits.getdata(img[0].replace('.fits','.segm.fits')) != 0
                box_shape = (56,56)
                filter_shape = (7,7)
                if 'box_shape' in self.config['background']:
                    box_shape = (self.config['background']['box_shape'][0],
                                 self.config['background']['box_shape'][1])
                if 'filter_shape' in self.config['background']:
                    filter_shape = (self.config['background']['filter_shape'][0],
                                    self.config['background']['filter_shape'][1])

                backgrdcorr.subtract_background_section(show=True, box_shape = box_shape,
                                                        filter_shape=filter_shape, mask=mask)
                log.debug('Writing %s' % img[1])
                backgrdcorr.ccd.write(img[1])

            except Exception, e:
                log.exception(e)
                continue
            else:
                self.config['objects'][img[2][1]][img[2][2]]['backgrd'].append(os.path.basename(img[1]))


    def naive_combine(self, image_type='flatcorr', overwrite=False):
        '''
        Combine all images of the objects in a single filter using a naive approach, whitout fixing astrometric

        :param overwrite:
        :return:
        '''

        for object in self.config['objects']:
            for filter in FILTERS:
                naive_name = 'naive_combine_%s_%s.fits' % (object.replace(' ', '_'),
                                                           filter)

                if filter not in self.config['objects'][object]:
                    log.debug('Object %s has no images in filter %s' % (object,
                                                                        filter))
                    continue
                if 'naive' in self.config['objects'][object][filter] and not overwrite:
                    log.warning('%s in %s already naively combined. Run with --overwrite to continue.' % (object,
                                                                                                          filter))
                    continue
                else:
                    self.config['objects'][object][filter]['naive'] = naive_name

                path = os.path.join(self.config['path'],
                                    self.config['objects'][object]['night'],
                                    object.replace(' ', '_'),
                                    filter)

                wpath = os.path.join(self.config['path'],
                                     self.config['objects'][object]['night'],
                                     object.replace(' ', '_'),
                                     filter) if 'wpath' not in self.config else os.path.join(
                    self.config['wpath'],
                    self.config['objects'][object]['night'],
                    object.replace(' ', '_'),
                    filter)

                if not os.path.exists(wpath):
                    log.debug('Working directory %s does not exists, creating.' % path)
                    mkpath(wpath)

                img_list = []
                for raw in self.config['objects'][object][filter][image_type]:
                    img_list.append(os.path.join(path, image_type, raw))
                imcombine(img_list, os.path.join(wpath, naive_name), overwrite=overwrite)

    def astrometry_scamp(self, overwrite=False):
        '''
        Perform astrometric calibration using local installation of astrometry.net, optimized for T80S.

        :return:
        '''

        for obj in self.config['objects']:
            img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                            overwrite=overwrite,
                                            getobject=[obj])

            ref_img = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                           overwrite=overwrite,
                                           getobject=[obj],
                                           getfilter='R')[0]
            log.debug('Solving reference image: %s' % ref_img[0])

            # Run sextractor on image
            sex = SExtractor()

            if 'sex-setup' in self.config:
                sex.setup(self.config['sex-setup'])
            # default params
            sex.config['CONFIG_FILE'] = self.config['sex-config']

            # ok, here we go!
            log.info('Running sextractor')
            # copyfile(rpath[0],wpath[0])
            log.debug('Processing %s...' % ref_img[0])
            sex.run(ref_img[0], updateconfig=False, clean=False)

            header = Scamp.solveField(self.config['sex-catalog-name'])

            # get astrometric solution
            wcs = fits.Header.fromtextfile(header)
            log.debug('Aligning images with reference')
            img1 = fits.getdata(ref_img[0])
            ix, iy = img1.shape
            register_region = [ix / 2 - 1000,ix / 2 + 1000, iy / 2 - 1000,iy / 2 + 1000]
            print self.config['register-region']
            if 'register-region' in self.config:
                register_region[0] = self.config['register-region'][0]
                register_region[1] = self.config['register-region'][1]
                register_region[2] = self.config['register-region'][2]
                register_region[3] = self.config['register-region'][3]
                log.debug('Using user-defined region %s' % register_region)
            img1 = img1[register_region[0]:register_region[1],register_region[2]:register_region[3]]

            # Todo: Paralellize this!
            for i, img in enumerate(img_list):
                try:
                    if not os.path.exists(os.path.dirname(img[1])):
                        log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                        mkpath(os.path.dirname(img[1]))
                    img2 = fits.getdata(img[0])[register_region[0]:register_region[1],
                           register_region[2]:register_region[3]]
                    result = image_registration.chi2_shift(img1, img2)
                    # tvec = result["tvec"].round(4)

                    img2 = fits.open(img[0])
                    wcs = fits.Header.fromtextfile(header)
                    log.debug('Offset: %.3f x %.3f' % (result[0], result[1]))
                    wcs['CRPIX1'] += result[0]
                    wcs['CRPIX2'] += result[1]
                    for key in wcs:
                        if key == 'HISTORY' or key == 'COMMENT':
                            continue
                        img2[0].header[key] = wcs[key]
                    log.debug('Writing %s' % img[1])
                    img2.writeto(img[1], clobber=overwrite)
                    # offset[0][i] = reg[0]
                    # offset[1][i] = reg[1]


                except SExtractorException, e:
                    log.exception(e)
                    return -1
                except Exception, e:
                    log.exception(e)
                    continue
                else:
                    self.config['objects'][img[2][1]][img[2][2]]['astrometry'].append(os.path.basename(img[1]))

    def register(self):
        '''
        Calculate offset between frames on the multiple filters using the first image in R as reference. Store the
        information on the configuration file.

        :return:
        '''

        for obj in self.config['objects']:
            img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='register',
                                            overwrite=True,
                                            getobject=[obj])

            ref_index = 0 if 'astrometry-reference-index' not in self.config['objects'][obj] else \
            self.config['objects'][obj]['astrometry-reference-index']
            ref_filter = 'R' if 'astrometry-reference-filter' not in self.config['objects'][obj] else \
            self.config['objects'][obj]['astrometry-reference-filter']
            ref_img = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                           overwrite=False,
                                           getobject=[obj],
                                           getfilter=ref_filter)[ref_index]

            img1 = fits.getdata(ref_img[0])
            log.debug('Aligning images with reference')
            ix, iy = img1.shape
            register_region = [ix / 2 - 1000,ix / 2 + 1000, iy / 2 - 1000,iy / 2 + 1000]

            if 'register-region' in self.config:
                register_region = self.config['register-region']
                log.debug('Using user-defined region %s' % register_region)

            img1 = img1[register_region[0]:register_region[1], register_region[2]:register_region[3]]
            # self.config['objects'][img[2][1]][img[2][2]]['register'] = {'reference': ref_img[0]}

            # Todo: Paralellize this!
            for i, img in enumerate(img_list):
                log.debug('Computing offset for %s' % img[0])
                if type(self.config['objects'][img[2][1]][img[2][2]]['register']) != type({}):
                    self.config['objects'][img[2][1]][img[2][2]]['register'] = {'reference': ref_img[0]}
                try:
                    # Todo: Try to guess a good grid size from the number of stars.
                    img2 = fits.getdata(img[0])[register_region[0]:register_region[1],
                           register_region[2]:register_region[3]]
                    result = image_registration.chi2_shift(img1, img2)
                    # tvec = result["tvec"].round(4)

                    img2 = fits.open(img[0])
                    # wcs = fits.Header.fromtextfile(header)
                    log.debug('Offset: %.3f x %.3f' % (result[0], result[1]))
                    # wcs['CRPIX1'] += result[0]
                    # wcs['CRPIX2'] += result[1]
                    self.config['objects'][img[2][1]][img[2][2]]['register'][os.path.basename(img[0])] = (
                        float(result[0]), float(result[1]))
                    # self.config['objects'][img[2][1]][img[2][2]]['register'].append(
                    #     (float(result[0]), float(result[1])))
                    # for key in wcs:
                    #     if key == 'HISTORY' or key == 'COMMENT':
                    #         continue
                    #     img2[0].header[key] = wcs[key]
                    # img2[0].header['CRPIX1'] += result[0]
                    # img2[0].header['CRPIX2'] += result[1]
                    #
                    # log.debug('Writing %s' % img[1])
                    # img2.writeto(img[1], clobber=overwrite)
                    # offset[0][i] = reg[0]
                    # offset[1][i] = reg[1]


                # except SExtractorException, e:
                #     log.exception(e)
                #     return -1
                except Exception, e:
                    log.exception(e)
                    continue
                    # break
                    # else:
                    #     self.config['objects'][img[2][1]][img[2][2]]['astrometry'].append(os.path.basename(img[1]))

    def astrometry(self, overwrite=False):
        '''
        Perform astrometric calibration using local installation of astrometry.net, optimized for T80S.

        :return:
        '''

        for obj in self.config['objects']:
            # img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
            #                                 overwrite=overwrite,
            #                                 getobject=[obj])

            ref_index = 0 if 'astrometry-reference-index' not in self.config['objects'][obj] else \
            self.config['objects'][obj]['astrometry-reference-index']
            ref_filter = 'R' if 'astrometry-reference-filter' not in self.config['objects'][obj] else \
            self.config['objects'][obj]['astrometry-reference-filter']
            ref_img = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                           overwrite=overwrite,
                                           getobject=[obj],
                                           getfilter=ref_filter)[ref_index]
            log.debug('Solving reference image: %s' % ref_img[0])

            sex = SExtractor()

            if 'sex-setup' in self.config:
                sex.setup(self.config['sex-setup'])
            # default params
            sex.config['CONFIG_FILE'] = self.config['sex-config']

            # ok, here we go!
            log.info('Running sextractor')
            # copyfile(rpath[0],wpath[0])
            log.debug('Processing %s...' % ref_img[0])
            sex.run(ref_img[0], updateconfig=False, clean=False)

            sex_table = ascii.read(self.config['sex-catalog-name'])

            bad_mask = np.zeros(len(sex_table), dtype=np.bool) == 0
            if 'bad-pixel-mask' in self.config:
                bad_pixel_mask = np.loadtxt(self.config['bad-pixel-mask'])

                for ibad in range(len(bad_pixel_mask)):
                    if bad_pixel_mask[ibad][0] == -1:
                        # Todo: bad line
                        log.debug('%f -> %f' % (bad_pixel_mask[ibad][1] - 1, bad_pixel_mask[ibad][1] + 1))
                        bad_mask = np.bitwise_and(bad_mask,
                                                  np.bitwise_or(
                                                      sex_table['Y_IMAGE'].data < np.float(
                                                          bad_pixel_mask[ibad][1] - 1.5),
                                                      sex_table['Y_IMAGE'].data > np.float(
                                                          bad_pixel_mask[ibad][1] + 1.5)))
                        log.debug('size %i' % len(bad_mask[bad_mask]))
                    elif bad_pixel_mask[ibad][1] == -1:
                        # bad column
                        log.debug('%f -> %f' % (bad_pixel_mask[ibad][0] - 1, bad_pixel_mask[ibad][0] + 1))
                        bad_mask = np.bitwise_and(bad_mask,
                                                  np.bitwise_or(
                                                      sex_table['X_IMAGE'].data < np.float(
                                                          bad_pixel_mask[ibad][0] - 1.5),
                                                      sex_table['X_IMAGE'].data > np.float(
                                                          bad_pixel_mask[ibad][0] + 1.5)))
                        log.debug('size %i' % len(bad_mask[bad_mask]))

            out_table = sex_table[bad_mask]

            # mean_mag = np.mean(out_table['MAG_AUTO'])
            # mask = out_table['MAG_AUTO'] < mean_mag
            # mean_mag_upper = out_table['MAG_AUTO'][mask]
            # mask = np.bitwise_and(mask,
            #                       out_table['MAG_AUTO'] > mean_mag_upper)
            # std_mag = np.std(out_table['MAG_AUTO'])
            # mag_arr = np.abs(out_table['MAG_AUTO'].data - mean_mag + 3*std_mag)
            sort_mag = np.argsort(out_table['MAG_AUTO'].data)[:400]

            axy_table = ref_img[1].replace('.fits', '_axy.fits')
            if not os.path.exists(os.path.dirname(axy_table)):
                log.debug('Path does not exists. Creating %s' % os.path.dirname(axy_table))
                os.mkdir(os.path.dirname(axy_table))
            elif os.path.exists(axy_table) and overwrite:
                os.remove(axy_table)

            out_table[sort_mag].write(axy_table, format='fits')

            AstrometryNet.solveField(ref_img[0],
                                     axy_table)

            # write astrometric solution to a new fits file
            hdu = fits.open(ref_img[0])
            wcs = fits.getheader(ref_img[1].replace('.fits', '_axy-out.wcs'))

            for key in wcs:
                if key != 'COMMENT' and key != 'HISTORY':
                    hdu[0].header[key] = wcs[key]

            log.debug('Writing file %s...' % ref_img[1])
            hdu.writeto(ref_img[1], clobber=overwrite)

            # continue

    def astrometry_align(self, overwrite=True):

        log.debug('Aligning images with reference')

        for obj in self.config['objects']:
            img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                            overwrite=overwrite,
                                            getobject=[obj])

            ref_img = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry',
                                           overwrite=overwrite,
                                           getobject=[obj],
                                           getfilter='R')[0]

            ref_wcs_name = ref_img[1].replace('.fits', '_axy-out.wcs')
            # Check that reference image wcs information exists
            if not os.path.exists(ref_wcs_name):
                log.error('Reference image wcs %s does not exists. Run astrometry before aligning!' % obj)
                continue
            log.debug('Reference WCS %s ...' % ref_wcs_name)
            wcs = fits.getheader(ref_wcs_name)
            # Todo: Paralellize this!
            nimg = 0
            for i, img in enumerate(img_list):
                try:
                    # Check that offset was calculated
                    reg_len = len(self.config['objects'][img[2][1]][img[2][2]]['register'])
                    flatcorr_len = len(self.config['objects'][img[2][1]][img[2][2]]['flatcorr'])

                    if ('register' not in self.config['objects'][img[2][1]][img[2][2]]):
                        log.error('Offsets for %s not properly set! Run register again.' % obj)
                        continue

                    if nimg >= flatcorr_len:
                        nimg = 0
                    # wcs = fits.Header.fromtextfile(header)
                    # log.debug('Offset: %.3f x %.3f' % (result[0],result[1]))
                    # wcs['CRPIX1'] += result[0]
                    # wcs['CRPIX2'] += result[1]
                    img2 = fits.open(img[0])
                    for key in wcs:
                        if key == 'HISTORY' or key == 'COMMENT':
                            continue
                        img2[0].header[key] = wcs[key]
                    crpix1 = float(wcs['CRPIX1'])
                    crpix2 = float(wcs['CRPIX2'])
                    dcrpix1 = float(
                        self.config['objects'][img[2][1]][img[2][2]]['register'][os.path.basename(img[0])][0])
                    dcrpix2 = float(
                        self.config['objects'][img[2][1]][img[2][2]]['register'][os.path.basename(img[0])][1])
                    img2[0].header['CRPIX1'] = crpix1 + dcrpix1
                    img2[0].header['CRPIX2'] = crpix2 + dcrpix2

                    log.debug('%s: %.3f(%.3f) x %.3f(%.3f) = %.3f(%.3f) x %.3f(%.3f)' % (os.path.basename(img[0]),
                                                                                         crpix1,
                                                                                         dcrpix1,
                                                                                         crpix2,
                                                                                         dcrpix2,
                                                                                         img2[0].header['CRPIX1'],
                                                                                         crpix1 + dcrpix1,
                                                                                         img2[0].header['CRPIX2'],
                                                                                         crpix2 + dcrpix2
                                                                                         ))
                    nimg += 1

                    if not os.path.exists(os.path.dirname(img[1])):
                        os.mkdir(os.path.dirname(img[1]))
                    if overwrite and os.path.exists(img[1]):
                        log.debug('Removing %s' % img[1])
                        os.remove(img[1])
                    log.debug('Writing %s' % img[1])
                    img2.writeto(img[1], clobber=overwrite)
                    # offset[0][i] = reg[0]
                    # offset[1][i] = reg[1]


                except SExtractorException, e:
                    log.exception(e)
                    return -1
                except Exception, e:
                    log.exception(e)
                    continue
                else:
                    self.config['objects'][img[2][1]][img[2][2]]['astrometry'].append(os.path.basename(img[1]))

                    # img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry', overwrite=overwrite)
                    #
                    # # Todo: Paralellize this!
                    # for i, img in enumerate(img_list):
                    #     try:
                    #         if not os.path.exists(os.path.dirname(img[1])):
                    #             log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    #             mkpath(os.path.dirname(img[1]))
                    #
                    #         # Run sextractor on image
                    #         sex = SExtractor()
                    #
                    #         if 'sex-setup' in self.config:
                    #             sex.setup(self.config['sex-setup'])
                    #         # default params
                    #         sex.config['CONFIG_FILE'] = self.config['sex-config']
                    #
                    #         # ok, here we go!
                    #         log.info('Running sextractor')
                    #         # copyfile(rpath[0],wpath[0])
                    #         log.debug('Processing %s...' % img[0])
                    #         sex.run(img[0], updateconfig=False, clean=False)
                    #
                    #         sex_table = ascii.read(self.config['sex-catalog-name'])
                    #
                    #         bad_mask = np.zeros(len(sex_table), dtype=np.bool) == 0
                    #         if 'bad-pixel-mask' in self.config:
                    #             bad_pixel_mask = np.loadtxt(self.config['bad-pixel-mask'])
                    #
                    #             for ibad in range(len(bad_pixel_mask)):
                    #                 if bad_pixel_mask[ibad][0] == -1:
                    #                     # Todo: bad line
                    #                     pass
                    #                 elif bad_pixel_mask[ibad][1] == -1:
                    #                     # bad column
                    #                     log.debug('%f -> %f' % (bad_pixel_mask[ibad][0] - 1, bad_pixel_mask[ibad][0] + 1))
                    #                     bad_mask = np.bitwise_and(bad_mask,
                    #                                               np.bitwise_or(
                    #                                                   sex_table['X_IMAGE'].data < np.float(
                    #                                                       bad_pixel_mask[ibad][0] - 1.5),
                    #                                                   sex_table['X_IMAGE'].data > np.float(
                    #                                                       bad_pixel_mask[ibad][0] + 1.5)))
                    #                     log.debug('size %i' % len(bad_mask[bad_mask]))
                    #
                    #         out_table = sex_table[bad_mask]
                    #
                    #         mean_mag = np.mean(out_table['MAG_AUTO'])
                    #         std_mag = np.std(out_table['MAG_AUTO'])
                    #         mag_arr = np.abs(out_table['MAG_AUTO'].data - mean_mag + std_mag)
                    #         sort_mag = np.argsort(mag_arr)[:400]
                    #
                    #         axy_table = img[1].replace('.fits', '_axy.fits')
                    #         out_table[sort_mag].write(axy_table, format='fits')
                    #
                    #         AstrometryNet.solveField(img[0],
                    #                                  axy_table)
                    #
                    #         # write astrometric solution to a new fits file
                    #         hdu = fits.open(img[0])
                    #         wcs = fits.getheader(img[1].replace('.fits', '_axy-out.wcs'))
                    #
                    #         for key in wcs:
                    #             if key != 'COMMENT' and key != 'HISTORY':
                    #                 hdu[0].header[key] = wcs[key]
                    #
                    #
                    #         log.debug('Writing file %s...' % img[1])
                    #         hdu.writeto(img[1])
                    #         hdu.close()
                    #
                    #     except SExtractorException, e:
                    #         log.exception(e)
                    #         return -1
                    #     except Exception, e:
                    #         log.exception(e)
                    #         continue
                    #     else:
                    #         self.config['objects'][img[2][1]][img[2][2]]['astrometry'].append(os.path.basename(img[1]))
                    #

    def coadd(self, overwrite=False):

        for obj in self.config['objects']:
            get_type = 'backgrd' if 'backgrd' in self.config['objects'][obj]['R'] else 'astrometry'
            log.debug('Working on %s images' % get_type)
            ref_img = self.get_target_list(get_file_type=get_type, write_file_type='coadd',
                                           overwrite=overwrite,
                                           getobject=[obj],
                                           getfilter='R')[0]
            ref_hdr = fits.getheader(ref_img[0])
            for fltr in FILTERS:
                img_list = self.get_target_list(get_file_type=get_type, write_file_type='coadd',
                                                overwrite=overwrite,
                                                getobject=[obj],
                                                getfilter=[fltr])
                if len(img_list) == 0:
                    log.warning('No images to combine for %s in %s' % (obj, fltr))
                    continue
                try:
                    log.debug('Coadding %i images...' % len(img_list))
                    wpath = os.path.join(self.config['path'],
                                         self.config['objects'][obj]['night'],
                                         obj.replace(' ', '_'),
                                         fltr) if 'wpath' not in self.config else os.path.join(self.config['wpath'],
                                                                                               self.config['objects'][
                                                                                                   obj]['night'],
                                                                                               obj.replace(' ', '_'),
                                                                                               fltr)
                    coadd_img = '%s_%s' % (obj.replace(" ", "-"), fltr)
                    swarp_coadd = Swarp()
                    with open(swarp_coadd.config['IMAGE_LIST'], 'w') as fp:
                        for img in img_list:
                            fp.write(img[0] + '\n')

                    swarp_coadd.config['CENTER'] = '%s, %s' % (ref_hdr["RA"],
                                                               ref_hdr["DEC"])
                    swarp_coadd.config['CENTER_TYPE'] = 'MANUAL'
                    if 'coadd-weight-map' in self.config and os.path.exists(self.config['coadd-weight-map']):
                        swarp_coadd.config['WEIGHT_IMAGE'] = self.config['coadd-weight-map']
                        swarp_coadd.config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
                        swarp_coadd.config['BLANK_BADPIXELS'] = 'N'
                    swarp_coadd.config['IMAGE_SIZE'] = 10000
                    swarp_coadd.config['IMAGEOUT_NAME'] = os.path.join(wpath, coadd_img + '.swarp.fits')
                    swarp_coadd.config['WEIGHTOUT_NAME'] = os.path.join(wpath, coadd_img + '.weight.fits')
                    if os.path.exists(swarp_coadd.config['IMAGEOUT_NAME']) and not overwrite:
                        log.warning('File already exists, run in overwrite mode!')
                        self.config['objects'][obj][fltr]['coadd'] = coadd_img
                        continue
                    swarp_coadd.run()
                    log.debug('Saved coadded image to %s' % (os.path.join(wpath, coadd_img)))
                except Exception, e:
                    log.exception(e)
                    continue
                else:
                    log.debug('%s' % coadd_img)
                    self.config['objects'][obj][fltr]['coadd'] = coadd_img

    def master_photometry(self, overwrite=False):
        '''
        Compute photometry for each field and generate a single catalog per field using coadded images.

        :return:
        '''

        filter_list = FILTERS
        object_list = self.config['objects']
        photometry_detect_filters = self.config[
            'photometry-detect-filters'] if 'photometry-detection-filters' in self.config else ['R']

        sex = SExtractor()
        sex.config['CONFIG_FILE'] = self.config['sex-config']

        for obj in object_list:
            phot_files = []
            detect_files = []
            for fltr in filter_list:
                rpath = os.path.join(self.config['path'],
                                     self.config['objects'][obj]['night'],
                                     obj.replace(' ', '_'),
                                     fltr)
                coadd_img = '%s_%s' % (obj.replace(" ", "-"), fltr)
                if os.path.exists(os.path.join(rpath, coadd_img + '.fits')):
                    log.debug('Found %s ...' % coadd_img)
                    phot_files.append(os.path.join(rpath, coadd_img + '.fits'))
                    if fltr in photometry_detect_filters:
                        detect_files.append(os.path.join(rpath, coadd_img + '.fits'))
                else:
                    log.warning('File %s not found... Expected to find at %s' % (coadd_img,
                                                                                 rpath))
            if len(phot_files) == 0:
                log.warning('No files to process for object %s. Skipping...' % obj)
                continue
            elif len(detect_files) == 0:
                log.warning('No detection file found for object %s. Skipping...' % obj)
                continue

            rpath = os.path.join(self.config['path'],
                                 self.config['objects'][obj]['night'],
                                 obj.replace(' ', '_'))

            wpath = rpath if 'wpath' not in self.config else os.path.join(self.config['wpath'],
                                                                          self.config['objects'][object]['night'],
                                                                          object.replace(' ', '_'))

            # Creating/selecting detection image

            detect_img = '%s.detection.fits' % (obj.replace(" ", "-"))
            if len(detect_files) > 1:
                if os.path.exists(os.path.join(wpath, detect_img)) and not overwrite:
                    log.warning('Detection image already exists.')
                else:
                    log.debug('Creating detection image, combining %i frames' % len(detect_files))
                    imcombine(detect_files, output=detect_img, overwrite=overwrite)
            else:
                detect_img = detect_files[0]


            # ok, here we go!
            log.info('Running sextractor')
            # copyfile(rpath[0],wpath[0])
            # log.debug('Processing %s...' % ref_img[0])
            for img in phot_files:
                sex.run('%s,%s' % (detect_img,
                                   img), updateconfig=False, clean=False,
                    path=self.config['sex-path'])

                if os.path.exists(self.config['sex-catalog-name']):
                    move(self.config['sex-catalog-name'],
                         os.path.join(wpath,img.replace('.fits', '.cat')))
                    move('test.aper.fits',
                         os.path.join(wpath,img.replace('.fits', '.aper.fits')))
                    move('test.segm.fits',
                         os.path.join(wpath,img.replace('.fits', '.segm.fits')))
                    move('test.backg.fits',
                         os.path.join(wpath,img.replace('.fits', '.backg.fits')))

    def single_photometry(self, objname, overwrite=False):

        img_list = self.get_target_list(get_file_type='astrometry', write_file_type=None,
                                        overwrite=overwrite,
                                        getobject=[objname])

        sex = SExtractor()
        # default params
        with open(self.config['master-photometry-sex-config'], 'r') as fp:
            sex_config = yaml.load(fp)
        for key in sex_config.keys():
            sex.config[key] = sex_config[key]

        sex.config['CONFIG_FILE'] = self.config['sex-config']

        for img in img_list:
            log.debug('Performing photometry in %s' % img[0])
            # copyfile(rpath[0],wpath[0])
            # log.debug('Processing %s...' % ref_img[0])
            # sex.config['CATALOG_NAME'] = img[0].replace('.fits', '.cat')
            # sex.config['CHECKIMAGE_TYPE'] = 'APERTURES,SEGMENTATION'
            # sex.config['CHECKIMAGE_NAME'] = '%s,%s' % (img[0].replace('.fits', '.apert.fits'),
            #                                            img[0].replace('.fits', '.segm.fits'))
            sex.run(img[0], updateconfig=False, clean=False,
                    path=self.config['sex-path'])
            if os.path.exists(self.config['sex-catalog-name']):
                move(self.config['sex-catalog-name'],
                     img[0].replace('.fits', '.cat'))
                move('test.aper.fits',
                     img[0].replace('.fits', '.aper.fits'))
                move('test.segm.fits',
                     img[0].replace('.fits', '.segm.fits'))
                move('test.backg.fits',
                     img[0].replace('.fits', '.backg.fits'))

    def prep_extinction(self,obs_type="EXTMONI"):
        '''
        Prepare catalog to measure zero point and extinction.

        :param obs_type:
        :return:
        '''

        # Determine which targets belongs to extintion monitor
        object_list = []
        cat = []

        for ext_object in self.config['objects']:
            if self.config['objects'][ext_object]['type'] == obs_type:
                # Checking that required information is present
                if 'plus-spec' not in self.config['objects'][ext_object] and \
                                'catalog' not in self.config['objects'][ext_object]['plus-spec'] and \
                                'filter-translation' not in self.config['objects'][ext_object]['plus-spec']:
                    raise IOError('Required information on %s not found. Set source catalog and filter translation '
                                  'information. Skipping...' % ext_object)

                log.debug('Adding %s to extinction catalog.' % ext_object)
                object_list.append(ext_object)

        # gather information on each target
        for ext_object in object_list:
            field_cat = {'calib_spec' : {},
                         'inst_spec' : {}}
            log.debug('Working on %s' % ext_object)
            # open spectrum
            log.debug(
                'Found %i source%s in catalog.' % (len(self.config['objects'][ext_object]['plus-spec']['catalog']),
                                                   '' if len(self.config['objects'][ext_object]['plus-spec'][
                                                                 'catalog']) < 2 else 's'))

            for i,source in enumerate(self.config['objects'][ext_object]['plus-spec']['catalog']):
                field_cat['calib_spec'][str(source['id'])] = {}
                field_cat['inst_spec'][str(source['id'])] = {}
                data = np.load(source['file'])
                for i in range(len(data)):
                    field_cat['calib_spec'][str(source['id'])][str(data['name'][i])] = {'mag': float(data['mag'][i])}
                    field_cat['inst_spec'][str(source['id'])][str(data['name'][i])] = {'mag': [],
                                                               'sigma':  [],
                                                               'secz': []}

                # field_cat['calib_spec'][str(source['id'])] = dict(
                #     [(line['name'], {'mag': float(line['mag']), 'sigma':  0.0, 'secz': 0.0}) for line
                #      in
                #      np.load(source['file'])])
                # field_cat['inst_spec'][str(source['id'])] = dict(
                #     [(line['name'], {'mag': float(line['mag']), 'sigma':  0.0, 'secz': 0.0}) for line
                #      in
                #      np.load(source['file'])])
                # print cat[i]['calib_spec']
                # obs_spec = calib_spec.copy()


            for flt in FILTERS:
                tflt = self.config['objects'][ext_object]['plus-spec']['filter-translation'][flt]

                # cat.append(cat_entry)
                img_list = self.get_target_list(get_file_type='astrometry', write_file_type=None,
                                                overwrite=False,
                                                getobject=[ext_object],
                                                getfilter=[flt])
                for img in img_list:
                    sexcat_name = img[0].replace('.fits', '.cat')
                    if not os.path.exists(sexcat_name):
                        raise IOError('Image/catalog pair not found for %s! Try processing this target with '
                                      'single-photometry or exclude it from the list!' % img[0])

                    log.debug('Found %s image/catalog in %s' % (ext_object,
                                                                os.path.basename(img[0])))

                    sex_table = ascii.read(sexcat_name)
                    sex_catalog = SkyCoord(ra=sex_table['ALPHA_J2000']*u.deg,
                                           dec=sex_table['DELTA_J2000']*u.deg)

                    hdr = fits.getheader(img[0])
                    for isrc, source in enumerate(self.config['objects'][ext_object]['plus-spec']['catalog']):
                        coord = SkyCoord("%s %s" % (source['RA'],source['DEC']),
                                         frame=FK5,
                                         unit=(u.hourangle, u.deg))
                        idx, d2d, d3d = coord.match_to_catalog_sky(sex_catalog) # Todo: Threshold separation
                        field_cat['inst_spec'][source['id']][tflt]['mag'].append(float(sex_table['MAG_AUTO'][idx]))
                        field_cat['inst_spec'][source['id']][tflt]['sigma'].append(float(sex_table['MAGERR_AUTO'][idx]))
                        field_cat['inst_spec'][source['id']][tflt]['secz'].append(float(hdr['AIRMASS']))
                        #
                        log.debug(
                            'Source %i @ (ra,dec): %s, %s (dist. %.2f )Mag. %s: (catalog/observed) '
                            '%.2f/(%.2f/%.2f/%.2f)' % (
                            isrc + 1, source['RA'], source['DEC'], d2d.value, flt,
                            float(field_cat['calib_spec'][source['id']][tflt]['mag']),
                            float(field_cat['inst_spec'][source['id']][tflt]['mag'][-1]),
                            float(field_cat['inst_spec'][source['id']][tflt]['sigma'][-1]),
                            float(field_cat['inst_spec'][source['id']][tflt]['secz'][-1])))
            # cat.append(field_cat)
            self.config['objects'][ext_object]['plus-spec']['ext-cat'] = field_cat
            # self.config['objects'][ext_object]['plus-spec']['catalog'] = cat

        # print cat

    def extinction(self,obs_type="EXTMONI"):
        '''
        Measure zero point and extinction.

        :param obs_type:
        :return:
        '''

        # Determine which targets belongs to extintion monitor
        object_list = []
        cat = []

        if 'extinction' not in self.config:
            self.config['extinction'] = {}

        for ext_object in self.config['objects']:
            if self.config['objects'][ext_object]['type'] == obs_type:
                # Checking that required information is present
                if 'plus-spec' not in self.config['objects'][ext_object] and 'ext-cat' not in \
                        self.config['objects'][ext_object]['plus-spec']:
                    raise IOError('Required information on %s not found. You may need to run prep-extinction before '
                                  'trying to calculate extinction.' % ext_object)


                log.debug('Adding %s to extinction catalog.' % ext_object)
                object_list.append(ext_object)

        # check that catalog was generated
        if len(object_list) == 0:
            raise IOError('No extinction monitor object found in configuration file!')

        # Build data structure

        import pylab as py

        for flt in FILTERS:
            ext_cat = np.array([],dtype=[('mag',np.float),
                                         ('secz',np.float)])
            for ext_object in object_list:
                tflt = self.config['objects'][ext_object]['plus-spec']['filter-translation'][flt]
                for isrc, source in enumerate(self.config['objects'][ext_object]['plus-spec']['catalog']):
                    for itr in range(len(self.config['objects'][ext_object]['plus-spec']['ext-cat']['inst_spec'][
                                             source['id']][tflt]['mag'])):
                        ext_cat = np.append(ext_cat,
                                            np.array((self.config['objects'][ext_object]['plus-spec']['ext-cat'][
                                                          'inst_spec'][
                                                          source['id']][tflt]['mag'][itr] -
                                                      self.config['objects'][ext_object]['plus-spec']['ext-cat'][
                                                          'calib_spec'][
                                                          source['id']][tflt]['mag'],
                                                      self.config['objects'][ext_object]['plus-spec']['ext-cat'][
                                                          'inst_spec'][
                                                          source['id']][tflt]['secz'][itr]),
                                                     dtype=[('mag', np.float),
                                                            ('secz', np.float)]))
            sol = np.polyfit(ext_cat['secz'],
                             ext_cat['mag'],
                             1)
            log.info('%4s Zero Point: %.4f / Ext. coef: %.4f' % (flt, sol[1], sol[0]))
            self.config['extinction'][flt] = {'m0' : float(sol[1]),
                                              'x' : float(sol[0])}
            # py.title(flt)
            # xx = np.linspace(1.0,2.5)
            # py.plot(ext_cat['secz'],
            #         ext_cat['mag'],
            #         'o')
            # py.plot(xx,
            #         np.polyval(sol,xx),
            #         '-')
            # py.show()

    def slr(self, objname, overwrite=False):

        if 'slr-config' not in self.config:
            raise IOError('Required information not found. Define SLR configuration.')

        import pylab as py

        for slr_group in self.config['slr-config']:
            # print slr_group
            group_data = {}
            mask=None
            for id in slr_group:
                if slr_group[id] in group_data:
                    continue

                rpath = os.path.join(self.config['path'],
                                     self.config['objects'][objname]['night'],
                                     objname.replace(' ', '_'),
                                     slr_group[id]) if 'wpath' not in self.config else os.path.join(
                    self.config['wpath'],
                    self.config['objects'][
                        objname]['night'],
                    objname.replace(' ', '_'),
                    slr_group[id])

                catalog = os.path.join(rpath, self.config['objects'][objname][slr_group[id]]['coadd'] + '.cat')
                if not os.path.exists(catalog):
                    raise IOError('Could not find catalog file %s for %s in %s. Try running master-photometry '
                                  'before SLR.' % (catalog, objname, slr_group[id]))

                log.debug('SLR process %s: %s' % (slr_group[id],
                                                  self.config['objects'][objname][slr_group[id]]['coadd']))
                group_data[slr_group[id]] = ascii.read(catalog)

                if mask is None:
                    mask = np.bitwise_and(np.bitwise_and(np.bitwise_and(group_data[slr_group[id]]['FLAGS'] != 0,
                                                         group_data[slr_group[id]]['MAG_AUTO'] < 99.),
                                          group_data[slr_group[id]]['MAG_AUTO'] / group_data[slr_group[id]][
                                              'MAGERR_AUTO'] > 100.),
                                          group_data[slr_group[id]]['CLASS_STAR'] < 0.25)
                else:
                    mask = np.bitwise_and(np.bitwise_and(np.bitwise_and(mask,
                                          np.bitwise_and(group_data[slr_group[id]]['FLAGS'] != 0,
                                          group_data[slr_group[id]]['MAG_AUTO'] < 99.)),
                                          group_data[slr_group[id]]['MAG_AUTO'] / group_data[slr_group[id]][
                                              'MAGERR_AUTO'] > 100.),
                                          group_data[slr_group[id]]['CLASS_STAR'] < 0.25)

            py.xlabel('%s - %s' % (slr_group[3],slr_group[4]))
            py.ylabel('%s - %s' % (slr_group[1],slr_group[2]))

            c0 = 0.
            c1 = 0.
            if 'extinction' in self.config:
                c1 = self.config['extinction'][slr_group[1]]['m0']
                c2 = self.config['extinction'][slr_group[2]]['m0']
                c3 = self.config['extinction'][slr_group[3]]['m0']
                c4 = self.config['extinction'][slr_group[4]]['m0']

            py.errorbar(x=group_data[slr_group[3]]['MAG_AUTO'][mask]-c3 - group_data[slr_group[4]]['MAG_AUTO'][mask] +c4 ,
                        y=group_data[slr_group[1]]['MAG_AUTO'][mask] -c1 - group_data[slr_group[2]]['MAG_AUTO'][mask] + c2,
                        xerr=group_data[slr_group[3]]['MAGERR_AUTO'][mask] + group_data[slr_group[4]]['MAGERR_AUTO'][
                            mask],
                        yerr=group_data[slr_group[1]]['MAGERR_AUTO'][mask] + group_data[slr_group[2]]['MAGERR_AUTO'][
                            mask],
                        fmt='ro',
                        capsize=0,
                        ecolor='0.85',
                        alpha=0.5)
            if 'isochrones' in self.config:
                with open(self.config['isochrones']['file']) as fp:
                    iso = IsochroneSet(fp)
                py.plot(iso.isochrones[0][self.config['isochrones']['filter-translation'][slr_group[3]]] -
                        iso.isochrones[0][self.config['isochrones']['filter-translation'][slr_group[4]]],
                        iso.isochrones[0][self.config['isochrones']['filter-translation'][slr_group[1]]] -
                        iso.isochrones[0][self.config['isochrones']['filter-translation'][slr_group[2]]],
                        'k-',
                        lw=2)

            py.show()

        # img_list = self.get_target_list(get_file_type='astrometry', write_file_type=None,
        #                                 overwrite=overwrite,
        #                                 getobject=[objname])
        #
        # sex = SExtractor()
        # # default params
        # with open(self.config['master-photometry-sex-config'], 'r') as fp:
        #     sex_config = yaml.load(fp)
        # for key in sex_config.keys():
        #     sex.config[key] = sex_config[key]
        #
        # sex.config['CONFIG_FILE'] = self.config['sex-config']
        #
        # for img in img_list:
        #     log.debug('Performing photometry in %s' % img[0])
        #     # copyfile(rpath[0],wpath[0])
        #     # log.debug('Processing %s...' % ref_img[0])
        #     # sex.config['CATALOG_NAME'] = img[0].replace('.fits', '.cat')
        #     # sex.config['CHECKIMAGE_TYPE'] = 'APERTURES,SEGMENTATION'
        #     # sex.config['CHECKIMAGE_NAME'] = '%s,%s' % (img[0].replace('.fits', '.apert.fits'),
        #     #                                            img[0].replace('.fits', '.segm.fits'))
        #     sex.run(img[0], updateconfig=False, clean=False,
        #             path=self.config['sex-path'])
        #     if os.path.exists(self.config['sex-catalog-name']):
        #         move(self.config['sex-catalog-name'],
        #              img[0].replace('.fits', '.cat'))
        #         move('test.aper.fits',
        #              img[0].replace('.fits', '.aper.fits'))
        #         move('test.segm.fits',
        #              img[0].replace('.fits', '.segm.fits'))
        #         move('test.backg.fits',
        #              img[0].replace('.fits', '.backg.fits'))