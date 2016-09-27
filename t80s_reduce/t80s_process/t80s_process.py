import os
import numpy as np
import logging
import yaml
from astropy.io import fits, ascii
import datetime
from distutils.dir_util import mkpath
from shutil import copyfile

from t80s_reduce.core.constants import *
from t80s_reduce.util.imcombine import imcombine
from t80s_reduce.util.imarith import imarith
from t80s_reduce.util.sextractor import SExtractor, SExtractorException
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
        with open(self._config_file, 'w') as fp:
            yaml.dump(self.config, fp, default_flow_style=False)

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

    def get_target_list(self, get_file_type='raw files', write_file_type=None, overwrite=False, getfilter=None):
        config = self.config
        img_list = []

        filter_list = getfilter if getfilter is not None else FILTERS

        for object in config['objects']:
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
                                'Checking files to process.')
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
                    img_list.append((os.path.join(rpath, get_file_type.replace(' ', '_'), raw),
                                     os.path.join(wpath, write_file_type.replace(' ', '_'), raw),
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
        elif key == 'master-flat':
            for filter in self.config['calibrations']['sky-flat']['filters']:
                if 'master' in self.config['calibrations']['sky-flat'][filter] and not overwrite:
                    log.warning('Master flat already created. Run in overwrite mode to continue.')
                    continue
                log.debug('Creating master flat in %s' % filter)
                imglist = self.get_flat_list(get_file_type='norm', getfilter=[filter])
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
                imarith(img[0], '-', self.masterbias, img[1], overwrite=overwrite)
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

    def astrometry(self, overwrite=False):
        '''
        Perform astrometric calibration using local installation of astrometry.net, optimized for T80S.

        :return:
        '''

        img_list = self.get_target_list(get_file_type='flatcorr', write_file_type='astrometry', overwrite=overwrite)

        # Todo: Paralellize this!
        for i, img in enumerate(img_list):
            try:
                if not os.path.exists(os.path.dirname(img[1])):
                    log.debug('Creating parent directory %s' % os.path.dirname(img[1]))
                    mkpath(os.path.dirname(img[1]))

                # Run sextractor on image
                sex = SExtractor()

                if 'sex-setup' in self.config:
                    sex.setup(self.config['sex-setup'])
                # default params
                sex.config['CONFIG_FILE'] = self.config['sex-config']

                # ok, here we go!
                log.info('Running sextractor')
                # copyfile(rpath[0],wpath[0])
                log.debug('Processing %s...' % img[0])
                sex.run(img[0], updateconfig=False, clean=False)

                sex_table = ascii.read(self.config['sex-catalog-name'])

                bad_mask = np.zeros(len(sex_table), dtype=np.bool) == 0
                if 'bad-pixel-mask' in self.config:
                    bad_pixel_mask = np.loadtxt(self.config['bad-pixel-mask'])

                    for ibad in range(len(bad_pixel_mask)):
                        if bad_pixel_mask[ibad][0] == -1:
                            # Todo: bad line
                            pass
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

                mean_mag = np.mean(out_table['MAG_AUTO'])
                std_mag = np.std(out_table['MAG_AUTO'])
                mag_arr = np.abs(out_table['MAG_AUTO'].data - mean_mag + std_mag)
                sort_mag = np.argsort(mag_arr)[:400]

                axy_table = img[1].replace('.fits', '_axy.fits')
                out_table[sort_mag].write(axy_table, format='fits')

                AstrometryNet.solveField(img[0],
                                         axy_table)

                # write astrometric solution to a new fits file
                hdu = fits.open(img[0])
                wcs = fits.getheader(img[1].replace('.fits', '_axy-out.wcs'))

                for key in wcs:
                    if key != 'COMMENT' and key != 'HISTORY':
                        hdu[0].header[key] = wcs[key]


                log.debug('Writing file %s...' % img[1])
                hdu.writeto(img[1])
                hdu.close()

            except SExtractorException, e:
                log.exception(e)
                return -1
            except Exception, e:
                log.exception(e)
                continue
            else:
                self.config['objects'][img[2][1]][img[2][2]]['astrometry'].append(os.path.basename(img[1]))

