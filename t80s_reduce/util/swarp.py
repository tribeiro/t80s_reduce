import __builtin__




import os
import copy
from subprocess import Popen
import logging
import time
from shutil import copyfile

log = logging.getLogger(__name__)


class Swarp:

    _SW_config = {'IMAGEOUT_NAME': {'value': 'coadd.fits', 'comment': 'Output filename'},
                  'WEIGHTOUT_NAME': {'value': 'coadd.weight.fits', 'comment': 'Output weight-map filename'},
                  'HEADER_ONLY': {'value': 'N', 'comment': 'Only a header as an output file (Y/N)?'},
                  'HEADER_SUFFIX': {'value': '.head', 'comment': 'Filename extension for additional headers'},
                  'WEIGHT_TYPE': {'value': 'NONE', 'comment': 'BACKGROUND,MAP_RMS,MAP_VARIANCE or MAP_WEIGHT'},
                  'WEIGHT_SUFFIX': {'value': '.weight.fits', 'comment': 'Suffix to use for weight-maps'},
                  'WEIGHT_IMAGE': {'value': '#',
                                   'comment': 'Weightmap filename if suffix not used (all or for each weight-map)'},
                  'COMBINE': {'value': 'Y', 'comment': 'Combine resampled images (Y/N)?'},
                  'COMBINE_TYPE': {'value': 'MEDIAN',
                                   'comment': 'MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED,CHI-OLD,CHI-MODE,CHI-MEAN,SUM,'
                                              'WEIGHTED_WEIGHT,MEDIAN_WEIGHT,AND,NAND,OR or NOR'},
                  'CELESTIAL_TYPE': {'value': 'NATIVE',
                                     'comment': 'NATIVE, PIXEL, EQUATORIAL,GALACTIC,ECLIPTIC, or SUPERGALACTIC'},
                  'PROJECTION_TYPE': {'value': 'TAN', 'comment': 'Any WCS projection code or NONE'},
                  'PROJECTION_ERR': {'value': '0.001',
                                     'comment': 'Maximum projection error (in output pixels), or 0 for no '
                                                'approximation'},
                  'CENTER_TYPE': {'value': 'ALL', 'comment': 'MANUAL, ALL or MOST'},
                  'CENTER': {'value': '00:00:00.0, +00:00:00.0', 'comment': 'Coordinates of the image center'},
                  'PIXELSCALE_TYPE': {'value': 'MEDIAN', 'comment': 'MANUAL,FIT,MIN,MAX or MEDIAN'},
                  'PIXEL_SCALE': {'value': '0.0', 'comment': 'Pixel scale'},
                  'IMAGE_SIZE': {'value': '0', 'comment': 'Image size (0 = AUTOMATIC)'},
                  'RESAMPLE': {'value': 'Y', 'comment': 'Resample input images (Y/N)?'},
                  'RESAMPLE_DIR': {'value': '.', 'comment': 'Directory path for resampled images'},
                  'RESAMPLE_SUFFIX': {'value': '.resamp.fits', 'comment': 'filename extension for resampled images'},
                  'RESAMPLING_TYPE': {'value': 'LANCZOS3',
                                      'comment': 'NEAREST,BILINEAR,LANCZOS2,LANCZOS3,LANCZOS4 (1 per axis) or FLAGS'},
                  'OVERSAMPLING': {'value': '0', 'comment': 'Oversampling in each dimension (0 = automatic)'},
                  'INTERPOLATE': {'value': 'N',
                                  'comment': 'Interpolate bad input pixels (Y/N)? (all or for each image)'},
                  'FSCALASTRO_TYPE': {'value': 'FIXED', 'comment': 'NONE,FIXED, or VARIABLE'},
                  'FSCALE_KEYWORD': {'value': 'FLXSCALE',
                                     'comment': 'FITS keyword for the multiplicative factor applied to each input '
                                                'image'},
                  'FSCALE_DEFAULT': {'value': '1.0', 'comment': 'Default FSCALE value if not in header'},
                  'GAIN_KEYWORD': {'value': 'GAIN', 'comment': 'FITS keyword for effect. gain (e-/ADU)'},
                  'GAIN_DEFAULT': {'value': '0.0', 'comment': 'Default gain if no FITS keyword found'},
                  'SUBTRACT_BACK': {'value': 'Y',
                                    'comment': 'Subtraction sky background (Y/N)? (all or for each image)'},
                  'BACK_TYPE': {'value': 'AUTO', 'comment': 'AUTO or MANUAL (all or for each image)'},
                  'BACK_DEFAULT': {'value': '0.0',
                                   'comment': 'Default background value in MANUAL (all or for each image)'},
                  'BACK_SIZE': {'value': '128', 'comment': 'Background mesh size (pixels) (all or for each image)'},
                  'BACK_FILTERSIZE': {'value': '3',
                                      'comment': 'Background map filter range (meshes) (all or for each image)'},
                  'VMEM_DIR': {'value': '.', 'comment': 'Directory path for swap files'},
                  'VMEM_MAX': {'value': '2047', 'comment': 'Maximum amount of virtual memory (MB)'},
                  'MEM_MAX': {'value': '256', 'comment': 'Maximum amount of usable RAM (MB)'},
                  'COMBINE_BUFSIZE': {'value': '256', 'comment': 'RAM dedicated to co-addition(MB)'},
                  'DELETE_TMPFILES': {'value': 'Y', 'comment': 'Delete temporary resampled FITS files (Y/N)?'},
                  'COPY_KEYWORDS': {'value': 'OBJECT',
                                    'comment': 'List of FITS keywords to propagate from the input to the output '
                                               'headers'},
                  'WRITE_FILEINFO': {'value': 'N',
                                     'comment': 'Write information about each input file in the output image header?'},
                  'WRITE_XML': {'value': 'Y', 'comment': 'Write XML file (Y/N)?'},
                  'XML_NAME': {'value': 'swarp.xml', 'comment': 'Filename for XML output'},
                  'VERBOSE_TYPE': {'value': 'NORMAL', 'comment': 'QUIET,LOG,NORMAL, or FULL'},
                  'NTHREADS': {'value': '0',
                               'comment': 'Number of simultaneous threads for the SMP version of SWarp 0 = automatic'},
                  "CONFIG_FILE":
                      {"comment": '[Extra key] name of the main configuration file',
                       "value": "swarp.config"},
                  "IMAGE_LIST":
                      {"comment": '[Extra key] name of the image list file',
                       "value": "swarp_images.lis"},

                  }

    _SW_config_special_keys = ["CONFIG_FILE", "IMAGE_LIST"]

    def __init__(self):
        """
        SExtractor class constructor.
        """

        self.config = (
            dict([(k, copy.deepcopy(Swarp._SW_config[k]["value"]))
                  for k in Swarp._SW_config.keys()]))

        # print self.config

        self.program = None
        self.version = None

    def update_config(self):
        """
        Update the configuration files according to the current
        in-memory SExtractor configuration.
        """

        # -- Write main configuration file
        main_f = __builtin__.open(self.config['CONFIG_FILE'], 'w')

        for key in self.config.keys():
            if (key in Swarp._SW_config_special_keys):
                continue

            value = str(self.config[key])

            print >>main_f, ("%-16s       %-16s # %s" %
                             (key, value, Swarp._SW_config[key]['comment']))

        main_f.close()

    def run(self, updateconfig=True, clean=False, path=None):
        """
        Run Swarp.

        If updateconfig is True (default), the configuration
        files will be updated before running Swarp.

        If clean is True (default: False), configuration files
        (if any) will be deleted after Swarp terminates.

        """

        if updateconfig:
            self.update_config()

        line = "swarp @%s -c %s" % (self.config['IMAGE_LIST'],
                                   self.config['CONFIG_FILE'])

        log.debug("RUN: %s" % line)
        # *** it would be nice to add a test here to check
        # whether astrometrynet is running OK, if not raise a new exception
        # like AstrometryNetInstallProblem
        t0 = time.time()
        solve = Popen(line.split())  # ,env=os.environ)
        solve.wait()
        log.debug('Swarp finished. Took %3.2f sec' % (time.time() - t0))

        # if solution failed, there will be no file .solved
        # wcs_filename = basefilename + ".head"
        if not os.path.exists(self.config['IMAGEOUT_NAME']):
            raise SwarpException("Result image %s not found!" % self.config['IMAGEOUT_NAME'])


        return True


class SwarpException(Exception):
    pass

