from subprocess import Popen
import os
import logging
import time
from shutil import copyfile

log = logging.getLogger(__name__)


class Swarp:

    _SW_config = '''# Default configuration file for SWarp 2.38.0
# EB 2016-09-27
#
#----------------------------------- Output -----------------------------------
IMAGEOUT_NAME          coadd.fits      # Output filename
WEIGHTOUT_NAME       coadd.weight.fits # Output weight-map filename

HEADER_ONLY            N               # Only a header as an output file (Y/N)?
HEADER_SUFFIX          .head           # Filename extension for additional headers

#------------------------------- Input Weights --------------------------------

WEIGHT_TYPE            NONE            # BACKGROUND,MAP_RMS,MAP_VARIANCE
                                       # or MAP_WEIGHT
WEIGHT_SUFFIX          .weight.fits    # Suffix to use for weight-maps
WEIGHT_IMAGE                           # Weightmap filename if suffix not used
                                       # (all or for each weight-map)

#------------------------------- Co-addition ----------------------------------

COMBINE                Y               # Combine resampled images (Y/N)?
COMBINE_TYPE           MEDIAN          # MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CLIPPED
                                       # CHI-OLD,CHI-MODE,CHI-MEAN,SUM,
                                       # WEIGHTED_WEIGHT,MEDIAN_WEIGHT,
                                       # AND,NAND,OR or NOR

#-------------------------------- Astrometry ----------------------------------

CELESTIAL_TYPE         NATIVE          # NATIVE, PIXEL, EQUATORIAL,
                                       # GALACTIC,ECLIPTIC, or SUPERGALACTIC
PROJECTION_TYPE        TAN             # Any WCS projection code or NONE
PROJECTION_ERR         0.001           # Maximum projection error (in output
                                       # pixels), or 0 for no approximation
CENTER_TYPE            MANUAL          # MANUAL, ALL or MOST
CENTER         %(ra)s, %(dec)s         # Coordinates of the image center
PIXELSCALE_TYPE        MEDIAN          # MANUAL,FIT,MIN,MAX or MEDIAN
PIXEL_SCALE            0.0             # Pixel scale
IMAGE_SIZE             %(size)s        # Image size (0 = AUTOMATIC)

#-------------------------------- Resampling ----------------------------------

RESAMPLE               Y               # Resample input images (Y/N)?
RESAMPLE_DIR           .               # Directory path for resampled images
RESAMPLE_SUFFIX        .resamp.fits    # filename extension for resampled images

RESAMPLING_TYPE        LANCZOS3        # NEAREST,BILINEAR,LANCZOS2,LANCZOS3
                                       # LANCZOS4 (1 per axis) or FLAGS
OVERSAMPLING           0               # Oversampling in each dimension
                                       # (0 = automatic)
INTERPOLATE            N               # Interpolate bad input pixels (Y/N)?
                                       # (all or for each image)

FSCALASTRO_TYPE        FIXED           # NONE,FIXED, or VARIABLE
FSCALE_KEYWORD         FLXSCALE        # FITS keyword for the multiplicative
                                       # factor applied to each input image
FSCALE_DEFAULT         1.0             # Default FSCALE value if not in header

GAIN_KEYWORD           GAIN            # FITS keyword for effect. gain (e-/ADU)
GAIN_DEFAULT           0.0             # Default gain if no FITS keyword found

#--------------------------- Background subtraction ---------------------------

SUBTRACT_BACK          Y               # Subtraction sky background (Y/N)?
                                       # (all or for each image)

BACK_TYPE              AUTO            # AUTO or MANUAL
                                       # (all or for each image)
BACK_DEFAULT           0.0             # Default background value in MANUAL
                                       # (all or for each image)
BACK_SIZE              128             # Background mesh size (pixels)
                                       # (all or for each image)
BACK_FILTERSIZE        3               # Background map filter range (meshes)
                                       # (all or for each image)

#------------------------------ Memory management -----------------------------

VMEM_DIR               .               # Directory path for swap files
VMEM_MAX               2047            # Maximum amount of virtual memory (MB)
MEM_MAX                256             # Maximum amount of usable RAM (MB)
COMBINE_BUFSIZE        256             # RAM dedicated to co-addition(MB)

#------------------------------ Miscellaneous ---------------------------------

DELETE_TMPFILES        Y               # Delete temporary resampled FITS files
                                       # (Y/N)?
COPY_KEYWORDS          OBJECT          # List of FITS keywords to propagate
                                       # from the input to the output headers
WRITE_FILEINFO         N               # Write information about each input
                                       # file in the output image header?
WRITE_XML              Y               # Write XML file (Y/N)?
XML_NAME               swarp.xml       # Filename for XML output
VERBOSE_TYPE           NORMAL          # QUIET,LOG,NORMAL, or FULL

NTHREADS               0               # Number of simultaneous threads for
                                       # the SMP version of SWarp
                                       # 0 = automatic
    '''

    # staticmethod allows to use a single method of a class
    def coadd(self, img_list, path, root, ra, dec, size):
        """
        @param: fullfilename entire path to image
        @type: str

        @param: findstarmethod (astrometry.net, sex)
        @type: str

        Does astrometry to image=fullfilename
        Uses either astrometry.net or sex(tractor) as its star finder
        """

        # pathname, filename = os.path.split(sexcat_filename)
        # pathname = pathname + "/"
        # basefilename, file_xtn = os.path.splitext(filename)

        config = os.path.join(path,root+"swarp.config")
        with open(config, "w") as fp:
            fp.write(self._SW_config % {"ra" : ra,
                                        "dec" : dec,
                                        "size" : size})
        img_list_str = ''
        for img in img_list:
            img_list_str += img
            img_list_str += ' '
        line = "swarp %s -c %s" % (img_list_str,
                                   config)

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
        coadd = 'coadd.fits'
        weight = 'coadd.weight.fits'
        copyfile(coadd, os.path.join(path,root+'.fits'))
        copyfile(coadd, os.path.join(path,root+'.weight.fits'))

        return True


class ScampException(Exception):
    pass


class NoSolutionScampException(Exception):
    pass
