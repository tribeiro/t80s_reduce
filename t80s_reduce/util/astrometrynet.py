from subprocess import Popen
import os
import logging
import time
from astropy.io import fits

log = logging.getLogger(__name__)

class AstrometryNet:
    # staticmethod allows to use a single method of a class
    @staticmethod
    def solveField(fullfilename, sexoutfilename):
        """
        @param: fullfilename entire path to image
        @type: str

        @param: findstarmethod (astrometry.net, sex)
        @type: str

        Does astrometry to image=fullfilename
        Uses either astrometry.net or sex(tractor) as its star finder
        """

        pathname, filename = os.path.split(sexoutfilename)
        pathname = pathname + "/"
        basefilename, file_xtn = os.path.splitext(filename)
        # *** enforce .fits extension
        if file_xtn != ".fits":
            raise ValueError("File extension must be .fits it was = %s\n" % file_xtn)

        # *** check whether the file exists or not
        if os.path.exists(fullfilename) == False:
            raise IOError("You selected image %s  It does not exist\n" % fullfilename)

        # version 0.23 changed behavior of --overwrite
        # I need to specify an output filename with -o
        outfilename = basefilename + "-out"

        header = fits.getheader(fullfilename)
        try:
            ra = header["CRVAL1"]  # expects to see this in image
        except:
            raise AstrometryNetException("Need CRVAL1 and CRVAL2 and CD1_1 on header")
        try:
            dec = header["CRVAL2"]
        except:
            raise AstrometryNetException("Need CRVAL1 and CRVAL2 and CD1_1 on header")
        width = header["NAXIS1"]
        height = header["NAXIS2"]
        radius = 10.0 * abs(header["CD1_1"]) * width

        wcs_filename = pathname + outfilename + ".wcs"

        # FIXME: fixed some parameters!
        line = "solve-field %s --no-plots --overwrite -o %s --x-column X_IMAGE --y-column Y_IMAGE " \
               "--sort-column MAG_AUTO --sort-ascending --width %d --height %d --ra %f --dec %f --radius 0.5 " \
               "-L 0.5 -H 0.6 -u app --over -t 6 --crpix-center" % (
                   sexoutfilename, outfilename, width, height, ra, dec)


        # when there is a solution astrometry.net creates a file with .solved
        # added as extension.
        is_solved = pathname + outfilename + ".solved"
        # if it is already there, make sure to delete it
        if os.path.exists(is_solved):
            os.remove(is_solved)
        log.debug("SOLVE %s" % line)
        # *** it would be nice to add a test here to check
        # whether astrometrynet is running OK, if not raise a new exception
        # like AstrometryNetInstallProblem
        log.debug('Starting solve-field...')
        t0 = time.time()
        solve = Popen(line.split())  # ,env=os.environ)
        solve.wait()
        log.debug('Solve field finished. Took %3.2f sec' % (time.time() - t0))
        # if solution failed, there will be no file .solved
        if (os.path.exists(is_solved) == False):
            raise NoSolutionAstrometryNetException(
                "Astrometry.net could not find a solution for image: %s %s" % (fullfilename, is_solved))

        return wcs_filename


class AstrometryNetException(Exception):
    pass


class NoSolutionAstrometryNetException(Exception):
    pass
