from subprocess import Popen
import os
import logging
import time

log = logging.getLogger(__name__)

class Scamp:

    # staticmethod allows to use a single method of a class
    @staticmethod
    def solveField(sexcat_filename):
        """
        @param: fullfilename entire path to image
        @type: str

        @param: findstarmethod (astrometry.net, sex)
        @type: str

        Does astrometry to image=fullfilename
        Uses either astrometry.net or sex(tractor) as its star finder
        """

        pathname, filename = os.path.split(sexcat_filename)
        pathname = pathname + "/"
        basefilename, file_xtn = os.path.splitext(filename)

        # FIXME: fixed some parameters!
        line = "scamp %s " % (sexcat_filename)

        log.debug("SOLVE %s" % line)
        # *** it would be nice to add a test here to check
        # whether astrometrynet is running OK, if not raise a new exception
        # like AstrometryNetInstallProblem
        t0 = time.time()
        solve = Popen(line.split())  # ,env=os.environ)
        solve.wait()
        log.debug('Solve field finished. Took %3.2f sec' % (time.time() - t0))
        # if solution failed, there will be no file .solved
        wcs_filename = basefilename + ".head"

        return wcs_filename


class ScampException(Exception):
    pass


class NoSolutionScampException(Exception):
    pass
