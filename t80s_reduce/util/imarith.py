import os
import numpy as np
from astropy.io import fits
import datetime as dt
import logging

log = logging.getLogger(__name__)

def imarith(im1,op,im2,output,nhdu1=0,nhdu2=0,overwrite=False):

    log.debug("Reading image 1: %s" % im1)
    data1 = fits.getdata(im1,ext = nhdu1)
    
    if os.path.exists(im2):
        log.debug("Reading image 2: %s" % im2)
        data2 = fits.getdata(im2, ext = nhdu2)
    else:
        data2 = np.float(im2)
        log.debug("Value: %s" % im2)

    if not op in '+-*/':
        raise IOError("Unrecognized operator '%s'. Must be one of +, -, * or /"% op)
    elif op == '+':
        output_data = data1 + data2
    elif op == '-':
        output_data = data1 - data2
    elif op == '*':
        output_data = data1 * data2
    elif op == '/':
        output_data = data1 / data2
    else:
        raise IOError("Unrecognized operator '%s'. Must be one of +, -, * or /"%(op))

    output_hdu = fits.PrimaryHDU(header=fits.getheader(im1),
                                 data=output_data)
    header_comments = ["IMARITH: %s" % dt.datetime.now(),
                       "IMARITH: %s %s %s"%(im1,op,im2)]

    for comment in header_comments:
        output_hdu.header["COMMENT"] = comment
    output_hdulist = fits.HDUList([output_hdu])

    log.info('Saving output to %s' % output)
    output_hdulist.writeto(output,clobber = overwrite)

    return 0
