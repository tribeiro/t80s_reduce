
import numpy as np
from astropy.io import fits
import datetime as dt
import logging

log = logging.getLogger(__name__)

def imcombine(input,output,nhdu=0,overwrite=False):

    nimages = len(input)
    hdu = fits.open(input[0])
    sizex, sizey = hdu[nhdu].data.shape
    img_stack = np.zeros((nimages,sizex,sizey),dtype=hdu[nhdu].data.dtype)

    header_comments = ["IMCOMBINE: %s"%dt.datetime.now(),
                       "IMCOMBINE: Combining %i images with median algorithm"%nimages,
                       "IMCOMBINE: IMAGE    MEAN    MIN    MAX    STD"]

    for i in range(nimages):
        log.debug('Reading in %s'%input[i])
        img_stack[i] += fits.getdata(input[i],ext = nhdu)
        header_comments.append("IMCOMBINE: %s  %.2f  %.2f  %.sf  %.2f"%(input[i],
                                                                        np.mean(img_stack[i]),
                                                                        np.min(img_stack[i]),
                                                                        np.max(img_stack[i]),
                                                                        np.std(img_stack[i])))

    log.info("Combining %i images..." % nimages)
    output_hdu = fits.PrimaryHDU(np.median(img_stack,axis=0))

    for comment in header_comments:
        output_hdu.header["COMMENT"] = comment
    output_hdulist = fits.HDUList([output_hdu])
    #hdu[0].data = np.median(img_stack,axis=0)

    log.info('Saving output to %s' % output)
    output_hdulist.writeto(output,clobber=overwrite)

    return 0