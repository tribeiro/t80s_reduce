import sys, os
import numpy as np
import logging
from datetime import datetime
import yaml
from t80s_reduce.core.constants import *

# from t80s_reduce.util.get_file_list import get_flat_list,get_target_list
from t80s_reduce.t80s_process.t80s_process import T80SProcess


def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description='''
Given a configuration file, perform a set of processing steps. This script is responsible for all the heavy lifting on
data reduction pre-processing steps. Using information on your configuration file it will process bias, flats and
science frames up to the point until you get a final processing image, ready for higher order routines.

This script will do:
 - combine raw bias into a master bias
 - subtract a master bias from science and flatfield images
 - apply overscan correction (TBD)
 - trim overscan regions
 - linearize images
 - produce master flats per filters
 - apply flatfield correction to science frames
 - convert from ADU to electrons

    ''')

    parser.add_argument('-c', '--config',
                        help='Configuration file. The file will be over-written at the end with the result from the '
                             'checks.',
                        type=str)
    parser.add_argument('--action',
                        help='Perform one of the processing options.',
                        choices=['master-bias', 'biascorr', 'overcorr', 'trim', 'linearize', 'norm-flat',
                                 'master-flat', 'flatcorr', 'adu2e', 'naive-combine', 'astrometry'])
    parser.add_argument("--overwrite", action="store_true",
                        help='Overwrite existing processed images.')

    parser.add_argument('--verbose', '-v', action='count')

    args = parser.parse_args(argv[1:])

    logging.basicConfig(
        format='%(asctime)s %(levelname)s %(name)s %(filename)s:%(lineno)d %(message)s',
        level=args.verbose)

    log = logging.getLogger(__name__)

    process = T80SProcess(args.config)

    if args.action == 'master-bias':
        log.debug('Combining bias into master bias')
        # check that a master bias is not present and, in case it is, if we are running in overwrite mode
        if process.masterbias is not None and not args.overwrite:
            log.error('Master bias already created. Run in overwrite mode to continue.')
            return -1

        # Get the list of bias frames
        process.imcombine('master-bias', overwrite=args.overwrite)

    elif args.action == 'biascorr':
        log.debug('Subtracting bias from images.')
        process.biascorr(args.overwrite)
    elif args.action == 'overcorr':
        log.debug('Applying overscan correction.')
        log.warning('Overscan correction not fully implemented yet.')
        process.overcorr(args.overwrite)
    elif args.action == 'linearize':
        log.debug('Linearizing images.')
        process.linearize(args.overwrite)
    elif args.action == 'trim':
        log.debug('Trimming images.')
        process.trim(args.overwrite)
    elif args.action == 'norm-flat':
        log.debug('Normalize flat field images.')
        process.normalize_flatfield(overwrite=args.overwrite)
    elif args.action == 'master-flat':
        log.debug('Creating master flat per filter.')
        process.imcombine('master-flat', overwrite=args.overwrite)
    elif args.action == 'flatcorr':
        log.debug('Applying master flat correction.')
        process.flatcorr(args.overwrite)
    elif args.action == 'adu2e':
        log.debug('Converting from adu to electrons.')
    elif args.action == 'naive-combine':
        log.debug('Combining images per filter with a naive algorithm, ignoring astrometric variations.')
        process.naive_combine(overwrite=args.overwrite)
    elif args.action == 'astrometry':
        log.debug('Astrometric calibration.')
        process.astrometry(overwrite=args.overwrite)
    else:
        log.error('No such option %s' % args.action)

    return 0
