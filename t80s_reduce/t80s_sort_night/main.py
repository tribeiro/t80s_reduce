import sys, os
import numpy as np
import logging
from pymongo import MongoClient
import time
import yaml


def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description='Sort files from an observing night and generate configuration'
                                                 'scripts.')

    parser.add_argument('--host',
                        help='The host of the mongo database. Mandatory if not ',
                        type=str,
                        default='192.168.20.118')

    parser.add_argument('-p', '--port',
                        help='The host of the mongo database. Mandatory if not ',
                        type=str,
                        default='27017')

    parser.add_argument('-n', '--night',
                        help='Observing night to work on. If not given consider yesterday.',
                        type=str)

    parser.add_argument('-o', '--output',
                        help='Output name.',
                        type=str)

    parser.add_argument('--verbose', '-v', action='count')

    args = parser.parse_args(argv[1:])

    logging.basicConfig(format='%(levelname)s:%(asctime)s::%(message)s',
                        level=args.verbose)

    log = logging.getLogger(__name__)

    client = MongoClient('%s:%s' % (args.host,
                                    args.port))

    # Query all images for that specific night
    today = time.localtime()
    night = '%i%i%i' % (today.tm_year,
                        today.tm_mon,
                        today.tm_mday)
    if args.night is not None:
        night = args.night

    log.info('Querying images from %s.' % night)
    import re
    regx = re.compile('.*', re.IGNORECASE)
    cursor = client.images.fits_keywords.find({'_night': night,
                                               'OBJECT': regx}).sort([('OBJECT', 1),
                                                                      ('FILENAME', 1), ])

    log.debug('Found %i entries.' % cursor.count())

    # object_list = np.unique(np.array([entry['OBJECT'] for entry in cursor]))
    infos = np.array([(entry['OBJECT'], entry['FILENAME'], entry['FILTER'],) for entry in cursor]).T
    object_list = np.unique(infos[0])

    targets = {'night' : night,
        'objects': {}}

    log.debug('Found %i objects.' % len(object_list))
    for obj in object_list:
        log.debug('%s' % obj)
        targets['objects'][str(obj)] = {'type': 'UNKNOWN',
                                        'night': night}
        mask = infos[0] == obj
        ufilters = np.unique(infos[2][mask])
        for filter in ufilters:
            targets['objects'][str(obj)][str(filter)] = {}
            mmask = np.bitwise_and(mask,
                                   infos[2] == filter)
            raw = [str(raw) for raw in list(infos[1][mmask])]
            targets['objects'][str(obj)][str(filter)]['raw files'] = raw


    log.debug('Querying for bias images.')
    # Querying for calibrations in this night. If there are calibration missing, you can search it with other script
    cursor = client.images.fits_keywords.find({'_night': night,
                                               'IMAGETYP': 'ZERO'}).sort([('OBJECT', 1),
                                                                          ('FILENAME', 1), ])

    targets['calibrations'] = {'bias': {},
                               'sky-flat': {},
                               }

    log.debug('Found %i bias frames this night.' % cursor.count())

    targets['calibrations']['bias']['night'] = night
    targets['calibrations']['bias']['raw files'] = [str(frame['FILENAME']) for frame in cursor]

    log.debug('Querying for sky-flat images.')
    cursor = client.images.fits_keywords.find({'_night': night,
                                               'IMAGETYP': 'sky-flat'}).sort([('FILENAME', 1), ])
    log.debug('Found %i sky-flat frames this night.' % cursor.count())

    infos = np.array([(str(entry['FILENAME']),
                       str(entry['FILTER'])) for entry in cursor]).T

    ufilters = [str(flt) for flt in np.unique(infos[1])]

    # filters = [str(flt) for flt in infos[1]]

    targets['calibrations']['sky-flat']['filters'] = ufilters

    for filter in ufilters:
        targets['calibrations']['sky-flat'][filter] = { 'night' : night }
        mask = infos[1] == filter
        raw = [str(ff) for ff in infos[0][mask]]
        targets['calibrations']['sky-flat'][filter]['raw files'] = raw

    log.info('Saving selection to %s' % args.output)
    with open(args.output, 'w') as fp:
        yaml.dump(targets, fp, default_flow_style = False)

    return 0

if(__name__ == "__main__"):
    args = sys.argv
    print(args)
    main(args)
