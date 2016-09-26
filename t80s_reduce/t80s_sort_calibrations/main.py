import sys, os
import numpy as np
import logging
from pymongo import MongoClient
import time
from datetime import datetime, timedelta
import yaml
from t80s_reduce.core.constants import *

def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description='Given a configuration file, look into the database for the required'
                                                 'calibration frames.')

    parser.add_argument('--host',
                        help='The host of the mongo database. Mandatory if not ',
                        type=str,
                        default='192.168.20.118')

    parser.add_argument('-p', '--port',
                        help='The host of the mongo database. Mandatory if not ',
                        type=str,
                        default='27017')

    parser.add_argument('-c', '--config',
                        help='Configuration file. The file will be over-written at the end with the result from the '
                             'checks.',
                        type=str)

    parser.add_argument('-o', '--output',
                        help='Output name. If not given overwrite configuration file.',
                        type=str)

    parser.add_argument('--verbose', '-v', action='count')

    args = parser.parse_args(argv[1:])

    logging.basicConfig(format='%(levelname)s:%(asctime)s::%(message)s',
                        level=args.verbose)

    log = logging.getLogger(__name__)

    client = MongoClient('%s:%s' % (args.host,
                                    args.port))

    # Open configuration file
    log.debug('Reading configuration file: %s' % args.config)
    with open(args.config,'r') as fp:
        config = yaml.load(fp)

    # Check that there are actions defined for the calibrations

    if 'actions' not in config['calibrations']:
        log.error('No actions defined on configuration file. Nothing to do.')
        exit(-1)

    if config['calibrations']['actions']['bias']['ok']:
        log.debug('Action bias ok. Skipping.')
    else:
        log.debug('Running action bias')

    for filter in config['calibrations']['actions']['flat']:
        if config['calibrations']['actions']['flat'][filter]['ok']:
            log.debug('Action flat in %s ok. Skipping.' % filter)

        else:
            log.debug('Running action flat in %s' % filter)
            night = config['night']
            if night not in config['calibrations']['actions']['flat'][filter]['avoid-nights']:
                cursor = client.images.skyflats.find({'night': night,
                                                           'filter': filter})
                log.debug('Found %i flats in %s' % (cursor.count(),
                                                    night))

                if cursor.count() < NFLATS:
                    config['calibrations']['actions']['flat'][filter]['avoid-nights'].append(night)

            shift_days = 1
            max_shift = 10
            year,mm,dd = int(night[:4]),int(night[4:6]),int(night[6:8])
            today = datetime(year,mm,dd)

            while True:
                dt = timedelta(days=shift_days)
                next_try = today+dt
                night = '%04i%02i%02i' % (next_try.year,
                                          next_try.month,
                                          next_try.day)

                if night not in config['calibrations']['actions']['flat'][filter]['avoid-nights']:
                    cursor = client.images.skyflats.find({'night': night,
                                               'filter': filter})
                    log.debug('Found %i flats in %s' % (cursor.count(),
                                                        night))
                    if cursor.count() >= NFLATS:
                        if filter not in config['calibrations']['sky-flat']['filters']:
                            config['calibrations']['sky-flat']['filters'].append(filter)
                        config['calibrations']['sky-flat'][filter] = {'night' : night,
                                                                      'raw files' : [str(item['FILENAME'])
                                                                                     for item in cursor]}
                        break
                    else:
                        config['calibrations']['actions']['flat'][filter]['avoid-nights'].append(night)

                else:
                    log.debug('Skipping night %s' % night)

                next_try = today-dt
                night = '%04i%02i%02i' % (next_try.year,
                                          next_try.month,
                                          next_try.day)
                if night not in config['calibrations']['actions']['flat'][filter]['avoid-nights']:
                    cursor = client.images.skyflats.find({'night': night,
                                               'filter': filter})
                    log.debug('Found %i flats in %s' % (cursor.count(),
                                                        night))
                    if cursor.count() >= NFLATS:
                        if filter not in config['calibrations']['sky-flat']['filters']:
                            config['calibrations']['sky-flat']['filters'].append(filter)
                        config['calibrations']['sky-flat'][filter] = {'night' : night,
                                                                      'raw files' : [str(item['FILENAME'])
                                                                                     for item in cursor]}
                        break
                    else:
                        config['calibrations']['actions']['flat'][filter]['avoid-nights'].append(night)

                else:
                    log.debug('Skipping night %s' % night)

                if shift_days > max_shift:
                    log.warning('Maximum number of iterations reached without finding required number of flats.')
                    break
                shift_days += 1

    log.info('Saving selection to %s' % args.output)
    with open(args.output, 'w') as fp:
        yaml.dump(config, fp, default_flow_style = False)

    return 0
