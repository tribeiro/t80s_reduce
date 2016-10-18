import sys, os
import numpy as np
import logging
from datetime import datetime
import yaml
from t80s_reduce.core.constants import *

def main(argv):
    import argparse

    parser = argparse.ArgumentParser(description='Given a configuration file, check the calibration frames are'
                                                 'ok. If not, it will write specific keywords that can be used'
                                                 'to select calibrations using another script.')

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

    # Open configuration file
    log.debug('Reading in configuration file: %s' % args.config)
    with open(args.config,'r') as fp:
        config = yaml.load(fp)

    # Check that 'calibrations' check exists
    if 'calibrations' not in config:
        # print config
        log.warning('No section "calibrations" in configuration file. Creating empty one.')
        config['calibrations'] = {'sky-flat' : {},
                                  'actions': {'comments':
                                                  ['%s - Creating section "actions".' % datetime.today().__str__()],
                                              'flat': {}}}

    else:
        # if calibrations exists check that 'actions' section is present and create one, otherwise
        if 'actions' in config['calibrations']:
            log.debug('Reprocessing calibration file.')
        else:
            log.debug('Apparently, first pass on this configuration file. Creating "actions" section.')
            config['calibrations']['actions'] = {'comments' :
                                                     ['%s - Creating section "actions".' % datetime.today().__str__()],
                                                 }
            config['calibrations']['actions']['bias'] = {'closest' : False, 'avoid-nights' : [ ], 'ok' : False}
            config['calibrations']['actions']['flat'] = {}


    # Check for bias
    if 'bias' not in config['calibrations']:
        log.warning('No bias frames in calibrations. Will mark to look for closest in time.')
        config['calibrations']['actions']['bias'] = {'closest' : True, 'avoid-nights' : [ ], 'ok' : False}
        config['calibrations']['actions']['comments'].append('%s - No bias frames in calibrations. '
                                                             'Look for closest in time.' %
                                                             datetime.today().__str__())
    elif not config['calibrations']['actions']['bias']['ok']:
        # check that there is "raw files" in bias and how many
        if "raw files" not in config['calibrations']['bias']:
            log.info('Bias section defined in configuration file but no raw file list defined.')
            # what to do?
        else:
            # compute number of bias frames.
            nbias = len(config['calibrations']['bias']['raw files'])
            if nbias < NBIAS:
                msg = 'Insuficient number of bias frames. Found %i, expect %i.' % (nbias, NBIAS)
                log.warning(msg)
                config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))
                config['calibrations']['actions']['bias']['avoid-nights'].append(
                    config['calibrations']['bias']['night'])
            else:
                msg = 'Bias ok'
                log.info(msg)
                config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))
                config['calibrations']['actions']['bias']['ok'] = True
    else:
        log.info('Bias ok.')

    # Check for flat fields

    # Build a list of all required filters
    required_filters = []
    for objects in config['objects']:
        obj_filters =  config['objects'][objects].keys()
        for filter in obj_filters:
            if filter not in required_filters and filter in FILTERS:
                required_filters.append(filter)
    log.info('%i required filters: %s' % (len(required_filters),
                                          required_filters))

    if 'filters' not in config['calibrations']['sky-flat']:
        config['calibrations']['sky-flat']['filters'] = required_filters

    nmissing = 0
    nok = 0
    nfault = 0
    for filter in required_filters:
        if filter in config['calibrations']['sky-flat']:
            if filter not in config['calibrations']['actions']['flat']:
                config['calibrations']['actions']['flat'][filter] = {}
                config['calibrations']['actions']['flat'][filter]['closest'] = False
                config['calibrations']['actions']['flat'][filter]['ok'] = False
                config['calibrations']['actions']['flat'][filter]['avoid-nights'] = []
                config['calibrations']['actions']['flat'][filter]['start-night'] = config['night']
            elif config['calibrations']['actions']['flat'][filter]['ok']:
                log.debug('Filter %s already checked.' % filter)
                nok += 1
                continue

            msg = 'Filter %s in calibrations' % filter
            log.debug(msg)
            config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))

            if len(config['calibrations']['sky-flat'][filter]['raw files']) < NFLATS:
                msg = 'Insufficient number of flats in filter %s (%i/%i)' % (
                    filter,
                    len(config['calibrations']['sky-flat'][filter]['raw files']),
                    NFLATS)
                log.info(msg)
                config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))

                config['calibrations']['actions']['flat'][filter]['closest'] = True
                config['calibrations']['actions']['flat'][filter]['ok'] = False
                config['calibrations']['actions']['flat'][filter]['avoid-nights'].append(
                    config['calibrations']['sky-flat'][filter]['night'])
                config['calibrations']['actions']['flat'][filter]['start-night'] = config['night']

                nfault += 1

            else:
                msg = 'Filter %s has %i flat frames. OK' % (
                    filter,
                    len(config['calibrations']['sky-flat'][filter]['raw files']))
                log.info(msg)
                config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))

                config['calibrations']['actions']['flat'][filter]['ok'] = True
                nok += 1

        else:
            nmissing += 1
            msg = 'Filter %s missing in calibrations' % filter
            log.debug(msg)
            config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))
            if filter in config['calibrations']['actions']['flat']:
                config['calibrations']['actions']['flat'][filter]['closest'] = True
                config['calibrations']['actions']['flat'][filter]['ok'] = False
                config['calibrations']['actions']['flat'][filter]['start-night'] = config['night']
                # config['calibrations']['actions']['flat'][filter]['avoid-nights'].append(
                #     config['calibrations']['sky-flat']['night'])
            else:
                config['calibrations']['actions']['flat'][filter] = {}
                config['calibrations']['actions']['flat'][filter]['closest'] = True
                config['calibrations']['actions']['flat'][filter]['ok'] = False
                config['calibrations']['actions']['flat'][filter]['avoid-nights'] = []
                config['calibrations']['actions']['flat'][filter]['start-night'] = config['night']

    msg = 'Flat-field %i missing, %i present but bad and %i ok.' % (nmissing, nfault, nok)
    log.debug(msg)
    config['calibrations']['actions']['comments'].append('%s - %s' % (datetime.today().__str__(), msg))

    log.info('Saving selection to %s' % args.output)
    with open(args.output, 'w') as fp:
        yaml.dump(config, fp, default_flow_style = False)

    return 0
