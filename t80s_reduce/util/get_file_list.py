import os
import logging
from t80s_reduce.core.constants import *

log = logging.getLogger(__name__)


def get_flat_list(config, args, get_file_type='raw files', write_file_type=None):
    img_list = []
    # get list of flats
    for filter in config['calibrations']['sky-flat']['filters']:
        path = os.path.join(config['path'],
                            config['calibrations']['sky-flat'][filter]['night'],
                            'flat',
                            filter)
        if (write_file_type is not None) and (args.overwrite) and \
                (write_file_type in config['calibrations']['sky-flat'][filter]):
            log.warning('Running in overwrite mode. Cleaning existing %s list.' % write_file_type)
            config['calibrations']['sky-flat'][filter][write_file_type] = []
        elif write_file_type is not None and write_file_type in config['calibrations']['sky-flat'][filter]:
            log.warning('Flat field frames on %s already processed but running in non-overwrite mode. '
                        'Skipping.' % filter)
            continue
        else:
            config['calibrations']['sky-flat'][filter][write_file_type] = []

        for flat in config['calibrations']['sky-flat'][filter][get_file_type]:
            img_list.append((os.path.join(path, flat),
                             os.path.join(path, 'b_' + flat),
                             ('flat', filter)))
        return img_list

def get_target_list(config, args, get_file_type='raw files', write_file_type=None):
    img_list = []

    for object in config['objects']:
        for filter in FILTERS:
            if filter not in config['objects'][object]:
                continue
            if (write_file_type is not None) and (args.overwrite) and \
                    (write_file_type in config['objects'][object][filter]):
                log.warning('Running in overwrite mode. Cleaning existing %s list.' % write_file_type)
                config['objects'][object][filter][write_file_type] = []
            elif (write_file_type is not None) and (write_file_type in config['objects'][object][filter]):
                log.warning('%s frames in %s already processed but running in non-overwrite mode. '
                            'Skipping.' % (object, filter))
                continue
            else:
                config['objects'][object][filter][write_file_type] = []

            path = os.path.join(config['path'],
                                config['objects'][object]['night'],
                                object.replace(' ', '_'),
                                filter)
            for raw in config['objects'][object][filter][get_file_type]:
                img_list.append((os.path.join(path, raw),
                                 os.path.join(path, 'b_' + raw),
                                 ('target', object, filter)))
    return img_list
