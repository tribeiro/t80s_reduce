
import os
import logging
from pymongo import MongoClient
import time
import argparse
import re
from numpy import array, unique, bitwise_and
import yaml
import glob
from astropy.io import fits

class SortData:

    def __init__(self, argv):
        parser = argparse.ArgumentParser(description='Sort files from an observing night or dataset and generate '
                                                     'configuration scripts.')

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

        parser.add_argument('-i', '--input',
                            help='Folder with images.',
                            type=str)

        parser.add_argument('--generate-move-script',
                            help='When inspecting a folder with data do reduce, generate a script to move the data into '
                                 'appropriate folders for processing.',
                            action='store_true')

        parser.add_argument('--generate-copy-script',
                            help='Same as "--generate-move-script" but use copy instead. This will duplicate the data!',
                            action='store_true')

        parser.add_argument('--work-path',
                            help='Choose a working path to store/move/copy the data. Default is either input or local '
                                 'directory, depending on running condition.',
                            type=str,
                            default=None)

        parser.add_argument('--verbose', '-v', action='count')

        self.args = parser.parse_args(argv[1:])

        logging.basicConfig(format='%(levelname)s:%(asctime)s::%(message)s',
                            level=self.args.verbose)

        self.log = logging.getLogger(__name__)

        self.client = MongoClient('%s:%s' % (self.args.host,
                                        self.args.port))
        today = time.localtime()

        self.night = self.args.night if self.args.night is not None else '%i%i%i' % (today.tm_year,
                                                                                     today.tm_mon,
                                                                                     today.tm_mday)

        if self.args.night is not None:
            self.night = self.args.night

        self.targets = {'night': self.night,
                        'objects': {},
                        'calibrations' : {'bias': {},
                                          'sky-flat': {}}
                        }

        self.move_cmd = ''

    def cursor(self, query, sort):
        search = {'_night': self.night}
        for key, value in query.items():
            search[key] = value

        cursor = self.client.images.fits_keywords.find(search).sort(sort)
        return cursor

    def queryImage(self):
        regx = re.compile('.*', re.IGNORECASE)

        cursor = self.cursor({'OBJECT':regx}, [('OBJECT', 1),
        ('FILENAME', 1), ])

        self.log.debug('Found %i entries.' % cursor.count())

        infos = array([(entry['OBJECT'], entry['FILENAME'], entry['FILTER'],)
                       for entry in cursor]).T

        object_list = unique(infos[0])

        self.log.debug('Found %i objects.' % len(object_list))
        for obj in object_list:
            self.log.debug('%s' % obj)
            self.targets['objects'][str(obj)] = {'type': 'UNKNOWN',
                                            'night': self.night}
            mask = infos[0] == obj
            ufilters = unique(infos[2][mask])
            for filter in ufilters:
                self.targets['objects'][str(obj)][str(filter)] = {}
                mmask = bitwise_and(mask,
                                       infos[2] == filter)
                raw = [str(raw) for raw in list(infos[1][mmask])]
                self.targets['objects'][str(obj)][str(filter)]['raw files'] = raw
        return self

    def queryNight(self):

        regx = re.compile('.*', re.IGNORECASE)

        cursor = self.cursor({'_night': self.night, 'OBJECT': regx},
                             [('OBJECT', 1),('FILENAME', 1), ]
                             )

        self.log.debug('Found %i entries.' % cursor.count())

        infos = array([(entry['OBJECT'], entry['FILENAME'], entry['FILTER'],)
                       for entry in cursor]).T

        object_list = unique(infos[0])

        self.log.debug('Found %i objects.' % len(object_list))
        for obj in object_list:
            self.log.debug('%s' % obj)
            self.targets['objects'][str(obj)] = {'type': 'UNKNOWN',
                                            'night': self.night}
            mask = infos[0] == obj
            ufilters = unique(infos[2][mask])
            for filter in ufilters:
                self.targets['objects'][str(obj)][str(filter)] = {}
                mmask = bitwise_and(mask,
                                       infos[2] == filter)
                raw = [str(raw) for raw in list(infos[1][mmask])]
                self.targets['objects'][str(obj)][str(filter)]['raw files'] = raw
        return self

    def queryBias(self):
        self.log.debug('Querying for bias images.')
        # Querying for calibrations in this night. If there are calibration missing, you can search it with other script
        cursor = self.cursor({'IMAGETYP': 'ZERO'}, [('OBJECT', 1),
                                   ('FILENAME', 1), ])

        self.targets['calibrations'] = {'bias': {},
                                   'sky-flat': {},
                                   }

        self.log.debug('Found %i bias frames this night.' % cursor.count())

        self.targets['calibrations']['bias']['night'] = self.night

        self.targets['calibrations']['bias']['raw files'] = [
                                    str(frame['FILENAME']) for frame in cursor]
        return self

    def querySkyFlat(self):
        self.log.debug('Querying for sky-flat images.')
        cursor = self.cursor({'IMAGETYP': 'sky-flat'}, [('FILENAME', 1), ])
        self.log.debug('Found %i sky-flat frames this night.' % cursor.count())

        infos = array([(str(entry['FILENAME']),
                           str(entry['FILTER'])) for entry in cursor]).T

        ufilters = [str(flt) for flt in unique(infos[1])]

        self.targets['calibrations']['sky-flat']['filters'] = ufilters

        for filter in ufilters:
            self.targets['calibrations']['sky-flat'][filter] = { 'night' :
                                                                self.night }
            mask = infos[1] == filter
            raw = [str(ff) for ff in infos[0][mask]]
            self.targets['calibrations']['sky-flat'][filter]['raw files'] = raw
        return self

    def saveFile(self):
        self.log.info('Saving selection to %s' % self.args.output)
        with open(self.args.output, 'w') as fp:
            yaml.dump(self.targets, fp, default_flow_style = False)

        if len(self.move_cmd) > 0 and (self.args.generate_copy_script or self.args.generate_move_script):
            self.log.info('Saving %s script to %s' % ('move' if self.args.generate_move_script else 'copy',
                                                      self.args.output+'_move.sh'))
            with open(self.args.output+'_move.sh', 'w') as fp:
                fp.write(self.move_cmd)


    def queryFromLocalImages(self):

        if(self.args.input is None):
            return self

        files = glob.glob(self.args.input + "/*.fits")

        if(len(files) == 0):
            self.log.debug("No objects found in %s" % self.args.input)
            raise NameError("No objects in the Directory %s")

        self.log.debug('Found %i objects,' % len(files))
        self.log.debug('in the directory %s objects.' % self.args.input)
        i = 0
        self.targets = {'night' : self.night,
            'objects': {},
                        'path' : self.args.input}
        teste = []
        raws = {}

        mv_dir = ''
        cmd = 'mv' if self.args.generate_move_script else 'cp'

        for obj in files:
            hdulist = fits.open(obj)
            header = hdulist[0].header
            raw = obj.split("/")[-1]

            if('OBJECT' in header.keys()):
                self.log.debug('%s' % header['OBJECT'])
                if(header['OBJECT'] not in teste):
                    teste.append(header['OBJECT'])
                    self.targets['objects'][header['OBJECT']] = {'type': 'UNKNOWN',
                                                'night': self.night}
                    self.targets['objects'][header['OBJECT']][
                        header["FILTER"]] = { 'raw files': [raw]}
                else:
                    if(header["FILTER"] in self.targets['objects'][header['OBJECT']].keys()):
                        self.targets['objects'][header['OBJECT']][
                        header["FILTER"]]['raw files'].append(raw)
                    else:
                        self.targets['objects'][header['OBJECT']][
                        header["FILTER"]] = { 'raw files': [raw]}
            else:
                continue

            if self.args.work_path is not None:
                new_mv_dir = os.path.join(self.args.work_path,
                                          self.night,
                                          header['OBJECT'].replace(' ', '_'),
                                          header['FILTER'],
                                          'raw_files/')
            else:
                new_mv_dir = os.path.join(self.args.input,
                                          self.night,
                                          header['OBJECT'].replace(' ', '_'),
                                          header['FILTER'],
                                          'raw_files/')
            if not os.path.exists(new_mv_dir) and new_mv_dir != mv_dir:
                self.move_cmd += 'mkdir -p %s\n' % new_mv_dir
                mv_dir = new_mv_dir

            if self.args.generate_move_script or self.args.generate_copy_script:
                self.move_cmd += '%s -v %s %s\n' % (cmd,
                                                    obj,
                                                    mv_dir)

        return self

