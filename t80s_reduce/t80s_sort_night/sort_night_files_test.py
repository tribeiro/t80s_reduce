#!/usr/bin/env python
# -*- Coding: UTF-8 -*-

__data = "27/09/2016"
__autor = "Eduardo S. Pereira"

import unittest
import sys
from sort_night_files import SortNightFiles

class TestSortNightFiles(unittest.TestCase):

    @classmethod
    def setUp(self):

        self.sort_night_file = SortNightFiles(["__main__",
                                    '-o', 'teste.yaml',
                                    '-i', '/home/eduardo/ProjetosEnv/t80/dados'
                                              ])

    def test_cursor(self):
        self.sort_night_file.cursor({'IMAGETYP': 'ZERO'}, [('OBJECT', 1),
        ('FILENAME', 1), ])

    def test_queryFromLocalImages(self):
        self.sort_night_file = SortNightFiles(["__main__",
                                    '-o', 'new.yaml',
                                    '-i', '/home/eduardo/ProjetosEnv/t80/dados'
                                    ])
        self.sort_night_file.queryFromLocalImages().saveFile()

    def test_saveFile(self):
        self.sort_night_file.target = {"test": "teste"}
        self.sort_night_file.saveFile()

    #def test_queryImage(self):
    #    self.sort_night_file.queryImage()

    #def test_querySkyFlat(self):
    #    self.sort_night_file.querySkyFlat()

    #def test_queryBias(self):
    #    self.sort_night_file.queryBias()


if(__name__ == "__main__"):
    unittest.main()
