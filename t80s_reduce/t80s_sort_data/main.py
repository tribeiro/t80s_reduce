# -*- Coding; UTF-0 -*-

__data = '27/09/2016'
__autor = 'Eduardo S. Pereira'

"""
Class to generate the download files given from a query in the astronomical database
"""

from t80s_reduce.util.sortdata import SortData

def main(argv):
    SortData(argv).queryFromLocalImages().saveFile()
