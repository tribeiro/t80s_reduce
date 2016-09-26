import logging
import logging.handlers
import sys
import os.path

# from chimera.core.constants import (SYSTEM_CONFIG_LOG_NAME,
#                                     SYSTEM_CONFIG_DIRECTORY,
#                                     MANAGER_DEFAULT_HOST,
#                                     MANAGER_DEFAULT_PORT)

# try to use faster (C)StringIO, use slower one if not available
try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

__all__ = ['setConsoleLevel']

try:
    if not os.path.exists(SYSTEM_CONFIG_DIRECTORY):
        os.mkdir(SYSTEM_CONFIG_DIRECTORY)
except Exception:
    pass

root = logging.getLogger("chimera")
root.setLevel(logging.DEBUG)
root.propagate = False

fmt = logging.Formatter(
    fmt='%(asctime)s.%(msecs)d %(origin)s %(levelname)s %(name)s %(filename)s:%(lineno)d %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

flt = ChimeraFilter()

consoleHandler = logging.StreamHandler(sys.stderr)
consoleHandler.setFormatter(fmt)
consoleHandler.setLevel(logging.WARNING)
consoleHandler.addFilter(flt)
root.addHandler(consoleHandler)

def setConsoleLevel(level):
    consoleHandler.setLevel(level)

try:
    fileHandler = logging.handlers.RotatingFileHandler(SYSTEM_CONFIG_LOG_NAME,
                                                       maxBytes=5 *
                                                       1024 * 1024,
                                                       backupCount=10)
    fileHandler.setFormatter(fmt)
    fileHandler.setLevel(logging.DEBUG)
    fileHandler.addFilter(flt)
    root.addHandler(fileHandler)
except Exception, e:
    root.warning("Couldn't start Log System FileHandler (%s)" % e)
