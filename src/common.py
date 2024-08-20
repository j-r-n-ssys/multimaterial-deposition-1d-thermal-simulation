"""Container for common functions."""

import logging as lg
import time

from datetime import datetime
from pathlib import Path

LOGGER_STR_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

LOGGER_DATETIME_FORMAT = '%Y-%m-%d %H:%M:%S.%f'


class Formatter(lg.Formatter):
    """Extension of the logging module formatter with prefered datetime formatting. """

    def formatTime(self, record, datefmt=None):  #pylint disable=invalid-name
        """Datetime formatter."""
        ct = self.converter(record.created)
        if datefmt is not None:
            # support %z and %f in datefmt (struct_time doesn't carry ms or tz)
            datefmt = datefmt.replace('%f', str(datetime.fromtimestamp(record.created).microsecond).zfill(6))
            datefmt = datefmt.replace('%z', time.strftime('%z'))
            s = time.strftime(datefmt, ct)
        else:
            t = time.strftime('%Y-%m-%d %H:%M:%S', ct)
        return s


logger = lg.getLogger('adhesion_modeling')

logger.setLevel(lg.DEBUG)

ch = lg.StreamHandler()

# create the console handler
ch.setLevel(lg.DEBUG)

formatter = Formatter(
    fmt=LOGGER_STR_FORMAT,
    datefmt=LOGGER_DATETIME_FORMAT,
)

# Add formatter to console handler
ch.setFormatter(formatter)

# add console handler to logger
logger.addHandler(ch)

if __name__ == '__main__':
    print(f'\nModule <{Path(__file__).name}> is not intended to be run as standalone module.')

    logger.debug('a')
