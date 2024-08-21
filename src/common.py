"""Container for common functions."""

import logging as lg

from datetime import datetime
from pathlib import Path

LOGGER_STR_FORMAT = '%(asctime)s - %(name)s - %(levelname)s - message: %(message)s'

LOGGER_DATETIME_FORMAT = '%Y-%m-%d %H:%M:%S.%f%Z'


class CustomFormatter(lg.Formatter):
    """Logging module formatter with custom time formatter. """

    def formatTime(self, record, datefmt=None):  #pylint disable=invalid-name
        """Customer datetime formatter."""

        dt = datetime.fromtimestamp(record.created)

        s = dt.strftime('%Y-%m-%d %H:%M:%S.%f')[:-3] if datefmt is None else dt.strftime(datefmt)

        return s


logger = lg.getLogger('adhesion_modeling')

logger.setLevel(lg.DEBUG)

ch = lg.StreamHandler()

# create the console handler
ch.setLevel(lg.DEBUG)

formatter = CustomFormatter(
    fmt=LOGGER_STR_FORMAT,
    datefmt=LOGGER_DATETIME_FORMAT,
)

# Add formatter to console handler
ch.setFormatter(formatter)

# add console handler to logger
logger.addHandler(ch)

if __name__ == '__main__':
    print(f'Module {Path(__file__).name} is not intended to be run as standalone module.')

    logger.debug('a')
