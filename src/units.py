"""Container for unit conversions."""

from pathlib import Path


def inch_to_millimeter(f: float, n: int = 1) -> float:
    """Convert an argument `f` in inches to millimeters. This function accepts
    an optional argument `n`, which specifies the power. For example, a value of 
    `n = 1` converts the argument from inches to millimeters. A value of `n = 3` 
    converts the argument from cubic inches to cubic millimeters. This function
    can also be used to do the reverse conversion, as, for example, `n = -1` 
    converts the argument from millimeteres to inches. """

    if not isinstance(f, (float, int)):
        raise TypeError('Argument f must be of type int or float).')
    elif not isinstance(n, int):
        raise TypeError('Argument n must be type int.')
    elif f == 0:
        raise ValueError('n = 0 has no effect.')

    return float(f) * (25.4**n)


if __name__ == '__main__':
    print(f'\nModule <{Path(__file__).name}> is not intended to be run as standalone module.')
