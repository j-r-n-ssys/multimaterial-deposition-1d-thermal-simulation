"""Collection of class definitions and functions for time-temperature superposition."""

import logging as lg
import numpy as np

from os.path import basename as get_module_fname

NUMERICAL_TYPES = (int, float)

CELSIUS_TO_KELVIN_OFFSET = 273.15

UNIVERSAL_GAS_CONSTANT = 8.31446261815324  # J/K-mol


class WilliamLandelFerryModel():
    """This object represents a William-Landel-Ferry time-temperature superposition model."""

    def __init__(self, c_1, c_2, t_ref: float, t_glass: float | None) -> None:
        """_summary_

        Args:
            c_1 (_type_): _description_
            c_2 (float): WLF model 
            t_ref (float): Time-temperature superposition reference temperature [degC]. 
            t_glass (float | None): Glass transition temperature [degC]
        """

        if not isinstance(c_1, NUMERICAL_TYPES):
            raise TypeError(f'WLF Model empirical constant C1 must be a numerical type, not {type(c_1)}.')
        elif not isinstance(c_2, NUMERICAL_TYPES):
            raise TypeError(f'WLF Model empirical constant C2 must be a numerical type, not {type(c_2)}.')
        elif c_1 <= 0 or c_2 <= 0:
            raise ValueError('WLF Model empirical constants C1 and C2 must be greater than zero.')

        if not isinstance(t_ref, NUMERICAL_TYPES):
            raise TypeError(f'Temperature T_ref mmust be a numerical type, not {type(t_ref)}.')
        elif t_ref < -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Reference temperature T_ref must be greater than absolute zero.')

        if not (t_glass is None or isinstance(t_glass, NUMERICAL_TYPES)):
            raise TypeError(f'Temperature T_ref mmust be a numerical type, not {type(t_glass)}.')
        elif t_glass is not None and t_glass < -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Glass transition temperature T_ref must be greater than absolute zero.')

        self.c_1 = c_1
        self.c_2 = c_2
        self.t_r = t_ref + CELSIUS_TO_KELVIN_OFFSET
        self.t_g = t_glass + CELSIUS_TO_KELVIN_OFFSET if t_glass is not None else None

    def get_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """_summary_

        Args:
            temp (float  |  np.ndarray): Temperature [degC]

        Returns:
            float  |  np.ndarray: _description_
        """

        if self.t_g is not None:
            if np.max(temp) > self.t_g + 100:
                lg.warning('Temperature exceeds WLF model validity upper bound.')
            elif np.min(temp) < self.t_g:
                lg.warning('Temperature exceeds WLF model validity lower bound.')

        temp = temp + CELSIUS_TO_KELVIN_OFFSET

        if isinstance(temp, NUMERICAL_TYPES):
            return 10**(-self.c_1 * (temp - self.t_r) / (self.c_2 + (temp - self.t_r)))
        elif isinstance(temp, np.ndarray):
            return 10**(-self.c_1 * np.divide((temp - self.t_r), (self.c_2 + (temp - self.t_r))))
        else:
            raise TypeError('Temperature argument must be a numerical type or numerical array.')


class ArrheniusModel():
    """This object represents a Arrhenius time-temperature superposition model."""

    def __init__(self, e_a: float, t_ref: float) -> None:

        self.e_a = e_a
        self.t_r = t_ref

    def get_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """_summary_

        Args:
            temp (float  |  np.ndarray): Temperature [degC]

        Returns:
            float  |  np.ndarray: _description_
        """

        return 10**((-self.e_a / 2.303 * UNIVERSAL_GAS_CONSTANT) * (1 / temp - 1 / self.t_r))


if __name__ == '__main__':
    if False:
        lg.warning('Module %s is not intended to be run as standalone module.', get_module_fname(__file__))

    print(' ')
    print(' ')

    wlf = WilliamLandelFerryModel(5.78, 182, 200, None)

    print(wlf.get_shift_factor(205))
