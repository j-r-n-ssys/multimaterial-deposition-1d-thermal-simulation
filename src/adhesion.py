"""Collection of class definitions and functions for time-temperature superposition."""

import logging as lg
import numpy as np

from os.path import basename as get_module_fname

NUMERICAL_TYPES = (int, float)

CELSIUS_TO_KELVIN_OFFSET = 273.15


class WilliamLandelFerryModel():
    """This model represents the William-Landel-Ferry time-temperature superposition model. 
    """

    def __init__(self, c_1, c_2, t_ref: float, t_glass: float) -> None:
        """_summary_

        Args:
            c_1 (_type_): _description_
            c_2 (float): WLF model 
            t_ref (float): Time-temperature superposition reference temperature [degC]. 
            t_glass (float): Glass transition temperature [degC]
        """

        if not isinstance(c_1, NUMERICAL_TYPES) or not isinstance(c_2, NUMERICAL_TYPES):
            raise TypeError('WLF Model empirical constant C1 and C2 must be a numerical type.')
        elif not isinstance(t_ref, NUMERICAL_TYPES) or not isinstance(t_glass, NUMERICAL_TYPES):
            raise TypeError('Temperatures T_ref and T_glass mmust be a numerical type.')
        elif c_1 <= 0 or c_2 <= 0:
            raise ValueError('WLF Model empirical constants C1 and C2 must be greater than zero.')
        elif t_ref < -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Temperatures T_ref and T_glass must be greater than absolute zero.')

        self.c_1 = c_1
        self.c_2 = c_2
        self.t_r = t_ref + CELSIUS_TO_KELVIN_OFFSET
        self.t_g = t_glass + CELSIUS_TO_KELVIN_OFFSET

    def get_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """_summary_

        Args:
            temp (float  |  np.ndarray): Temperature [degC]

        Returns:
            float  |  np.ndarray: _description_
        """

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


# wlf_model = lambda C_1, C_2, T_ref, T: -C_1 * (T - T_ref) / (C_2 + (T - T_ref))

if __name__ == '__main__':
    lg.warning('Module %s is not intended to be run as standalone module.', get_module_fname(__file__))
