"""Module to contain the common material class."""

import logging as lg
import math
import numpy as np

from abc import abstractmethod
from os.path import basename as get_module_fname
from types import NoneType

NUMERICAL_TYPES = (int, float)

CELSIUS_TO_KELVIN_OFFSET = 273.15

UNIVERSAL_GAS_CONSTANT = 8.31446261815324  # J/K-mol


class AbstractAdhesionModel():

    @abstractmethod
    def get_shift_factor(self) -> (float | np.ndarray):  #pylint: disable=unused-argument
        ...


class WilliamLandelFerryModel(AbstractAdhesionModel):
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

        self._c_1 = float(c_1)
        self._c_2 = float(c_2)
        self._t_r = float(t_ref) + CELSIUS_TO_KELVIN_OFFSET
        self._t_g = float(t_glass) + CELSIUS_TO_KELVIN_OFFSET if t_glass is not None else None

    def get_shift_factor(self, temp: (int | float | np.ndarray)) -> (float | np.ndarray):
        """Get time-temperature position shift factor at temperature.

        Args:
            temp (float  |  np.ndarray): Temperature [degC].

        Returns:
            float  |  np.ndarray: Time-temperature superposition shift factor.
        """

        if not isinstance(temp, (int, float, np.ndarray)):
            raise TypeError(f'temp must be a numerical type or array, not a {type(temp)}')

        if self._t_g is not None:
            if np.max(temp) > self._t_g + 100:
                lg.warning('Temperature exceeds WLF model validity upper bound.')
            elif np.min(temp) < self._t_g:
                lg.warning('Temperature exceeds WLF model validity lower bound.')

        temp = temp + CELSIUS_TO_KELVIN_OFFSET

        if isinstance(temp, NUMERICAL_TYPES):
            return 10**(-self._c_1 * (temp - self._t_r) / (self._c_2 + (temp - self._t_r)))
        elif isinstance(temp, np.ndarray):
            return 10**(-self._c_1 * np.divide((temp - self._t_r), (self._c_2 + (temp - self._t_r))))
        else:
            raise TypeError('Temperature argument must be a numerical type or numerical array.')


class ArrheniusModel(AbstractAdhesionModel):
    """This object represents a Arrhenius time-temperature superposition model."""

    def __init__(self, e_a: float, t_ref: float) -> None:
        """Instance init.

        Args:
            e_a (float): Activation energy.
            t_ref (float): Reference temperature [degC].
        """

        if not isinstance(e_a, NUMERICAL_TYPES):
            raise TypeError(f'e_a must be a numerical type, not {type(e_a)}')
        elif e_a < 0:
            raise ValueError('e_a must be a positive value.')

        if not isinstance(t_ref, NUMERICAL_TYPES):
            raise TypeError(f't_ref must be a numerical type, not a {type(t_ref)}')
        elif t_ref + CELSIUS_TO_KELVIN_OFFSET < 0:
            raise ValueError('t_ref must be greater than absolute zero.')

        self._e_a = e_a
        self._t_r = t_ref + CELSIUS_TO_KELVIN_OFFSET

    def get_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """_summary_

        Args:
            temp (float  |  np.ndarray): Temperature [degC].

        Returns:
            float  |  np.ndarray: 
        """

        

        f = (-self._e_a / 2.303 * UNIVERSAL_GAS_CONSTANT)

        if isinstance(temp, float):
            return 10**(f * (1 / (temp + CELSIUS_TO_KELVIN_OFFSET) - 1 / self._t_r))
        elif isinstance(temp, np.ndarray):
            return np.power(10, f * (np.divide(1, temp + CELSIUS_TO_KELVIN_OFFSET) - (1 / self._t_r)))
        else:
            raise TypeError('Temperature argument must be a numerical type or numerical array.')


ADHESION_MODELS = (WilliamLandelFerryModel, ArrheniusModel)


class Material():
    """This object represents a material."""

    def __init__(self, name: str, density: float, thermal_conductivity: float, specific_heat_capacity: float,
                 glass_transition: float, melt_transition: float, extrusion_temp: float, emissivity: float):
        """Init."""

        if not isinstance(name, str):
            raise TypeError('Name must be a string.')

        self.long_name = str(name)

        if not isinstance(density, NUMERICAL_TYPES):
            raise TypeError('Density must be type float.')

        if density < 0.0:
            raise ValueError('Density must be greater than zero.')

        self.density = float(density)  # [kg/m3]

        if not isinstance(thermal_conductivity, NUMERICAL_TYPES):
            raise TypeError('Thermal conductivity must be type float.')

        if thermal_conductivity < 0.0:
            raise ValueError('Thermal conducitivty must be greater than zero.')

        self.thermal_conductivity = float(thermal_conductivity)  # [W/m-K]

        if not isinstance(specific_heat_capacity, NUMERICAL_TYPES):
            raise TypeError('Specific heat capacity must be type float.')

        if specific_heat_capacity < 0.0:
            raise ValueError('Specific heater capacity must be greater than zero.')

        self.specific_heat_capacity = float(specific_heat_capacity)

        if not isinstance(glass_transition, NUMERICAL_TYPES):
            raise TypeError('Glass transition temperature must be type float.')

        self.glass_transition = float(glass_transition)  # [K]

        if not isinstance(melt_transition, (NoneType, int, float)):
            raise TypeError('Melt temperature must be type float.')

        self.melt_transition = float(melt_transition) if not isinstance(melt_transition, NoneType) else None

        if not isinstance(extrusion_temp, NUMERICAL_TYPES):
            raise TypeError('Extrusion temperature must be type float.')

        self.extrusion_temp = float(extrusion_temp)  # [K]

        if not isinstance(emissivity, (NoneType, int, float)):
            raise TypeError('Emissivity must be type float.')

        if not isinstance(emissivity, NoneType) and (emissivity < 0.0 or emissivity > 1.0):
            raise ValueError('Emissivity must be between 0 and 1.')

        self.emmisivity = float(emissivity) if not isinstance(emissivity, NoneType) else None  # [0.0-1.0]

        self._adhesion_model: WilliamLandelFerryModel | ArrheniusModel | None = None

    @property
    def volumetric_heat_capacity(self) -> float:
        """Volumetric heat capacity [J/cu.m-K]."""
        return self.density * self.specific_heat_capacity  # [J/m3-K]

    @property
    def thermal_effusivity(self) -> float:
        """Thermal effusivity []. """
        return math.sqrt(self.thermal_conductivity * self.volumetric_heat_capacity)

    @property
    def thermal_diffusivity(self) -> float:
        """Thermal diffusivity [sq.m/s]."""
        return self.thermal_conductivity / self.volumetric_heat_capacity  # [m2/s]

    @property
    def adhesion_model(self) -> (WilliamLandelFerryModel | ArrheniusModel):
        """Adhesion model instance."""
        return self._adhesion_model

    @adhesion_model.setter
    def adhesion_model(self, adhesion_model: WilliamLandelFerryModel | ArrheniusModel):
        if not isinstance(adhesion_model, ADHESION_MODELS):
            raise type(f'adhesion_model must be an adhesion model instance, not {type(adhesion_model)}')

        self._adhesion_model = adhesion_model

    @property
    def has_adhesion_model(self) -> bool:
        return self._adhesion_model is not None


def calculate_bond(material: Material, time_arr: np.ndarray, temp_arr: np.ndarray) -> float:
    """Calculate the time-temperature-strength integral.

    Args:
        material (Material): Material.
        time_arr (np.ndarray): Array of time values.
        temp_arr (np.ndarray): Array of temperatures.
        dt (float): Time step.

    Returns:
        float: Relative bond strength.
    """

    if not isinstance(material, Material):
        raise TypeError(f'material must be a Material instance, not {type(material)}.')
    elif not material.has_adhesion_model:
        raise ValueError('Material instance does not have an adhesion model.')

    if not isinstance(time_arr, np.ndarray):
        raise TypeError(f'time_arr must be an array, not {type(time_arr)}.')
    elif len(time_arr) < 2:
        raise IndexError('time_arr must have at least two values.')

    if not isinstance(temp_arr, np.ndarray):
        raise TypeError(f'temp_arr must be an array, not {type(temp_arr)}.')
    elif len(temp_arr) != len(time_arr):
        raise ValueError('temp_arr must have the same number of elements as time_arr.')

    shift_factors = material.adhesion_model.get_shift_factor(temp_arr)

    z = np.divide(1.0, 2.0 * np.pi * shift_factors)

    z[temp_arr < material.glass_transition] = 0

    res = abs(np.sum(np.trapz(y=z, x=time_arr)))**0.25

    if res > 1:
        lg.info('Estimated relative weld strength Adhesion ratio exceeds 1 ( R = %s).', res)

    return min(res, 1)


# https://thermtest.com/thermal-resources/materials-database

ABS = Material('ABS', 1040, 0.209, 1506, 105, None, 260, None)

QSR = Material('QSR', 1180, 0.1, 1950, 165, None, 295, None)

QSR.adhesion_model = WilliamLandelFerryModel(5.78, 182, 200, None)

F375M = Material('F375M Sinterable', 6174.2, 10.6, 942.8, -30.0, 170, 235, None)

F375M.adhesion_model = ArrheniusModel(5.78, 200.0)

if __name__ == '__main__':
    lg.warning('Module %s is not intended to be run as standalone module.', get_module_fname(__file__))
