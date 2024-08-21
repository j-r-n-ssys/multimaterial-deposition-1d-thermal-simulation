"""Module to contain the common material class."""

import logging as lg
import math
import numpy as np

from abc import abstractmethod
from pathlib import Path

NUMERICAL_TYPES = (int, float)

CELSIUS_TO_KELVIN_OFFSET = 273.15


class AdhesionModelBase():
    """This abstract class is a base class for adhesion models."""

    @property
    def t_ref(self) -> float:
        """Model activation energy."""
        return self._t_r

    @t_ref.setter
    def t_ref(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Reference temperature must be a numerical type, not a {type(value)}.')
        elif value <= -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Reference temperature must be greater than absolute zero.')

        self._t_r = float(value) + CELSIUS_TO_KELVIN_OFFSET

    @property
    def relaxation_freq_at_crossover(self) -> float:
        """Relaxation time at reference temperature."""
        return self._omega_r

    @relaxation_freq_at_crossover.setter
    def relaxation_freq_at_crossover(self, value: float) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Relaxation frequency at crossover must be a numerical type, not a {type(value)}.')
        elif value <= 0:
            raise ValueError('Relaxation frequency at crossover must be greater than zero.')

        self._omega_r = value

    @abstractmethod
    def calc_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        ...

    def calc_relaxation_time(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """Calculate the relaxation time at temperature

        Args:
            temp (float  |  np.ndarray): Temperature [degC].

        Returns:
            float  |  np.ndarray: Relaxation time(s).
        """

        if not isinstance(temp, (int, float, np.ndarray)):
            raise TypeError(f'Temperature must be a numerical type or a numpy array, not {type(temp)}.')

        if isinstance(temp, NUMERICAL_TYPES):
            if temp < -CELSIUS_TO_KELVIN_OFFSET:
                raise ValueError('Temperature must be greater than absolute zero')

        if isinstance(temp, np.ndarray):
            if temp[temp < -CELSIUS_TO_KELVIN_OFFSET].any():
                raise ValueError('All temperature values must be greater than absolute zero')

        a_t = self.calc_shift_factor(temp)

        return 2 * np.pi * a_t / self._omega_r


class WilliamLandelFerryModel(AdhesionModelBase):
    """This object represents a William-Landel-Ferry time-temperature superposition horizontal shift factor estimation 
    model."""

    def __init__(
        self,
        c_1: float,
        c_2: float,
        t_ref: float,
        omega_r: float,
    ) -> None:
        """Init.

        Args:
            c_1 (float): WLF horizontal shift factor model empircal constant 1.
            c_2 (float): WLF horizontal shift factor model empircal constant 2.
            t_ref (float): WLF horizontal shift factor model reference temperature [degC]. 
            omega_r (float): Relaxation angular frequency at tan(delta) = 0 [rad/s].
        """

        # Type/value execptions are handled by the property.
        self.c_1 = c_1

        # Type/value execptions are handled by the property.
        self.c_2 = c_2

        # Type/value execptions are handled by the property.
        self.t_ref = t_ref

        # Type/value execptions are handled by the property.
        self._omega_r = omega_r

    @property
    def c_1(self) -> float:
        """Model empirical constant C1. """
        return self._c_1

    @c_1.setter
    def c_1(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'C1 must be a numerical type, not a {type(value)}')
        elif value <= 0:
            raise ValueError('C1 must be greater than zero.')

        self._c_1 = float(value)

    @property
    def c_2(self) -> float:
        """Model empirical constant C2."""
        return self._c_2

    @c_2.setter
    def c_2(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'C2 must be a numerical type, not a {type(value)}')
        elif value <= 0:
            raise ValueError('C2 must be greater than zero.')

        self._c_2 = float(value)

    def calc_shift_factor(self, temp: (int | float | np.ndarray)) -> (float | np.ndarray):
        """Calculate the time-temperature position horizontal shift factor at temperature.

        Args:
            temp (float  |  np.ndarray): Temperature [degC].

        Returns:
            float  |  np.ndarray: TTS horizontal shift factor.
        """

        if not isinstance(temp, (int, float, np.ndarray)):
            raise TypeError(f'Temperature must be a numerical type or a numpy array, not {type(temp)}.')

        if isinstance(temp, NUMERICAL_TYPES):
            if temp <= -CELSIUS_TO_KELVIN_OFFSET:
                raise ValueError('Temperature must be greater than absolute zero')

        if isinstance(temp, np.ndarray):
            if temp[temp <= -CELSIUS_TO_KELVIN_OFFSET].any():
                raise ValueError('All temperature values must be greater than absolute zero')

        if not isinstance(temp, (int, float, np.ndarray)):
            raise TypeError(f'temp must be a numerical type or array, not a {type(temp)}')

        temp = temp + CELSIUS_TO_KELVIN_OFFSET

        if isinstance(temp, NUMERICAL_TYPES):
            return 10**(-self._c_1 * (temp - self._t_r) / (self._c_2 + (temp - self._t_r)))
        elif isinstance(temp, np.ndarray):
            return 10**(-self._c_1 * np.divide((temp - self._t_r), (self._c_2 + (temp - self._t_r))))
        else:
            raise TypeError('Temperature argument must be a numerical type or numerical array.')


class ArrheniusModel(AdhesionModelBase):
    """This object represents an Arrhenius time-temperature superposition horizontal shift factor estimation model."""

    UNIVERSAL_GAS_CONSTANT = 8.31446261815324  # J/K-mol

    def __init__(
        self,
        e_a: float,
        t_ref: float,
        omega_r: float,
    ) -> None:
        """Init.

        Args:
            e_a (float): Activation energy.
            t_ref (float): Reference temperature [degC].
            omega_r (float): Relaxation angular frequency at tan(delta) = 0 [rad/s].
        """

        # Type/value execptions are handled by the property.
        self.e_a = e_a

        # Type/value execptions are handled by the property.
        self.t_ref = t_ref

        # Type/value execptions are handled by the property.
        self._omega_r = omega_r

    @property
    def e_a(self) -> float:
        """Model activation energy."""
        return self._e_a

    @e_a.setter
    def e_a(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Activation energy must be a numerical type, not {type(value)}.')
        elif value < 0:
            raise ValueError('Activation energy must be a positive value.')

        self._e_a = float(value)

    def calc_shift_factor(self, temp: (float | np.ndarray)) -> (float | np.ndarray):
        """Calculate the time-temperature position horizontal shift factor at temperature.

        Args:
            temp (float  |  np.ndarray): Temperature [degC].

        Returns:
            float  |  np.ndarray: TTS horizontal shift factor(s)
        """

        if not isinstance(temp, (int, float, np.ndarray)):
            raise TypeError(f'Temperature must be a numerical type or a numpy array, not {type(temp)}.')

        if isinstance(temp, NUMERICAL_TYPES):
            if temp <= -CELSIUS_TO_KELVIN_OFFSET:
                raise ValueError('Temperature must be greater than absolute zero')

        if isinstance(temp, np.ndarray):
            if temp[temp <= -CELSIUS_TO_KELVIN_OFFSET].any():
                raise ValueError('All temperature values must be greater than absolute zero')

        f = -self._e_a / (2.303 * self.UNIVERSAL_GAS_CONSTANT)

        if isinstance(temp, (int, float)):
            return 10**(f * (1 / (temp + CELSIUS_TO_KELVIN_OFFSET) - 1 / self._t_r))
        elif isinstance(temp, np.ndarray):
            return np.power(10, f * (np.divide(1, temp + CELSIUS_TO_KELVIN_OFFSET) - (1 / self._t_r)))
        else:
            raise TypeError(f'Temperature argument must be a numerical type or numerical array, not a {type(temp)}')


ADHESION_MODELS = (WilliamLandelFerryModel, ArrheniusModel)


class Material():
    """This object represents a material."""

    def __init__(
        self,
        name: str,
        density: float,
        thermal_conductivity: float,
        specific_heat_capacity: float,
        glass_transition: float,
        melt_transition: float,
        extrusion_temp: float,
        emissivity: float,
        adhesion_model: ArrheniusModel | WilliamLandelFerryModel | None = None,
    ):
        """Init."""

        self.name = name

        self.density = density

        self.thermal_conductivity = thermal_conductivity

        self.specific_heat_capacity = specific_heat_capacity

        self.glass_transition = glass_transition

        self.melt_transition = melt_transition

        self.extrusion_temp = extrusion_temp

        self.emmisivity = emissivity

        self._adhesion_model = adhesion_model

    @property
    def name(self) -> str:
        """Name."""
        return self._long_name

    @name.setter
    def name(self, value) -> None:
        if not isinstance(value, str):
            raise TypeError('Name must be a string.')

        self._long_name = value

    @property
    def density(self) -> float:
        """Density [kg/cu.m]"""
        return self._density

    @density.setter
    def density(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Density must be a numeric type, not a {type(value)}')
        elif float(value) <= 0.0:
            raise ValueError('Density must be greater than zero.')

        self._density = float(value)

    @property
    def thermal_conductivity(self) -> float:
        """Thermal conducitivity [W/m-K]."""
        return self._thermal_conductivity

    @thermal_conductivity.setter
    def thermal_conductivity(self, value) -> None:
        if not isinstance(value, float):
            raise TypeError(f'Thermal conducitivity must be a numeric type, not a {type(value)}')
        elif float(value) <= 0:
            raise ValueError('Thermal conducitivity must be greater than zero.')

        self._thermal_conductivity = float(value)

    @property
    def specific_heat_capacity(self) -> float:
        """Specific heat capacity []."""
        return self._specific_heat_capacity

    @specific_heat_capacity.setter
    def specific_heat_capacity(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Specific heat capacity must be type float, not a {type(value)}')
        elif float(value) <= 0.0:
            raise ValueError('Specific heat capacity must be greater than zero.')

        self._specific_heat_capacity = float(value)

    @property
    def glass_transition(self) -> float:
        """Glass transition [degC]."""
        return self._glass_transition

    @glass_transition.setter
    def glass_transition(self, value) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Glass transition temperature must be a numeric type, not a {type(value)}')
        elif value + CELSIUS_TO_KELVIN_OFFSET < 0:
            raise ValueError('Glass transition temperature must be greater than absolute zero.')

        self._glass_transition = float(value)

    @property
    def melt_transition(self) -> float | None:
        """Melt transition [degC]."""
        return self._melt_transition

    @melt_transition.setter
    def melt_transition(self, value) -> None:
        if value is None:
            self._melt_transition = None
            return

        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Melt transition temperature must be a numeric type, not a {type(value)}')
        elif value <= -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Melt transition temperature must be greater than absolute zero.')

        self._melt_transition = float(value)

    @property
    def extrusion_temp(self) -> float:
        """Extruder temperature setpoint [degC]."""
        return self._extrusion_temp

    @extrusion_temp.setter
    def extrusion_temp(self, value: float) -> None:
        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Extruder temperature setpoint must be a numeric type, not a {type(value)}')
        elif value <= -CELSIUS_TO_KELVIN_OFFSET:
            raise ValueError('Extruder temperature setpoint must be greater than absolute zero.')

        self._extrusion_temp = float(value)

    @property
    def emissivity(self) -> float:
        """Emissivity [0.0-1.0]."""
        return self._emmisivity

    @emissivity.setter
    def emissivity(self, value: float) -> None:
        if value is None:
            self._emmisivity = None
            return

        if not isinstance(value, NUMERICAL_TYPES):
            raise TypeError(f'Emissivity must be a numeric type, not a {type(value)}')
        elif not (0.0 <= value <= 1.0):
            raise ValueError('Emissivity must be between 0 and 1.')

        self._emmisivity = float(value)

    @property
    def volumetric_heat_capacity(self) -> float:
        """Volumetric heat capacity [J/cu.m-K]."""
        return self._density * self._specific_heat_capacity  # [J/m3-K]

    @property
    def thermal_effusivity(self) -> float:
        """Thermal effusivity []. """
        return math.sqrt(self._thermal_conductivity * self.volumetric_heat_capacity)

    @property
    def thermal_diffusivity(self) -> float:
        """Thermal diffusivity [sq.m/s]."""
        return self._thermal_conductivity / self.volumetric_heat_capacity  # [m2/s]

    @property
    def adhesion_model(self) -> (WilliamLandelFerryModel | ArrheniusModel | None):
        """Adhesion model instance."""
        return self._adhesion_model

    @adhesion_model.setter
    def adhesion_model(self, adhesion_model: WilliamLandelFerryModel | ArrheniusModel | None):
        if adhesion_model is None:
            adhesion_model = None
            return

        if not isinstance(adhesion_model, ADHESION_MODELS):
            raise type(f'Adhesion model must be an adhesion model instance, not {type(adhesion_model)}')

        self._adhesion_model = adhesion_model

    @property
    def has_adhesion_model(self) -> bool:
        return self._adhesion_model is not None

    def __str__(self) -> str:
        return self._long_name


def calculate_healing(material: Material, time_arr: np.ndarray, temp_arr: np.ndarray) -> float:
    """Calculate healing ratio using Yang-Pitchumani non-isothermal healing model. 

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

    relaxation_times = material.adhesion_model.calc_relaxation_time(temp_arr)

    z = np.divide(1.0, relaxation_times)

    z[temp_arr < material.glass_transition] = 0

    res = abs(np.sum(np.trapz(y=z, x=time_arr)))**0.25

    if res > 1:
        #lg.info('Estimated relative weld strength adhesion ratio exceeds 1 ( R = %s).', res)
        pass

    floor_res = min(res, 1)

    return floor_res


# https://thermtest.com/thermal-resources/materials-database

ABS = Material(
    name='ABS',
    density=1040,
    thermal_conductivity=0.209,
    specific_heat_capacity=1506,
    glass_transition=105,
    melt_transition=None,
    extrusion_temp=260,
    emissivity=None,
)

QSR = Material(
    name='QSR',
    density=1180,
    thermal_conductivity=0.1,
    specific_heat_capacity=1950,
    glass_transition=165,
    melt_transition=None,
    extrusion_temp=295,
    emissivity=None,
    adhesion_model=WilliamLandelFerryModel(
        c_1=5.78,
        c_2=182,
        t_ref=200,
        omega_r=4.608,
    ),
)

F375M = Material(
    name='F375M Sinterable',
    density=6174.2,
    thermal_conductivity=10.6,
    specific_heat_capacity=942.8,
    glass_transition=-30.0,
    melt_transition=170,
    extrusion_temp=235,
    emissivity=None,
    adhesion_model=ArrheniusModel(
        e_a=48656,
        t_ref=200.0,
        omega_r=76.923,
    ),
)

if __name__ == '__main__':
    print(f'\nModule {Path(__file__).name} is not intended to be run as standalone module.')
