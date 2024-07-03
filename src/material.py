"""Module to contain the common material class."""

import math

from types import NoneType

class Material():
    """This object represents a material."""

    def __init__(self, name: str, density: float, thermal_conductivity: float, specific_heat_capacity: float,
                 glass_transition: float, melt_transition: float, extrusion_temp: float, emissivity: float):
        """Init."""

        if not isinstance(name, str):
            raise TypeError('Name must be a string.')

        self.long_name = str(name)

        if not isinstance(density, (int, float)):
            raise TypeError('Density must be type float.')

        if density < 0.0:
            raise ValueError('Density must be greater than zero.')

        self.density = float(density)  # [kg/m3]

        if not isinstance(thermal_conductivity, (int, float)):
            raise TypeError('Thermal conductivity must be type float.')

        if thermal_conductivity < 0.0:
            raise ValueError('Thermal conducitivty must be greater than zero.')

        self.thermal_conductivity = float(thermal_conductivity)  # [W/m-K]

        if not isinstance(specific_heat_capacity, (int, float)):
            raise TypeError('Specific heat capacity must be type float.')

        if specific_heat_capacity < 0.0:
            raise ValueError('Specific heater capacity must be greater than zero.')

        self.specific_heat_capacity = float(specific_heat_capacity)

        if not isinstance(glass_transition, (int, float)):
            raise TypeError('Glass transition temperature must be type float.')

        self.glass_transition = float(glass_transition)  # [K]

        if not isinstance(melt_transition, (NoneType, int, float)):
            raise TypeError('Melt temperature must be type float.')

        self.melt_transition = float(melt_transition) if not isinstance(melt_transition, NoneType) else None

        if not isinstance(extrusion_temp, (int, float)):
            raise TypeError('Extrusion temperature must be type float.')

        self.extrusion_temp = float(extrusion_temp)  # [K]

        if not isinstance(emissivity, (NoneType, int, float)):
            raise TypeError('Emissivity must be type float.')

        if not isinstance(emissivity, NoneType) and (emissivity < 0.0 or emissivity > 1.0):
            raise ValueError('Emissivity must be between 0 and 1.')

        self.emmisivity = float(emissivity) if not isinstance(emissivity, NoneType) else None  # [0.0-1.0]

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
    
# https://thermtest.com/thermal-resources/materials-database
    
ABS = Material('ABS', 1040, 0.209, 1506, 105, None, 260, None)

QSR = Material('QSR', 1180, 10.6, 943, 165, None, 295, None)

F375M = Material('F375M Sinterable', 6212.8, 10.6, 943, 30, 170, 235, None)
