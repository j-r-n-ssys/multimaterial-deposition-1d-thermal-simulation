import json
import math
import numpy as np

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8


class Material():

    def __init__(self, name: str, density: float, thermal_conductivity: float, specific_heat_capacity: float,
                 glass_transition: float, melt_transition: float, extrusion_temp: float):
        self.name: str = name
        self.density: float = density  # [kg/m3]
        self.thermal_conductivity: float = thermal_conductivity  # [W/m-K]
        self.specific_heat_capacity: float = specific_heat_capacity
        self.glass_transition: float = glass_transition  # [K]
        self.melt_transition: (float | None) = melt_transition  # [K]
        self.extrusion_temp: float = extrusion_temp  # [K]
        self.emmisivity: float = 0  #[0.0-1.0]

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


class Later():

    def __init__(self, width: float, height: float):
        """_summary_

        Args:
            width (float): [m]
            height (float): [m]
        """
        self.width: float = width
        self.height: float = height

    @property
    def area(self) -> float:
        """Approximate layer area. """
        return self.width * self.height


def calculate_interface_temperature(M1: Material, M2: Material, T1: float, T2: float) -> float:
    """Calculate interface temperature between two semi-infinite bodies using 
    thermal effusivity. 

    Args:
        M1 (Material): Body 1 material parameter object. 
        M2 (Material): Body 2 material parameter object. 
        T1 (float): Body 1 temperature at interface [degC]. 
        T2 (float): Body 2 temperature at interface [degC]. 

    Returns:
        float: Interface temperature [degC].
    """

    if isinstance(T1) != float or isinstance(T2) != float:
        raise TypeError('T1 and T2 must be type float.')

    if isinstance(M1) != Material or isinstance(M2) != Material:
        raise TypeError('M1 and M2 must be type Material.')

    # Handle no temperature change case.
    if T1 == T2:
        return T1

    # Store M1s effusivity.
    e1 = M1.thermal_effusivity

    if e1 == 0:
        raise ValueError('M1 effusivity is zero.')

    # Store M2s effusivity.
    e2 = M2.thermal_effusivity

    if e2 == 0:
        raise ValueError('M2 effusivity is zero.')

    num = (e1 * T1) + (e2 * T2)

    den = e1 + e2

    Ti = num / den

    return Ti


def inch_to_millimeter(f: float, n: int = 1) -> float:
    """Convert an argument `f` in inches to millimeters. This function accepts
    an optional argument `n`, which specifies the power. For example, a value of 
    `n = 1` converts the argument from inches to millimeters. A value of `n = 3` 
    converts the argument from cubic inches to cubic millimeters. This function
    can also be used to do the reverse conversion, as, for example, `n = -1` 
    converts the argument from millimeteres to inches. """

    if isinstance(f) != float or isinstance(f) != int:
        raise TypeError('Argument f must be of type int or float).')

    if isinstance(n) != int:
        raise TypeError('Argument n must be type int.')

    if f == 0:
        raise ValueError('n = 0 has no effect.')

    res = float(f) * (25.4**n)

    return res


LAYER_CNT = 10

NODES_PER_LAYER_CNT = 10

node_cnt = LAYER_CNT * NODES_PER_LAYER_CNT

node_spacing = 0

time_spacing = 0

mat = Material('tst', 1000, 0.200, 2000, 180, None, 300)
