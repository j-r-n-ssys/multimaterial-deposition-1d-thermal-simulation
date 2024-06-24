import json
import math
import numpy as np

from types import NoneType

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8

FLOAT64 = np.float64


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

    if not isinstance(T1, (int, float)) or not isinstance(T2, (int, float)):
        raise TypeError('T1 and T2 must be type float.')

    if not isinstance(M1, Material) or not isinstance(M2, Material):
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

    if not isinstance(f, (float, int)):
        raise TypeError('Argument f must be of type int or float).')

    if not isinstance(n, int):
        raise TypeError('Argument n must be type int.')

    if f == 0:
        raise ValueError('n = 0 has no effect.')

    res = float(f) * (25.4**n)

    return res


def mk_conduction_matrix(M1: Material, M2: Material) -> np.ndarray:
    """---IN DEVELOPMENT---"""

    def mk_interface_kernal(D1, K1, D2, K2, dZ2):

        Yp2 = (2 * K1 + 3 * K2) * D1 - D2 * K1

        Yp1 = 4 * D2 * K1 - 2 * (K1 + 3 * K2) * D1

        Ym1 = 4 * D1 * K2 - 2 * (3 * K1 + K2) * D2

        Ym2 = (3 * K1 + 2 * K2) * D2 - D1 * K2

        den = (6 * (K1 + K2) * dZ2)

        fda = np.array([Yp2, Yp1, 0, Ym1, Ym2], dtype=FLOAT64) / den  # finite difference approximation

        return fda

    dZ = DELTA_Z

    dZ2 = dZ**2

    D1 = M1.thermal_diffusivity

    K1 = M1.thermal_conductivity

    D2 = M2.thermal_diffusivity

    K2 = M2.thermal_conductivity

    Z = np.zeros(shape=[node_cnt, node_cnt], dtype=FLOAT64)

    ARR = np.array([1, -2, 1], dtype=FLOAT64)

    for i in range(0, node_cnt):
        match i:

            case 0:
                Z[i, i:i + 3] = (D1 / dZ2) * ARR  # fwd

            case i if (i > 0) & (i < NODES_PER_LAYER_CNT - 1):
                Z[i, i - 1:i + 2] = (D1 / dZ2) * ARR

            case i if i == NODES_PER_LAYER_CNT - 1:
                Z[i, i - 2:i - 2 + 5] = mk_interface_kernal(D1, K1, D2, K2, dZ2)

            case i if (i > NODES_PER_LAYER_CNT - 1) & (i < node_cnt - 1):
                Z[i, i - 1:i + 2] = (D2 / dZ2) * ARR
            case i if i == node_cnt - 1:

                Z[i, i - 2:i + 1] = (D2 / dZ2) * ARR  # bkwd

            case _:
                raise ValueError('Illegal index reached.')
        print(i)

    print(Z)


T_AMB = 200

T_HOT = 400

LAYER_CNT = 2

NODES_PER_LAYER_CNT = 5

node_cnt = LAYER_CNT * (NODES_PER_LAYER_CNT - 1) + 1

DELTA_Z = 0.001

M1 = Material('QSR', 1180, 0.1, 1950, 165, None, 295, None)

T = np.array([T_AMB] * node_cnt, dtype=FLOAT64)

res = np.zeros(shape=[node_cnt, int(1e6)], dtype=FLOAT64)

# Apply the hot temperature.
T[0:NODES_PER_LAYER_CNT - 1] = T_HOT

# Set the interface temperature
T[NODES_PER_LAYER_CNT - 1] = calculate_interface_temperature(M1, M1, T_HOT, T_AMB)

print(T)

mk_conduction_matrix(M1, M1)
