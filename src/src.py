import json
import math

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8


class Material():

    def __init__(self, name, density, thermal_conductivity, specific_heat_capacity, glass_transition, melt_transition,
                 extrusion_temp):
        """_summary_

        Args:
            name (_type_): _description_
            density (_type_): _description_
            thermal_conductivity (_type_): _description_
            specific_heat_capacity (_type_): _description_
            glass_transition (_type_): _description_
            melt_transition (_type_): _description_
            extrusion_temp (_type_): _description_
        """
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


def calculate_interface_temperature(mat_1: Material, mat_2: Material, temp_1: float, temp_2: float) -> float:

    num = mat_1.thermal_effusivity * temp_1 + mat_2.thermal_effusivity * temp_2

    den = mat_1.thermal_effusivity + mat_2.thermal_effusivity

    return num / den


def inch_to_millimeter(arg: float, n: int = 1) -> float:
    """Convert an argument in inches to millimeters. This function accepts an 
    optional argument `n`, which specifies the power. For example, a value of 
    `n = 1` converts the argument from inches to millimeters. A value of `n = 3` 
    converts the argument from cubic inches to cubic millimeters. This function
    can also be used to do the reverse conversion, as, for example, `n = -1` 
    converts the argument from millimeteres to inches. """

    if isinstance(arg) != int:
        raise TypeError('arg must be type int.')

    if arg == 0:
        raise ValueError('n = 0 will have no effect.')

    res = arg * (25.4**n)

    return res


def import_materials(fpath: str) -> dict[Material]:

    return []

layer_count = 10

nodes_per_layer = 10

node_spacing = layer_count / nodes_per_layer

time_spacing = 0

mat = Material('tst', 1000, 0.200, 2000, 180, None, 300)
