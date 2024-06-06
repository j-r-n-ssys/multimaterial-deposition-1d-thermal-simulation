import math

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8


class Material():

    def __init__(self, name, density, thermal_conductivity, specific_heat_capacity, glass_transition, melt_transition,
                 extrusion_temp):
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


mat = Material('tst', 1000, 0.200, 2000, 180, None, 300)
