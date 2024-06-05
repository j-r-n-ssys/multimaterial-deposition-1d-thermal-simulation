import math

class Material():

    def __init__(self, name, density, thermal_conductivity, specific_heat_capacity, glass_transition, melt_transition, extrusion_temp):
        self.name: str = name
        self.density: float = density
        self.thermal_conductivity: float = thermal_conductivity
        self.specific_heat_capacity:float = specific_heat_capacity
        self.volumetric_heat_capacity:float = self.density * self.specific_heat_capacity
        self.thermal_effusivity: float = math.sqrt(self.thermal_conductivity * self.volumetric_heat_capacity)
        self.glass_transition: float = glass_transition
        self.melt_transition: float = melt_transition
        self.extrusion_temp: float = extrusion_temp
        self.thermal_diffusivity: float = self.thermal_conductivity / self.volumetric_heat_capacity
        


# notes
# density
# heat capacity
# specific heat capacity
# melt temperature
# thermal conducitivity
# thermal effusicity
