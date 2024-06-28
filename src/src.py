import json
import math
import matplotlib.pyplot as plt
import numpy as np

from types import NoneType

import hickson

from material import Material, QSR, F375M

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8

FLOAT64 = np.float64


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

    if LAYER_CNT < 2:
        raise ValueError('At least 2 layers are needed.')

    if NODES_PER_LAYER < 3:
        raise ValueError('Node count must be greater than or equal to 3.')

    # Node-to-node distance.
    h0 = DELTA_Z

    # Node-to-interface distance.
    h1 = h2 = h0 / 2

    # Store params for convenience.
    D1 = M1.thermal_diffusivity
    K1 = M1.thermal_conductivity
    D2 = M2.thermal_diffusivity
    K2 = M2.thermal_conductivity

    # Initialize the conduction coefficient matrix.
    coeff = np.zeros(shape=[NODE_CNT, NODE_CNT], dtype=FLOAT64)

    # Initialize a convenience array for second order finite difference approximation matrix.
    ARR = np.array([1, -2, 1], dtype=FLOAT64)

    for i in range(0, NODE_CNT):
        match i:
            case 0:
                coeff[i, i:i + 3] = (D1 / h0**2) * ARR  # fwd

            case i if i > 0 and i < NODES_PER_LAYER - 1:
                coeff[i, i - 1:i + 2] = (D1 / h0**2) * ARR

            case i if i == NODES_PER_LAYER - 1:
                C1, C2 = hickson.jump_match_coeff(K1, D1, K2, D1, h0, h1, h2, 1e10)
                coeff[i, i - 1:i + 3] = C1
                coeff[i + 1, i - 1:i + 3] = C2

            case i if i == NODES_PER_LAYER:
                pass  # This nodes coefficients were already set in i == NODES_PER_LAYER - 1.

            case i if i > NODES_PER_LAYER - 2 and i < NODE_CNT - 1:
                coeff[i, i - 1:i + 2] = (D2 / h0**2) * ARR

            case i if i == NODE_CNT - 1:
                pass
                #coeff[i, i - 2:i + 1] = (D2 / h0**2) * ARR

            case _:
                raise ValueError('Illegal index.')

    return coeff


def mk_convection_matrix(M1: Material, M2: Material) -> np.array:
    """This function generates a Nx1 matrix of convection coefficients. 

    Args:
        M1 (Material): _description_
        M2 (Material): _description_

    Returns:
        np.array: Convection FD coefficients. 
    """

    coeff = np.zeros(shape=NODE_CNT, dtype=FLOAT64)

    coeff[0:NODES_PER_LAYER - 1] = CONVECTION_COEFF / (M1.volumetric_heat_capacity * LAYER_THICKNESS)

    coeff[NODES_PER_LAYER:-1] = CONVECTION_COEFF / (M2.volumetric_heat_capacity * LAYER_THICKNESS)


T_AMB = 23

T_HOT = QSR.extrusion_temp

LAYER_THICKNESS = 0.001

LAYER_CNT = 10

NODES_PER_LAYER = 5

NODE_CNT = LAYER_CNT * NODES_PER_LAYER

print(NODE_CNT)

DELTA_Z = LAYER_THICKNESS / NODES_PER_LAYER

print(DELTA_Z)

CONVECTION_COEFF = 10

dt = 0.001

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

# Apply the hot temperature.
T[0:NODES_PER_LAYER] = T_HOT

res = np.zeros(shape=[NODE_CNT, int(5e3)], dtype=FLOAT64)

res[:, 0] = T

np.set_printoptions(precision=6, linewidth=400)

c_k = mk_conduction_matrix(F375M, QSR)

c_c = mk_convection_matrix(F375M, QSR)

for i in range(1, int(5e3)):

    T = res[:, i - 1]

    res[:, i] = T + (dt * np.dot(c_k, T.T))
    
    

plt.plot(res[0, :])
plt.plot(res[1, :])
plt.plot(res[2, :])
plt.plot(res[3, :])
plt.plot(res[4, :])
plt.plot(res[5, :])
plt.plot(res[6, :])
plt.plot(res[7, :])
plt.plot(res[8, :])
plt.plot(res[9, :])
plt.plot([QSR.extrusion_temp, QSR.extrusion_temp])
plt.show()
