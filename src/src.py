import logging as lg
import math
import matplotlib.pyplot as plt
import numpy as np

import hickson

from material import Material, QSR, F375M

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8

FLOAT64 = np.float64

lg.basicConfig(level=lg.INFO)


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
                C1, C2 = hickson.jump_match_coeff(K1, D1, K2, D1, h0, h1, h2, CONTACT_TRANSFER_COEFF)
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

    coeff[0:NODES_PER_LAYER - 1] = CONVECTION_COEFF / (M1.volumetric_heat_capacity * DELTA_Z)

    coeff[NODES_PER_LAYER:-1] = CONVECTION_COEFF / (M2.volumetric_heat_capacity * DELTA_Z)

    return coeff


def solve_system(M1: Material, M2: Material, temp: np.ndarray) -> np.ndarray:
    """_summary_

    Args:
        M1 (Material): _description_
        M2 (Material): _description_
        temp (np.ndarray): _description_

    Returns:
        np.ndarray: _description_
    """

    # Calculate the time step; see eq. 4 in Basgul1 et al.
    time_step = 0.9 * (0.5 * (DELTA_Z**2) / max(M1.thermal_diffusivity, M2.thermal_diffusivity))

    lg.info(f'time step = {time_step}')

    # Calculate the total number of time steps.
    time_steps = int(math.ceil(TIME / time_step))

    lg.info(f'time step cnt = {time_steps}')

    lg.info(f'simulation time = {float(time_steps) * time_step}')

    # Make the conduction coefficient matrix.
    coeff_k = mk_conduction_matrix(M1, M2)

    # Make the convection coefficient matrix.
    coeff_h = mk_convection_matrix(M1, M2)

    res = np.zeros(shape=[NODE_CNT, time_steps], dtype=FLOAT64)

    # Store the initial temperature profile.
    T = temp

    # Store the initial temperature state in the result array.
    res[:, 0] = T

    for i in range(1, time_steps):

        T = res[:, i - 1]

        res[:, i] = T + time_step * (np.dot(coeff_k, T.T) - coeff_h * (T - T_AMB))

    return res


print(' ')
print(' ')

T_AMB = 23

T_HOT = QSR.extrusion_temp

LAYER_THICKNESS = inch_to_millimeter(0.007)

LAYER_CNT = 10

NODES_PER_LAYER = 400

lg.info(f'layer thickness = {LAYER_THICKNESS}\nlayer count = {LAYER_CNT}\nnodes per layer = {NODES_PER_LAYER}')

NODE_CNT = LAYER_CNT * NODES_PER_LAYER

DELTA_Z = LAYER_THICKNESS / NODES_PER_LAYER

lg.info(f'node spacing = {DELTA_Z}')

CONVECTION_COEFF = 50

CONTACT_TRANSFER_COEFF = 1e10

TIME = 50.000

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

# Apply the hot temperature.
T[0:NODES_PER_LAYER] = T_HOT

res = solve_system(QSR, QSR, T)

print()

# plt.plot(res[NODES_PER_LAYER - 1 - 4, :])
# plt.plot(res[NODES_PER_LAYER - 1 - 3, :])
# plt.plot(res[NODES_PER_LAYER - 1 - 2, :])
# plt.plot(res[NODES_PER_LAYER - 1 - 1, :])
# plt.plot(res[NODES_PER_LAYER - 1 - 0, :])
# plt.plot(res[NODES_PER_LAYER - 1 + 1, :])
# plt.plot(res[NODES_PER_LAYER - 1 + 2, :])
# plt.plot(res[NODES_PER_LAYER - 1 + 3, :])
# plt.plot(res[NODES_PER_LAYER - 1 + 4, :])
# plt.plot(res[NODES_PER_LAYER - 1 + 5, :])

interface_temp = hickson.calc_interface_temp(QSR.thermal_effusivity, res[NODES_PER_LAYER - 1, :],
                                             QSR.thermal_effusivity, res[NODES_PER_LAYER, :])

plt.plot(interface_temp)

# plt.plot([QSR.extrusion_temp, QSR.extrusion_temp])
plt.show()
