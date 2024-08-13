"""_missing_docstring_"""

import logging as lg
import math
import matplotlib.pyplot as plt
import numpy as np

import hickson

from material import Material, calculate_healing, QSR, F375M

from units import inch_to_millimeter
from util import pad_kv_pair_str

STEFAN_BOLTZMANN_CONSTANT = 5.670374419 * 1e-8

FLOAT64 = np.float64

lg.basicConfig(level=lg.INFO)


def mk_conduction_matrix(m_1: Material, m_2: Material) -> np.ndarray:  #pylint: disable=redefined-outer-name
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
    d_1 = m_1.thermal_diffusivity
    k_1 = m_1.thermal_conductivity
    d_2 = m_2.thermal_diffusivity
    k_2 = m_2.thermal_conductivity

    # Initialize the conduction coefficient matrix.
    coeff = np.zeros(shape=[NODE_CNT, NODE_CNT], dtype=FLOAT64)

    # Initialize a convenience array for second order finite difference approximation matrix.
    fd_arr = np.array([1, -2, 1], dtype=FLOAT64)

    for i in range(0, NODE_CNT):
        match i:
            case 0:
                coeff[i, i:i + 3] = (d_1 / h0**2) * fd_arr  # fwd

            case i if i > 0 and i < NODES_PER_LAYER - 1:
                coeff[i, i - 1:i + 2] = (d_1 / h0**2) * fd_arr

            case i if i == NODES_PER_LAYER - 1:
                c_1, c_2 = hickson.jump_match_coeff(k_1, d_1, k_2, d_1, h0, h1, h2, CONTACT_TRANSFER_COEFF)
                coeff[i, i - 1:i + 3] = c_1
                coeff[i + 1, i - 1:i + 3] = c_2

            case i if i == NODES_PER_LAYER:
                pass  # This nodes coefficients were already set in i == NODES_PER_LAYER - 1.

            case i if i > NODES_PER_LAYER - 2 and i < NODE_CNT - 1:
                coeff[i, i - 1:i + 2] = (d_2 / h0**2) * fd_arr

            case i if i == NODE_CNT - 1:
                pass  # coeff[i, i - 2:i + 1] = (d_2 / h0**2) * fd_arr

            case _:
                raise ValueError('Illegal index.')

    return coeff


def mk_convection_matrix(m_1: Material, m_2: Material) -> np.array:  #pylint: disable=redefined-outer-name
    """This function generates a Nx1 matrix of convection coefficients. 

    Args:
        M1 (Material): _description_
        M2 (Material): _description_

    Returns:
        np.array: Convection FD coefficients. 
    """

    cv_1 = m_1.volumetric_heat_capacity

    cv_2 = m_2.volumetric_heat_capacity

    coeff = np.zeros(shape=NODE_CNT, dtype=FLOAT64)

    coeff[0:NODES_PER_LAYER - 1] = CONVECTION_COEFF / (cv_1 * DELTA_Z)

    coeff[NODES_PER_LAYER:-1] = CONVECTION_COEFF / (cv_2 * DELTA_Z)

    return coeff


def solve_system(
    m_1: Material,
    m_2: Material,
    temp: np.ndarray,
    time_step: float | None = None
) -> tuple[np.ndarray, np.ndarray]:  #pylint: disable=redefined-outer-name, disable=line-too-long
    """Solve the 1D heat transfer problem.

    Args:
        M1 (Material): Material 1 instance.
        M2 (Material): Material 2 instance.
        temp (np.ndarray): Temperature array [degC].
        time_step (float | none): Simulation time step [s]. If None, a timestep will be calculated using a simulation
        stability criteria. Defauls to None. 

    Returns:
        tuple[np.ndarray, np.ndarray]: Tuple of time array, time-temperature array. 
    """

    max_thermal_diffusivity = max(m_1.thermal_diffusivity, m_2.thermal_diffusivity)

    if time_step is None:
        # Calculate the time step; see eq. 4 in Basgul1 et al.
        time_step = 0.5 * (0.5 * DELTA_Z**2 / max_thermal_diffusivity)
    else:
        if max_thermal_diffusivity * time_step / DELTA_Z**2 > 0.5:
            lg.warning('Time step is too large to guarantee simulation stability: %s > %s', time_step,
                       0.5 * DELTA_Z**2 / max_thermal_diffusivity)

    lg.info(pad_kv_pair_str('time step', time_step))

    # Calculate the total number of time steps.
    time_steps = int(math.ceil(TIME_F / time_step))

    lg.info(pad_kv_pair_str('time step cnt', time_steps))

    lg.info(pad_kv_pair_str('simulation time', float(time_steps) * time_step))

    # Make the n x n conduction coefficients matrix.
    coeff_k = mk_conduction_matrix(m_1, m_2)

    # Make the n x 1 convection coefficients matrix.
    coeff_h = mk_convection_matrix(m_1, m_2)

    temp_arr = np.zeros(shape=[NODE_CNT, time_steps], dtype=FLOAT64)  #pylint disable=redefined_outer_scope

    # Store the initial temperature state in the result array.
    temp_arr[:, 0] = temp

    for i in range(1, time_steps):

        temp = temp_arr[:, i - 1]

        temp_arr[:, i] = temp + time_step * (np.dot(coeff_k, temp.T) - coeff_h * (temp - T_AMB))

    time_arr = np.linspace(0, TIME_F, time_steps)

    return time_arr, temp_arr


def prep_next_temp_profile(curr_profile: np.ndarray, next_maerial: Material) -> np.ndarray:
    """Calculate the next temperature profile. 

    Args:
        curr_profile (np.ndarray): Current temperature profile. 

    Returns:
        np.ndarray: _description_
    """

    curr_profile[NODES_PER_LAYER:-1] = curr_profile[0:-NODES_PER_LAYER - 1]

    curr_profile[0:NODES_PER_LAYER] = next_maerial.extrusion_temp

    return curr_profile


def calc_interface_temperature(
        m_1: Material,
        m_2: Material,
        t_arr: np.ndarray,
        algorithm: str = 'effusivity') -> np.ndarray:  #pylint: disable:line-too-long, unused-argument
    """Calculate the interface temperature.

    Args:
        m_1 (Material): _description_
        m_2 (Material): _description_
        t_arr (np.ndarray): _description_
        algorithm (str, optional): Interface calculation algorithm: 'average' or 'effusivity'. Defaults to 'effusivity'.

    Returns:
        np.ndarray: Interface temperature array.
    """

    t_1 = t_arr[NODES_PER_LAYER - 1, :]

    t_2 = t_arr[NODES_PER_LAYER + 0, :]

    match algorithm:
        case 'average':
            t_int = (t_1 + t_2) / 2
        case 'effusivity':
            e_1 = m_1.thermal_effusivity
            e_2 = m_2.thermal_effusivity
            t_int = hickson.calc_interface_temp(e_1, t_1, e_2, t_2)
        case _:
            raise ValueError('Invalid interface temperature calculcation algorithm.')

    return t_int


print(' ')
print(' ')

LAYER_THICKNESS = 1e-3 * inch_to_millimeter(0.010)

LAYER_WIDTH = 1e-3 * inch_to_millimeter(0.030)

LAYER_CNT = 10

NODES_PER_LAYER = 20

pad_kv_pair_str('layer thickness', LAYER_THICKNESS)

lg.info(pad_kv_pair_str('layer thickness', LAYER_THICKNESS))
lg.info(pad_kv_pair_str('layer count', LAYER_CNT))
lg.info(pad_kv_pair_str('nodes per layer', NODES_PER_LAYER))

NODE_CNT = LAYER_CNT * NODES_PER_LAYER

DELTA_Z = LAYER_THICKNESS / NODES_PER_LAYER

lg.info(pad_kv_pair_str('node spacing', DELTA_Z))

CONVECTION_COEFF = 50

CONTACT_TRANSFER_COEFF = 1e10

TIME_F = 10.000

T_AMB = 115

m_top = QSR

m_bot = QSR

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

T = prep_next_temp_profile(T, m_top)

t, res = solve_system(m_top, m_bot, T)

interface_temp = calc_interface_temperature(m_top, m_bot, res, algorithm='average')

lg.info('weld strength = %s', calculate_healing(m_top, t, interface_temp))

plt.semilogx(t, interface_temp, label='S-S')

m_top = F375M

m_bot = F375M

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

T[0:NODES_PER_LAYER] = m_top.extrusion_temp

t, res = solve_system(m_top, m_bot, T)

interface_temp = calc_interface_temperature(m_top, m_bot, res, algorithm='average')

lg.info('weld strength = %s', calculate_healing(m_top, t, interface_temp))

plt.semilogx(t, interface_temp, label='M-M')

m_top = F375M

m_bot = QSR

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

T[0:NODES_PER_LAYER] = m_top.extrusion_temp

t, res = solve_system(m_top, m_bot, T, time_step=2.5e-6)

interface_temp = calc_interface_temperature(m_top, m_bot, res, algorithm='average')

lg.info('weld strength = %s', calculate_healing(m_top, t, interface_temp))

plt.semilogx(t, interface_temp, label='M-S')

m_top = QSR

m_bot = F375M

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

T[0:NODES_PER_LAYER] = m_top.extrusion_temp

t, res = solve_system(m_top, m_bot, T)

interface_temp = calc_interface_temperature(m_top, m_bot, res, algorithm='average')

lg.info('weld strength = %s', calculate_healing(m_top, t, interface_temp))

plt.semilogx(t, interface_temp, label='S-M')

plt.legend()

plt.xlabel('Time since deposition [s]')

plt.ylabel('Interface temperature [degC]')

plt.xlim([1e-3, 1e1])

plt.ylim([0, 300])

plt.show()
