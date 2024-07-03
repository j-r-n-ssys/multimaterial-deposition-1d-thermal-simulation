"""_missing_docstring_"""

import logging as lg
import math
import matplotlib.pyplot as plt
import numpy as np

import hickson

from material import Material, QSR, F375M

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
                
                print(d_1 / h0**2)

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
                pass
                #coeff[i, i - 2:i + 1] = (D2 / h0**2) * ARR

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


def solve_system(m_1: Material, m_2: Material, temp: np.ndarray, time_step:float = None) -> tuple[np.ndarray, np.ndarray]:  #pylint: disable=redefined-outer-name, disable=line-too-long
    """_summary_

    Args:
        M1 (Material): _description_
        M2 (Material): _description_
        temp (np.ndarray): _description_

    Returns:
        np.ndarray: _description_
    """
    
    if time_step is None:
        time_step = (0.5 * (DELTA_Z**2) / max(m_1.thermal_diffusivity, m_2.thermal_diffusivity)) # Calculate the time step; see eq. 4 in Basgul1 et al.
    else:
        
        
        
        if max(m_1.thermal_diffusivity, m_2.thermal_diffusivity) * time_step / DELTA_Z**2 > 0.5:
            lg.warning('Time step is too large to guarentee simulation stability.')
            print(time_step * max(m_1.thermal_diffusivity, m_2.thermal_diffusivity) / DELTA_Z**2)

    lg.info(pad_kv_pair_str('time step', time_step))

    # Calculate the total number of time steps.
    time_steps = int(math.ceil(TIME_F / time_step))

    lg.info(pad_kv_pair_str('time step cnt', time_steps))

    lg.info(pad_kv_pair_str('simulation time', float(time_steps) * time_step))

    # Make the conduction coefficient matrix.
    coeff_k = mk_conduction_matrix(m_1, m_2)

    # Make the convection coefficient matrix.
    coeff_h = mk_convection_matrix(m_1, m_2)

    temp_arr = np.zeros(shape=[NODE_CNT, time_steps], dtype=FLOAT64)  #pylint disable=redefined_outer_scope

    # Store the initial temperature state in the result array.
    temp_arr[:, 0] = temp

    for i in range(1, time_steps):

        temp = temp_arr[:, i - 1]

        temp_arr[:, i] = temp + time_step * (np.dot(coeff_k, temp.T) - coeff_h * (temp - T_AMB))

    time_arr = np.linspace(TIME_0, TIME_F, time_steps)

    return time_arr, temp_arr


def prep_next_temp_profile(curr_profile: np.ndarray) -> np.ndarray:
    """_summary_

    Args:
        curr_profile (np.ndarray): _description_

    Returns:
        np.ndarray: _description_
    """

    curr_profile[NODES_PER_LAYER:-1] = curr_profile[0:-NODES_PER_LAYER - 1]



print(' ')
print(' ')

T_AMB = 115

T_HOT = QSR.extrusion_temp

LAYER_THICKNESS = 1e-3 * inch_to_millimeter(0.010)

LAYER_CNT = 10

NODES_PER_LAYER = 10

pad_kv_pair_str('layer thickness', LAYER_THICKNESS)

lg.info(pad_kv_pair_str('layer thickness', LAYER_THICKNESS))
lg.info(pad_kv_pair_str('layer count', LAYER_CNT))
lg.info(pad_kv_pair_str('nodes per layer', NODES_PER_LAYER))

NODE_CNT = LAYER_CNT * NODES_PER_LAYER

DELTA_Z = LAYER_THICKNESS / NODES_PER_LAYER

lg.info(pad_kv_pair_str('node spacing', DELTA_Z))

CONVECTION_COEFF = 200

CONTACT_TRANSFER_COEFF = 1e10

TIME_0 = 0.0

TIME_F = 10.000

T = np.array([T_AMB] * NODE_CNT, dtype=FLOAT64)

# Apply the hot temperature.
T[0:NODES_PER_LAYER] = T_HOT

m_1 = QSR

m_2 = QSR

t, res = solve_system(m_1, m_2, T)

print()

interface_temp = hickson.calc_interface_temp(m_1.thermal_effusivity, res[NODES_PER_LAYER - 1, :],
                                             m_2.thermal_effusivity, res[NODES_PER_LAYER, :])

prep_next_temp_profile(res[:, -1])

plt.semilogx(t, res[0,:])
plt.semilogx(t, interface_temp)
plt.semilogx(t, np.array([m_1.glass_transition] * len(t), dtype=FLOAT64))
plt.semilogx(t, np.array([T_AMB] * len(t), dtype=FLOAT64))
# plt.plot([QSR.extrusion_temp, QSR.extrusion_temp])

plt.show()

print(QSR.thermal_diffusivity)

print(F375M.thermal_diffusivity)