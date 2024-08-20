"""This module is a collection of finite difference coefficient generation 
functions derived by Hickson et al.

R. I. Hickson, S. I. Barry, G. N. Mercer, and H. S. Sidhu, “Finite 
difference schemes for multilayer diffusion,” Mathematical and computer 
modelling, vol. 54, no. 1-2, pp. 210-220, Jul. 2011, doi: 
https://doi.org/10.1016/j.mcm.2011.02.003."""

import numpy as np

from os.path import basename as get_module_fname
from pathlib import Path

FLOAT64 = np.float64

DEBUG = False

NUMERICAL_TYPES = (int, float)


def cond_match_coeff(k_1: float, d_1: float, k_2: float, d_2: float, h: float) -> np.ndarray:
    """Calculate the second order finite difference (FD) approximation 
    coefficients at the two-body interface using conductivity matching. A 
    fundamental assumption is that the interface is co-located with the middle
    node of this FD approximation. 

    Args:
        k_1 (float): Body 1 thermal conductivity [W/m-K].
        d_1 (float): Body 1 thermal diffusivity [sq.m/s].
        k_2 (float): Body 2 thermal conductivity [W/m-K]. 
        d_2 (float): Body 2 thermal diffusivity [sq.m/s].
        h (float): Finite difference node-to-node distance [m].

    Returns:
        np.ndarray: Finite difference (FD) approximation coefficients.
    """

    if not isinstance(k_1, NUMERICAL_TYPES) or not isinstance(k_1, NUMERICAL_TYPES):
        raise TypeError('Conducitivities K1 and K2 must be a numerical type.')
    elif k_1 <= 0 or k_2 <= 0:
        raise ValueError('Conducitivties K1 and K2 must be greater than zero.')

    if not isinstance(d_1, NUMERICAL_TYPES) or not isinstance(d_2, NUMERICAL_TYPES):
        raise TypeError('Diffusivities D1 and D2 must be a numerical type.')
    elif d_1 <= 0 or d_2 <= 0:
        raise ValueError('Diffusivities D1 and D2 must be greater than zero.')

    if not isinstance(h, NUMERICAL_TYPES):
        raise TypeError('Node-to-node distance h must be a numerical type.')
    elif h <= 0:
        raise ValueError('Node-to-node distance h must be greater than zero.')

    # Initialize the 5 x 1 coefficient matrix.
    coeff = np.zeros(shape=5, dtype=FLOAT64)

    # Calculate
    den = 6 * (k_1 + k_2) * h**2

    print(den)

    coeff[0] = ((2 * k_1 + 3 * k_2) * d_1 - d_2 * k_1) / den

    coeff[1] = (4 * d_2 * k_1 - 2 * (k_1 + 3 * k_2) * d_1) / den

    coeff[3] = (4 * d_1 * k_2 - 2 * (3 * k_1 + k_2) * d_2) / den

    coeff[4] = ((3 * k_1 + 2 * k_2) * d_2 - d_1 * k_2) / den

    if DEBUG:
        print(coeff)

    return coeff


def jump_match_coeff(k_1: float, d_1: float, k_2: float, d_2: float, h0: float, h1: float, h2: float,
                     H: float) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the second order finite difference (FD) approximation 
    coefficients on both sides of a two-body interface using jump matching. A 
    fundamental assumption is that the interface is not co-located with a node.

    Args:
        k_1 (float): Body 1 thermal conductivity [W/m-K]. 
        d_1 (float): Body 1 thermal diffusivity [sq.m/s].
        k_2 (float): Body 2 thermal conductivity [W/m-K]. 
        d_2 (float): Body 2 thermal diffusivity [sq.m/s]. 
        h0 (float): Node-to-node distance [m].
        h1 (float): Body 1 nearest-node-to-interface distance [m]. 
        h2 (float): Body 2 interface-to-nearest-node distance [m]
        H (float): Contact transfer coefficient [0, +inf]. Examples: H = 0: the interface is a perfect insulator, H -> +inf, the interface is in perfect contact. 

    Returns:
        tuple[np.ndarray, np.ndarray]: FD approximation coeffieints are node n-1 and node n+1
    """

    #
    #   h1
    #   D_i         K_i         K1
    #
    #   h2
    #   D_i+1       K_i+1       K2
    #
    #

    if DEBUG:
        print('')
        print('')

    if not isinstance(k_1, NUMERICAL_TYPES) or not isinstance(k_1, NUMERICAL_TYPES):
        raise TypeError('Conducitivities K1 and K2 must be a numerical type.')
    elif k_1 <= 0 or k_2 <= 0:
        raise ValueError('Conducitivties K1 and K2 must be greater than zero.')

    if not isinstance(d_1, NUMERICAL_TYPES) or not isinstance(d_2, NUMERICAL_TYPES):
        raise TypeError('Diffusivities D1 and D2 must be a numerical type.')
    elif d_1 <= 0 or d_2 <= 0:
        raise ValueError('Diffusivities D1 and D2 must be greater than zero.')

    if not isinstance(h0, NUMERICAL_TYPES):
        raise TypeError('Node-to-node distnce h0 must be a numerical type.')
    elif h0 <= 0:
        raise ValueError('Node-to-node distance h0 must be greater than zero.')

    if not isinstance(h1, NUMERICAL_TYPES) or not isinstance(h2, NUMERICAL_TYPES):
        raise TypeError('Node-to-interface distance h1 and h2 must be a numerical type.')
    elif h1 <= 0 or h2 <= 0:
        raise ValueError('Node-to-interface distances h1 and h2 must be greater than zero.')

    if (h1 + h2) != h0:
        raise ValueError('The sume of node-to-interface distances h1 and h2 must be equal node-to-node distance h0')

    if not isinstance(H, NUMERICAL_TYPES):
        raise TypeError('Contact transfer coefficient H must be a numerical type.')
    elif H < 0:
        raise TypeError('Contact transfer coefficient H must be greater than or equal to zero.')

    arr_m = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_2 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_2 = np.zeros(shape=4, dtype=FLOAT64)

    lam = h0 * (((h0 + 2 * h1) * (h0 + h2) * h2 * k_1 + (h0 + h1) * (h0 + 2 * h2) * h1 * k_2) * H + (h0 + 2 * h1) *
                (h0 + 2 * h2) * k_1 * k_2)

    #      0       -2        1
    #      1       -1        1
    #      2        1        2
    #      3        2        2

    arr_m_1[0] = (-1 * ((h0 + h2) * (h2 * H + k_2) + h2 * k_2) * h1**2 * k_1) / lam

    arr_m_1[1] = (((h0 + h2) * (h2 * H + k_2) + h2 * k_2) * (h0 + h1)**2 * k_1) / lam

    arr_m_1[2] = (h1 * (h0 + h1) * (h0 + h2)**2 * H * k_2) / lam

    arr_m_1[3] = (-1 * h1 * h2**2 * (h0 + h1) * H * k_2) / lam

    if DEBUG:
        print(f'arr_m_1: {arr_m_1}\n')

    den = (h0 * h1 * (h0 + h1))

    arr_m_2[0] = h1 / den

    arr_m_2[1] = -1 * (h0 + h1) / den

    if DEBUG:
        print(f'arr_m_2: {arr_m_2}\n')

    arr_m = d_1 * (h0 * arr_m_1 / den + arr_m_2)

    if DEBUG:
        print(f'arr_m:   {arr_m}\n')

    #      0       -2        1
    #      1       -1        1
    #      2        1        2
    #      3        2        2

    arr_p_1[0] = (-1 * h1**2 * h2 * (h0 + h2) * H * k_1) / lam

    arr_p_1[1] = (h2 * (h0 + h1)**2 * (h0 + h2) * H * k_1) / lam

    arr_p_1[2] = (((h0 + h1) * (h1 * H + k_1) + h1 * k_1) * (h0 + h2)**2 * k_2) / lam

    arr_p_1[3] = (-1 * ((h0 + h1) * (h1 * H + k_1) + h1 * k_1) * h2**2 * k_2) / lam

    if DEBUG:
        print(f'arr_p_1: {arr_p_1}\n')

    den = (h0 * h2 * (h0 + h2))

    arr_p_2[2] = -1 * (h0 + h2) / den

    arr_p_2[3] = h2 / den

    if DEBUG:
        print(f'arr_p_2: {arr_p_2}\n')

    arr_p = d_2 * (h0 * arr_p_1 / den + arr_p_2)

    if DEBUG:
        print(f'arr_p:   {arr_p}\n')

    return arr_m, arr_p


def calc_interface_temp(e_1: float, t_1: float, e_2: float, t_2: float) -> float:
    """Calculate the temperature at the interface of two semi-infinite bodies.

    Args:
        e_1 (float): Body 1 effusivity [W-sqrt(s)/sq.m-K].
        t_1 (float): Body 1 temperature [degC].
        e_2 (float): Body 2 effusivity [W-sqrt(s)/sq.m-K].
        t_2 (float): Body 2 temperature [degC].

    Returns:
        float: Interface temperature.
    """

    if not isinstance(e_1, NUMERICAL_TYPES) or not isinstance(e_2, NUMERICAL_TYPES):
        raise TypeError('Effusivities e_1 and e_2 must be a numerical type')
    #elif not isinstance(t_1, NUMERICAL_TYPES) or not isinstance(t_2, NUMERICAL_TYPES):
    #    raise TypeError('Temperatures t_1 and t_2 must be a numerical type')
    elif e_1 <= 0 or e_2 <= 0:
        raise ValueError('Effusivities E1 and E2 must be greater than zero.')

    return (e_1 * t_1 + e_2 * t_2) / (e_1 + e_2)


if __name__ == '__main__':
    print(f'\nModule <{Path(__file__).name}> is not intended to be run as standalone module.')
