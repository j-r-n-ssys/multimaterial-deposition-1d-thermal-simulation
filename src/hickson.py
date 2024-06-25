"""This module contains functions from Hickson et al., 2011."""

import numpy as np

FLOAT64 = np.float64

DEBUG = True


def conductivity_match_coeff(D1: float, K1: float, D2: float, K2: float, dZ: float) -> np.ndarray:
    """From 'Finite Difference Schemes for Multilayer Diffusion' by Hickson et al., 2011. 
    Calculate the second order finite difference (FD) approximation coefficients at the two-body interface using
    conductivity matching as derived in in 'Finite Difference Schemes for Multilayer Diffusion' by Hickson et al.
    A relevant assumption is that the interface is at a node. 

    Args:
        D1 (float): Body 1 thermal diffusivity. 
        K1 (float): Body 1 thermal conductivity. 
        D2 (float): Body 2 thermal diffusivity. 
        K2 (float): Body 2 thermal conductivity. 
        dZ (float): Finite difference node spacing. 

    Returns:
        np.ndarray: Finite difference (FD) approximation coefficients.
    """

    dZ2 = dZ**2

    Yp2 = (2 * K1 + 3 * K2) * D1 - D2 * K1

    Yp1 = 4 * D2 * K1 - 2 * (K1 + 3 * K2) * D1

    Ym1 = 4 * D1 * K2 - 2 * (3 * K1 + K2) * D2

    Ym2 = (3 * K1 + 2 * K2) * D2 - D1 * K2

    den = (6 * (K1 + K2) * dZ2)

    coeff = np.array([Yp2, Yp1, 0, Ym1, Ym2], dtype=FLOAT64) / den  # finite difference approximation

    return coeff


def jump_match_coeff(K1, K2, h0, h1, h2, H) -> tuple[np.ndarray, np.ndarray]:
    """From 'Finite Difference Schemes for Multilayer Diffusion' by Hickson et al., 2011. 

    Args:
        K1 (_type_): _description_
        K2 (_type_): _description_
        h0 (_type_): _description_
        h1 (_type_): _description_
        h2 (_type_): _description_
        H (_type_): _description_

    Returns:
        tuple[np.ndarray, np.ndarray]: _description_
    """

    #   0       -2
    #   1       -1
    #   2        1
    #   3        2
    #
    #   h1
    #   D_i         K_i         K1
    #
    #   h2
    #   D_i+1       K_i+1       K2
    #
    #

    print(' ')

    if not isinstance(K1, (int, float)) or not isinstance(K1, (int, float)):
        raise TypeError('Conducitivity K1 and K2 must be a numerical type.')

    if K1 <= 0 or K2 <= 0:
        raise ValueError('Conducitivty K1 and K2 must be greater than zero.')

    if not isinstance(h0, (int, float)) or not isinstance(h1, (int, float)) or not isinstance(h2, (int, float)):
        raise TypeError('Node spacing h0, h1, and h2 must be a numerical type.')

    if h0 <= 0:
        raise ValueError('Node-to-node distance h0 must be greater than zero.')

    if h1 <= 0 or h2 <= 0:
        raise ValueError('Node-to-interface distances h1 and h2 must be greater than zero.')

    if h1 + h2 != h0:
        raise ValueError('The sume of node-to-interface distances h1 and h2 must be equal node-to-node distance h0')

    if not isinstance(K1, (int, float)) or not isinstance(K1, (int, float)):
        raise TypeError('Conducitivity K1 and K2 must be a numerical type.')

    arr_m = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_2 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_2 = np.zeros(shape=4, dtype=FLOAT64)

    lam = h0 * (((h0 + 2 * h1) * (h0 + h2) * h2 * K1 + (h0 + h1) * (h0 + 2 * h2) * h1 * K2) * H + (h0 + 2 * h1) *
                (h0 + 2 * h2) * K1 * K2)

    arr_m_1[0] = (-1 * ((h0 + h2) * (h2 * H + K2) + h2 * K2) * h1**2 * K1) / lam

    arr_m_1[1] = (((h0 + h2) * (h2 * H + K2) + h2 * K2) * (h0 + h1)**2 * K1) / lam

    arr_m_1[2] = (h1 * (h0 + h1) * (h0 + h2)**2 * H * K2) / lam

    arr_m_1[3] = (-1 * h1 * h2**2 * (h0 + h1) * H * K2) / lam

    if DEBUG:
        print(f'arr_m_1: {arr_m_1}')

    den = (h0 * h1 * (h0 + h1))

    arr_m_2[0] = h1 / den

    arr_m_2[1] = -1 * (h0 + h1) / den

    if DEBUG:
        print(f'arr_m_2: {arr_m_2}')

    arr_m = h0 * arr_m_1 / den + arr_m_2

    if DEBUG:
        print(f'arr_m:   {arr_m}')

    arr_p_1[0] = (-1 * h1**2 * h2 * (h0 + h2) * H * K1) / lam

    arr_p_1[1] = (h2 * (h0 + h1)**2 * (h0 + h2) * H * K1) / lam

    arr_p_1[2] = (((h0 + h1) * (h1 * H + K1) + h1 * K1) * (h0 + h2)**2 * K2) / lam

    arr_p_1[3] = (-1 * ((h0 + h1) * (h1 * H + K1) + h1 * K1) * h2**2 * K2) / lam

    if DEBUG:
        print(f'arr_p_1: {arr_p_1}')

    den = (h0 * h2 * (h0 + h2))

    arr_p_2[2] = -1 * (h0 + h2) / den

    arr_p_2[3] = h2 / den

    if DEBUG:
        print(f'arr_p_2: {arr_p_2}')

    arr_p = h0 * arr_p_1 / den + arr_p_2

    if DEBUG:
        print(f'arr_p:   {arr_p}')


if __name__ == '__main__':

    jump_match_coeff(4, 4, 2, 1, 1, 10**15)
