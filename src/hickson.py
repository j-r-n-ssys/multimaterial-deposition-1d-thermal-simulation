"""R. I. Hickson, S. I. Barry, G. N. Mercer, and H. S. Sidhu, “Finite 
difference schemes for multilayer diffusion,” Mathematical and computer 
modelling, vol. 54, no. 1-2, pp. 210-220, Jul. 2011, doi: 
https://doi.org/10.1016/j.mcm.2011.02.003."""

import numpy as np

FLOAT64 = np.float64

DEBUG = False


def cond_match_coeff(K1: float, D1: float, K2: float, D2: float, h: float) -> np.ndarray:
    """Calculate the second order finite difference (FD) approximation 
    coefficients at the two-body interface using conductivity matching. A 
    fundamental assumption is that the interface is co-located with the middle
    node of this FD approximation. 

    Args:
        K1 (float): Body 1 thermal conductivity [W/m-K]. 
        D1 (float): Body 1 thermal diffusivity [sq.m/s].
        K2 (float): Body 2 thermal conductivity [W/m-K]. 
        D2 (float): Body 2 thermal diffusivity [sq.m/s]. 
        dZ (float): Finite difference node-to-node distance [m].

    Returns:
        np.ndarray: Finite difference (FD) approximation coefficients.
    """

    if not isinstance(K1, (int, float)) or not isinstance(K1, (int, float)):
        raise TypeError('Conducitivities K1 and K2 must be a numerical type.')

    if K1 <= 0 or K2 <= 0:
        raise ValueError('Conducitivties K1 and K2 must be greater than zero.')

    if not isinstance(D1, (int, float)) or not isinstance(D2, (int, float)):
        raise TypeError('Diffusivities D1 and D2 must be a numerical type.')

    if D1 <= 0 or D2 <= 0:
        raise ValueError('Diffusivities D1 and D2 must be greater than zero.')

    if not isinstance(h, (int, float)):
        raise TypeError('Node-to-node distance h must be a numerical type.')

    if h <= 0:
        raise ValueError('Node-to-node distance h must be greater than zero.')

    # Initialize the 5 x 1 coefficient matrix.
    coeff = np.zeros(shape=5, dtype=FLOAT64)

    # Calculate
    den = (6 * (K1 + K2) * h**2)

    print(den)

    coeff[0] = ((2 * K1 + 3 * K2) * D1 - D2 * K1) / den

    coeff[1] = (4 * D2 * K1 - 2 * (K1 + 3 * K2) * D1) / den

    coeff[3] = (4 * D1 * K2 - 2 * (3 * K1 + K2) * D2) / den

    coeff[4] = ((3 * K1 + 2 * K2) * D2 - D1 * K2) / den

    if DEBUG:
        print(coeff)

    return coeff


def jump_match_coeff(K1: float, D1: float, K2: float, D2: float, h0: float, h1: float, h2: float,
                     H: float) -> tuple[np.ndarray, np.ndarray]:
    """Calculate the second order finite difference (FD) approximation 
    coefficients on both sides of a two-body interface using jump matching. A 
    fundamental assumption is that the interface is not co-located with a node.

    Args:
        K1 (float): Body 1 thermal conductivity [W/m-K]. 
        D1 (float): Body 1 thermal diffusivity [sq.m/s].
        K2 (float): Body 2 thermal conductivity [W/m-K]. 
        D2 (float): Body 2 thermal diffusivity [sq.m/s]. 
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

    if not isinstance(K1, (int, float)) or not isinstance(K1, (int, float)):
        raise TypeError('Conducitivities K1 and K2 must be a numerical type.')

    if K1 <= 0 or K2 <= 0:
        raise ValueError('Conducitivties K1 and K2 must be greater than zero.')

    if not isinstance(D1, (int, float)) or not isinstance(D2, (int, float)):
        raise TypeError('Diffusivities D1 and D2 must be a numerical type.')

    if D1 <= 0 or D2 <= 0:
        raise ValueError('Diffusivities D1 and D2 must be greater than zero.')

    if not isinstance(h0, (int, float)):
        raise TypeError('Node-to-node distnce h0 must be a numerical type.')

    if h0 <= 0:
        raise ValueError('Node-to-node distance h0 must be greater than zero.')

    if not isinstance(h1, (int, float)) or not isinstance(h2, (int, float)):
        raise TypeError('Node-to-interface distance h1 and h2 must be a numerical type.')

    if h1 <= 0 or h2 <= 0:
        raise ValueError('Node-to-interface distances h1 and h2 must be greater than zero.')

    if (h1 + h2) != h0:
        raise ValueError('The sume of node-to-interface distances h1 and h2 must be equal node-to-node distance h0')

    if not isinstance(H, (int, float)):
        raise TypeError('Contact transfer coefficient H must be a numerical type.')

    if H < 0:
        raise TypeError('Contact transfer coefficient H must be greater than or equal to zero.')

    arr_m = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_m_2 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_1 = np.zeros(shape=4, dtype=FLOAT64)

    arr_p_2 = np.zeros(shape=4, dtype=FLOAT64)

    lam = h0 * (((h0 + 2 * h1) * (h0 + h2) * h2 * K1 + (h0 + h1) * (h0 + 2 * h2) * h1 * K2) * H + (h0 + 2 * h1) *
                (h0 + 2 * h2) * K1 * K2)

    #      0       -2        1
    #      1       -1        1
    #      2        1        2
    #      3        2        2

    arr_m_1[0] = (-1 * ((h0 + h2) * (h2 * H + K2) + h2 * K2) * h1**2 * K1) / lam

    arr_m_1[1] = (((h0 + h2) * (h2 * H + K2) + h2 * K2) * (h0 + h1)**2 * K1) / lam

    arr_m_1[2] = (h1 * (h0 + h1) * (h0 + h2)**2 * H * K2) / lam

    arr_m_1[3] = (-1 * h1 * h2**2 * (h0 + h1) * H * K2) / lam

    if DEBUG:
        print(f'arr_m_1: {arr_m_1}\n')

    den = (h0 * h1 * (h0 + h1))

    arr_m_2[0] = h1 / den

    arr_m_2[1] = -1 * (h0 + h1) / den

    if DEBUG:
        print(f'arr_m_2: {arr_m_2}\n')

    arr_m = D1 * (h0 * arr_m_1 / den + arr_m_2)

    if DEBUG:
        print(f'arr_m:   {arr_m}\n')

    #      0       -2        1
    #      1       -1        1
    #      2        1        2
    #      3        2        2

    arr_p_1[0] = (-1 * h1**2 * h2 * (h0 + h2) * H * K1) / lam

    arr_p_1[1] = (h2 * (h0 + h1)**2 * (h0 + h2) * H * K1) / lam

    arr_p_1[2] = (((h0 + h1) * (h1 * H + K1) + h1 * K1) * (h0 + h2)**2 * K2) / lam

    arr_p_1[3] = (-1 * ((h0 + h1) * (h1 * H + K1) + h1 * K1) * h2**2 * K2) / lam

    if DEBUG:
        print(f'arr_p_1: {arr_p_1}\n')

    den = (h0 * h2 * (h0 + h2))

    arr_p_2[2] = -1 * (h0 + h2) / den

    arr_p_2[3] = h2 / den

    if DEBUG:
        print(f'arr_p_2: {arr_p_2}\n')

    arr_p = D2 * (h0 * arr_p_1 / den + arr_p_2)

    if DEBUG:
        print(f'arr_p:   {arr_p}\n')

    return arr_m, arr_p


if __name__ == '__main__':

    jump_match_coeff(4, 1, 4, 1, 2, 1, 1, 0)

def calc_interface_temp(E1:float, T1:float, E2:float, T2:float) -> float:
    """Calculate the temperature at the interface of two semi-infinite 

    Args:
        E1 (float): Body 1 effusivity [W-sqrt(s)/sq.m-K].
        T1 (float): Body 1 temperature [degC].
        E2 (float): Body 2 effusivity [W-sqrt(s)/sq.m-K].
        T2 (float): Body 2 temperature [degC].

    Raises:
        TypeError: _description_
        ValueError: _description_

    Returns:
        float: _description_
    """

    if not isinstance(E1, (int,float)) or not isinstance(E1, (int,float)) or not isinstance(E1, (int,float)) or not isinstance(E1, (int,float)):
        raise TypeError('Temperatures T1 and T2 and effusivities E1 and E2 must be a numerical type')
    
    if E1 <= 0 or E2 <= 0:
        raise ValueError('Effusivities E1 and E2 must be greater than zero.')
    
    return (E1 * T1 + E2 * T2) /(E1 + E2)