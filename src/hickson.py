"""This module contains functions from Hickson et al., 2011."""

import numpy as np

FLOAT64 = np.float64

def conductivity_match_coeff(D1: float, K1: float, D2: float, K2: float, dZ: float) -> np.ndarray:
    """Calculate the second order finite difference (FD) approximation coefficients at the two-body interface using
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