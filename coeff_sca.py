import scipy.special as spe
import numpy as np
from .bessel import *

def calc_coeff_int(n,k_int, k_ext, a):
    """
    Calculate the scattered field coefficient a_n and b_n (from Absorption an scattering of light by small particle, C. F. ohren, D. R. Huffman)

    Params :
    --------
                n     : int  | mode
                k_int : real | wave vector of the internal medium
                k_ext : real | wave vector of the external medium
                a     : real | radius of the sphere


    Returns :
    ---------

                a_n :   complex | coefficient for the ...
                b_n :   complex | coefficient for the ...


    """

    x = k_ext*a
    m = k_int/k_ext

    a_n = ( m*rj(m*x)*d_rj(x)  - rj(x)*d_rj(m*x) )/...
          ( m*rj(m*x)*d_rh(x) - m*rh(x)*d_rj(m*x) )
    
    b_n = ( rj(m*x)*d_rj(x) - m*rj(x)*d_rj(x) )/...
          ( rj(m*x)*d_rh(m*x) - m*d_rj(m*x) )


    return 