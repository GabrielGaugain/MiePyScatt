import scipy.special as spe
import numpy as np

def calc_coeff_int(k_int, k_ext, a):
    """
    Calculate the internal field coefficient c_n and d_n (from Absorption an scattering of light by small particle, C. F. ohren, D. R. Huffman)

    Params :
    --------
                n     : int  | mode
                k_int : real | wave vector of the internal medium
                k_ext : real | wave vector of the external medium
                a     : real | radius of the sphere


    Returns :
    ---------

                c_n :   real | coefficient for the ...
                d_n :   real | coefficient for the ...


    """

    ## usefull parameters for calculations
    x = k_int*a
    m = k_int/k_ext

    ## calculation of the coefficient from boundary conditions (see [1])
    c_n =  ( m*rj(x)*d_rh(x)  -  m*rh(x)*d_rj(x)   )/ \
           ( rj(m*x)*d_rh(x)  -  m*rh(x)*d_rj(m*x) ) 


    d_n = ( m*rj(x)*d_rh(x)   -  m*rh(x)*d_rj(x) )/ \
          ( m*rj(m*x)*d_rh(x) -  rh(x)d_rj(m*x)  )


    return (c_n, d_n)