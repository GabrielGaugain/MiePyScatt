import scipy.special as spe
import numpy as np
from bessel import rj,rh,d_rh, d_rj

def calc_coeff_int(n,k_int, k_ext, a):
      """
      Calculate the scattered field coefficient a_n and b_n 
      Calculation derived from [1] C. F. ohren, D. R. Huffman, Absorption and scattering of light by small particle

      Here, rj (resp rh) is the Riccati Bessel's function of the first kind (resp third kind)
      and d_rj (resp d_rh) is its derivative. These functions are calculated in the bessel.py file.
      
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

      ## Usefull reduced parameters for further calculation
      x = k_ext*a
      m = k_int/k_ext

      ## Calculation of the scattered fierld coefficients for the nth mode (based on [1])
      a_n = ( m* rj(n,m*x) * d_rj(n,x) -    rj(n,x) * d_rj(n,m*x)   )/  \
            ( m* rj(n,m*x) * d_rh(n,x) -  m*rh(n,x) * d_rj(n,m*x)   )
      
      b_n = (  rj(n,m*x) * d_rj(n,x)   -  m* rj(n,x) * d_rj(n,x)    )/  \
            (  rj(n,m*x) * d_rh(n,m*x) -  m* rh(n,x) * d_rj(n,m*x)  )


      return (a_n,b_n)