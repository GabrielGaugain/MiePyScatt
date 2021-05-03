import scipy.special as spe
import numpy as np
from bessel import rj,rh,d_rh, d_rj


##################################################################################
###       Script for the caclulation of coefficient for the external field     ###
##################################################################################



def scat_coeff(n, x, m):
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

      ## Calculation of the scattered fierld coefficients for the nth mode (based on [1])
      a_n = ( m* rj(n,m*x) * d_rj(n,x) -    rj(n,x) * d_rj(n,m*x)   )/  \
            ( m* rj(n,m*x) * d_rh(n,x) -    rh(n,x) * d_rj(n,m*x)   )
      
      b_n = (  rj(n,m*x) * d_rj(n,x) -  m* rj(n,x) * d_rj(n,m*x)    )/  \
            (  rj(n,m*x) * d_rh(n,x) -  m* rh(n,x) * d_rj(n,m*x)  )


      return (a_n,b_n)






##################################################################################
###       Script for the caclulation of coefficient for the internal field     ###
##################################################################################

def int_coeff(n, x, m):
      """
      Calculate the internal field coefficient c_n and d_n 
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

                  c_n :   real | coefficient for the ...
                  d_n :   real | coefficient for the ...


      """


      ## calculation of the coefficient from boundary conditions (see [1])
      c_n = ( m*rj(n,x) *  d_rh(n,x)  -  m* rh(n,x) * d_rj(n,x)    )/  \
            (   rj(n,m*x)* d_rh(n,x)  -  m* rh(n,x) * d_rj(n,m*x)  ) 


      d_n = ( m*   rj(n,x )  * d_rh(n,x) - m* rh(n,x) * d_rj(n,x)    )/  \
            (m**2* rj(n,m*x) * d_rh(n,x) -    rh(n,x) * d_rj(n,m*x)  )


      return (c_n, d_n)