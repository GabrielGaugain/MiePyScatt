import numpy as np
import scipy.special as spe

from coeff_int import calc_coeff_int
from coeff_sca import calc_coeff_scat


#############################################################
###        Script for field calculations                  ###
#############################################################





def E_internal(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=10):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)
    pi_n, tau_n = spe.lpmn(1,N, cos(theta))

    pi_n = pi_n/np.sin(theta)     


    E_i = np.zeros( np.size(r), np.size(theta), np.size(phi) )

    for i in range(1,N+1):

        E_n = E_0* np.i**(n) * (2*n + 1)/( n * (n+1) )

        a_n, b_n = calc_coeff_int(n, k_int, k_ext, a)

        

        ##end for 

    return


def E_scattered(r, theta, phi, N):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)
    pi_n, tau_n = spe.lpmn(1,N, cos(theta))
    
    pi_n = pi_n/np.sin(theta)     




    return








