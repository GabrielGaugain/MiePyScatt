import numpy as np
import scipy.special as spe

from coeff_int import calc_coeff_int
from coeff_sca import calc_coeff_scat
from spherical_vectors import M_o1n, M_e1n, N_e1n, N_o1n


#############################################################
###        Script for field calculations                  ###
#############################################################





def E_internal(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=10):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)

    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"

    pi_n, tau_n = spe.lpmn(1,N, np.cos(theta))

    pi_n = pi_n/np.sin(theta)     


    E_i = np.zeros( (np.size(r), 3) )

    for n in range(1,N+1):

        E_n = E_0* np.i**(n) * (2*n + 1)/( n * (n+1) )

        c_n, d_n = calc_coeff_int(n, k_int, k_ext, a)

        E_i = E_i + E_n* ( c_n * M_o1n(n, 1, k_int*r, theta, phi, pi_n, tau_n)  \
                  +        d_n * N_e1n(n, 1, k_int*r, theta, phi, pi_n, tau_n) )

        ##end for 

    return E_i


def E_scattered(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=10):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)

    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"
   
    pi_n, tau_n = spe.lpmn(1,N, np.cos(theta))
    
    pi_n = pi_n/np.sin(theta)     


    E_s = np.zeros( (np.size(r),3) )

    for n in range(1,N+1):

        E_n = E_0* np.i**(n) * (2*n + 1)/( n * (n+1) )

        a_n, b_n = calc_coeff_scat(n, k_int, k_ext, a)

        E_s = E_s  + E_n * ( a_n * M_o1n(n, 3, k_int*r, theta, phi, pi_n, tau_n)    \
                   +         b_n * N_e1n(n, 3, k_int*r, theta, phi, pi_n, tau_n)  )

        ##end for 

    return E_s








