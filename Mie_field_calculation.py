import numpy as np
import scipy.special as spe

from coeff_int import calc_coeff_int
from coeff_sca import calc_coeff_scat
from spherical_vectors import M_o1n, M_e1n, N_e1n, N_o1n


#############################################################
###        Script for field calculations                  ###
#############################################################

N_max = 10

def E_internal(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=N_max):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)

    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"

    pi_n, tau_n = spe.lpmn(1,N, np.cos(theta))

    #pi_n = pi_n/np.sin(theta)     

    E_i = np.zeros( (np.size(r), 3) )

    for n in range(1,N+1):

        E_n = E_0* np.i**(n) * (2*n + 1)/( n * (n+1) )

        c_n, d_n = calc_coeff_int(n, k_int, k_ext, a)

        E_i = E_i + E_n* (       c_n * M_o1n(n, 1, k_int*r, theta, phi, pi_n, tau_n)  \
                  -        np.i* d_n * N_e1n(n, 1, k_int*r, theta, phi, pi_n, tau_n) )

        ##end for 

    return E_i


def E_scattered(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=N_max):

    ## calculation of the associated legendre polynomials of the first order
    ## from degre 0 to N, and their first derivative (tau_n)

    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"
   
    pi_n, tau_n = angle_functions(N, np.cos(theta))

    print("angle function - OK ")

    E_s = np.zeros( (np.size(r),3) )

    for n in range(1,N+1):

        E_n = E_0* np.i**(n) * (2*n + 1)/( n * (n+1) )

        a_n, b_n = calc_coeff_scat(n, k_int, k_ext, a)

        E_s = E_s  + E_n * ( np.i* a_n * N_e1n(n, 3, k_ext*r, theta, phi, pi_n[:,n-1], tau_n[:,n])    \
                   -               b_n * M_o1n(n, 3, k_ext*r, theta, phi, pi_n[:,n-1], tau_n[:,n])  )

        ##end for 

    return E_s




def angle_functions(N, mu):
    """
    Recoded angle functions since lpmn function of the scipy special package is not vectorized
    """

    print( "In angle function ...")
    pi_n ,tau_n = np.zeros( (mu.size, N) ) ,  np.zeros( (mu.size, N) ) 

    pi_n[:,0] = 1    
    pi_n[:,1] = 3*mu

    tau_n[:,0] = mu 
    tau_n[:,1] = 6*mu**2 -3

    for n in range(2,N):
        m = n+1

        pi_n[:,n] = ( 2*m -1)/(m-1)* mu*pi_n[:,n-1] - m/(m-1)*pi_n[:,n-2]
        
        tau_n[:,n] = m*mu*pi_n[:,n] - (m+1)*pi_n[:,n-1]

    return pi_n,tau_n



