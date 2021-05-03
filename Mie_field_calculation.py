import numpy as np
import scipy.special as spe

from coeff_int import calc_coeff_int
from coeff_sca import calc_coeff_scat
from spherical_vectors import M_o1n, M_e1n, N_e1n, N_o1n
import coordinates as coor

#############################################################
###        Script for field calculations                  ###
#############################################################

N_max = 10

def E_internal(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=N_max):
    """
    Computation of the field in a dielectric sphere for an incident OPP.
    The calculations are based on [1] and use spherical harmonics decomposition.
    
    Params :
    --------
                r, theta, phi : array[float]  | coordinates of the grid on which the field is computed
                k_int : real | wave vector of the inner medium
                k_ext : real | wave vector of the outer medium
                a     : real | radius of the sphere
                E_0   : real | amplitude of the incident field (default : 1 V/m )
                N     : int  | number of iterations/harmonics computed (default : N_max)

    Returns :
    ---------
                E_i :   array[complex], Nx3 | Internal field at each point of the grid, in each direction (spherical coordinates)

    """


    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"

    pi_n, tau_n = angle_functions(N, np.cos(theta))

    #pi_n = pi_n/np.sin(theta)     

    E_i = np.zeros( (np.size(r), 3), dtype=complex )

    for n in range(1,N+1):

        E_n = E_0* 1j**n * (2*n + 1)/( n * (n+1) )

        c_n, d_n = calc_coeff_int(n, k_int, k_ext, a)

        E_i = E_i + E_n* (     c_n * M_o1n(n, 1, k_int*r, theta, phi, pi_n[:,n-1], tau_n[:,n-1])  \
                  -        1j* d_n * N_e1n(n, 1, k_int*r, theta, phi, pi_n[:,n-1], tau_n[:,n-1]) )

        ##end for 

    return E_i


def E_scattered(r, theta, phi, k_int, k_ext, a, E_0 = 1, N=N_max):
    """
    Computation of the scattered field by a dielectric sphere for an incident OPP.
    The calculations are based on [1] and use spherical harmonics decomposition.
    
    Params :
    --------
                r, theta, phi : array[float]  | coordinates of the grid on which the field is computed
                k_int : real | wave vector of the inner medium
                k_ext : real | wave vector of the outer medium
                a     : real | radius of the sphere
                E_0   : real | amplitude of the incident field (default : 1 V/m )
                N     : int  | number of iterations/harmonics computed (default : N_max)

    Returns :
    ---------
                E_s :   array[complex], Nx3 | scattered field at each point of the grid, in each direction (spherical coordinates)

    """

    assert (np.size(r)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"
   
    pi_n, tau_n = angle_functions(N, np.cos(theta))

    E_s = np.zeros( (np.size(r),3), dtype=complex )

    for n in range(1,N+1):

        E_n = E_0* 1j**(n) * (2*n + 1)/( n * (n+1) )

        a_n, b_n = calc_coeff_scat(n, k_int, k_ext, a)

        E_s = E_s  + E_n * ( 1j* a_n * N_e1n(n, 3, k_ext*r, theta, phi, pi_n[:,n-1], tau_n[:,n-1])    \
                   -             b_n * M_o1n(n, 3, k_ext*r, theta, phi, pi_n[:,n-1], tau_n[:,n-1])  )

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




def calc_E_field(X, k_int, k_ext, a, E_0=1, N=N_max):
    """
    Computation of the total electric field for the particular case of a sphere in an electric field (incident OPP).
    The calculations are based on [1] and use spherical harmonics decomposition.
    
    Params :
    --------
                X     : grid | mesh coordinates X = [x, y, z] in cartesian coordinates
                k_int : real | wave vector of the inner medium
                k_ext : real | wave vector of the outer medium
                a     : real | radius of the sphere
                E_0   : real | amplitude of the incident field (default : 1 V/m )
                N     : int  | number of iterations/harmonics computed (default : N_max)

    Returns :
    ---------
                E     :   array[complex], Nx3 | Electric field at each point of the grid, in each direction in cartesian coordinates

    """    
    (r, theta, phi) = coor.cart2sph(X[0],X[1],X[2])

    E = np.zeros((np.size(r),3), dtype=complex)

    r_int, theta_int, phi_int = r[r<a], theta[r<a] , phi[r<a]    
    r_ext, theta_ext, phi_ext = r[r>=a], theta[r>=a] , phi[r>=a]

    print("Fields calculations ...")
    E_int = E_internal(r_int,theta_int,phi_int, k_int,k_ext,a)
    E_sca = E_scattered(r_ext,theta_ext,phi_ext, k_int,k_ext,a)
    print("Fields computed !")

    ## transform back the fields into cartesian coordinates
    
    
    for i in range(theta_int.size):
        R= coor.Mat_sph2cart(theta_int[i],phi_int[i])
        E_int[i,:] = R@E_int[i,:]

    for i in range(theta_ext.size):
        R= coor.Mat_sph2cart(theta_ext[i],phi_ext[i])
        E_sca[i,:] = R@E_sca[i,:]


    """    
    R = coor.Mat_sph2cart(theta,phi)

    print("E size : ", E_int.shape )
    print("R size : ", R.shape)

    E_int = R[r<a,:,:]@np.expand_dims(E_int, axis=0)
    print("E size : ", E_int.shape )

    E_int = E_int[:,0,:]

    E_sca = R[r>=a,:,:]@np.expand_dims(E_sca, axis=0)
    E_sca = E_sca[:,0,:]
    """
    ## total field => + OPP
    E[r<a,:] = E_int
    E[r>=a,:] = E_sca 
    E[r>=a,0] =  E[r>=a,0] + np.exp(1j*k_ext * r[r>=a]*np.cos(theta[r>=a]))

    return E

