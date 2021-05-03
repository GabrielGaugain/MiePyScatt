import numpy as np 
import scipy.special as spe
from bessel import hn1, d_rh, d_rj

from scipy.special import lpmn 
#############################################################
###        Script defining spherical vectors M and N      ###
#############################################################
#
#
#   order 1 => j_n + rho = k_1*r (inside)
#   order 3 => h_n + rho = k  *r (outside)
#
#



def M_o1n(n,order,rho,theta,phi, pi_n, tau_n):

    assert (np.size(rho)==np.size(theta))&(np.size(theta)==np.size(phi)), "r, phi and theta should have the same size !!!"

    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
    elif(order==3):
        z_n = hn1(n,rho)
    else :
        assert (order==3)|(order==1), "Order must be 1 (Bessel) or 3 (Hankel)"

    m_o1n = np.zeros((np.size(theta),3 ), dtype=complex)


    m_o1n[:,1] =  np.cos(phi) * pi_n  * z_n       # e_theta

    m_o1n[:,2] = -np.sin(phi) * tau_n * z_n       # e_phi

    return m_o1n



def M_e1n(n, order, rho, theta, phi, pi_n, tau_n):

    assert (np.size(rho)==np.size(theta))&(np.size(theta)==np.size(phi)), \
            "r, phi and theta should have the same size !!!"

    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
    elif(order==3):
        z_n = hn1(n,rho)
    else :
        assert (order==3)|(order==1), "Order must be 1 (Bessel) or 3 (Hankel)"

    m_e1n = np.zeros((np.size(theta),3 ), dtype=complex)

    m_e1n[:,1] = - np.sin(phi) * pi_n  * z_n        # e_theta

    m_e1n[:,2] = - np.cos(phi) * tau_n * z_n        # e_phi


    return m_e1n


def N_o1n(n, order, rho, theta,phi, pi_n, tau_n):

    assert (np.size(rho)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"



    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
        d_zn = d_rj(n,rho)

    elif(order==3):
        z_n = hn1(n, rho)
        d_zn = d_rh(n,rho)
    else :
        assert (order==3)|(order==1), "Order must be 1 (Bessel) or 3 (Hankel)"


    n_o1n = np.zeros((np.size(theta),3 ), dtype=complex)

    n_o1n[:,0] = np.sin(phi)*np.sin(theta)* n*(n+1)* pi_n * z_n / rho   # e_r

    n_o1n[:,1] = np.sin(phi) * tau_n * d_zn /rho                         # e_theta

    n_o1n[:,2] = np.cos(phi) * pi_n *  d_zn/ rho                         # e_phi

    return n_o1n



def N_e1n(n, order, rho, theta,phi, pi_n, tau_n):

    assert (np.size(rho)==np.size(theta))&(np.size(theta)==np.size(phi)) , "r, phi and theta should have the same size !!!"

    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
        d_zn = d_rj(n,rho)

    elif(order==3):
        z_n = hn1(n, rho)
        d_zn = d_rh(n,rho)
    else :
        assert (order==3)|(order==1), "Order must be 1 (Bessel) or 3 (Hankel)"

    n_e1n = np.zeros((np.size(theta),3 ), dtype=complex)

    n_e1n[:,0] = np.cos(phi)*np.sin(theta)* n*(n+1)* pi_n * z_n / rho   # e_r

    n_e1n[:,1] = np.cos(phi) * tau_n * d_zn /rho                         # e_theta

    n_e1n[:,2] =-np.sin(phi) * pi_n  * d_zn/ rho                         # e_phi

    return n_e1n




def angle_functions(N, mu):
    """
    Recoded angle functions since lpmn function of the scipy special package is not vectorized
    """

    print( "In angle function ...")
    pi_n ,tau_n = np.zeros( (mu.size, N+1) ) ,  np.zeros( (mu.size, N+1) ) 

    pi_n[:,0] = 0    
    pi_n[:,1] = 1

    tau_n[:,0] = 0 
    tau_n[:,1] = mu

    for n in range(2,N+1):
        
        pi_n[:,n]  = (2*n -1)/(n-1) * mu * pi_n[:,n-1]   -  n/(n-1) * pi_n[:,n-2]
        
        tau_n[:,n] = n * mu * pi_n[:,n] - (n+1) * pi_n[:,n-1]

    return pi_n,tau_n


if __name__ == "__main__":
    

    print(M_o1n(1,1, 0.2, 10,10,1,1))