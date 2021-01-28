import numpy as np 
import scipy.special as spe
from bessel import hn1, d_rh, d_rj

import scipy.special.lpmn as Pmn
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

    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
    elif(order==3):
        z_n = hn1(n,rho)
    else :
        print("you have to enter order 1 or 3, "+order+" is not accepted")

    m_o1n = np.zeros((np.size(theta),3 ))


    m_o1n[:,2] = np.cos(phi) * pi_n  * zn       # e_theta

    m_01n[:,3] = np.sin(phi) * tau_n * zn       # e_phi

    return m_01n



def M_e1n(n, rho, theta, phi, pi_n, tau_n):

    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
    elif(order==3):
        z_n = hn1(n,rho)
    else :
        print("you have to enter order 1 or 3, "+order+" is not accepted")

    m_o1n = np.zeros((np.size(theta),3 ))


    m_o1n[:,2] = - np.sin(phi) * pi_n  * zn         # e_theta

    m_01n[:,3] = - np.cos(phi) * tau_n * zn         # e_phi


    return


def N_o1n(n,theta,phi, pi_n, tau_n):


    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
        d_r = d_rj(n,rho)

    elif(order==3):
        z_n = hn1(n, rho)
        d_r = d_rh(n,rho)
    else :
        print("you have to enter order 1 or 3, "+order+" is not accepted")

    m_o1n = np.zeros((np.size(theta),3 ))


    m_o1n[:,1] = np.sin(phi)* n*(n+1)*sin(theta)* pi_n * zn / rho   # e_r

    m_01n[:,2] = np.sin(phi) * tau_n * d_r /rho                     # e_theta

    m_01n[:,3] = np.cos(phi) * tau_n * d_r/ rho                     # e_phi

    return



def N_e1n(n,theta,phi, pi_n, tau_n):


    if (order ==1 ):
        z_n = spe.spherical_jn(n,rho)
        d_r = d_rj(n,rho)

    elif(order==3):
        z_n = hn1(n, rho)
        d_r = d_rh(n,rho)
    else :
        print("you have to enter order 1 or 3, "+order+" is not accepted")

    m_o1n = np.zeros((np.size(theta),3 ))


    m_o1n[:,1] = np.cos(phi)* n*(n+1)*sin(theta)* pi_n * zn / rho   # e_r

    m_01n[:,2] = np.cos(phi) * tau_n * d_r /rho                     # e_theta

    m_01n[:,3] = np.sin(phi) * tau_n * d_r/ rho                     # e_phi

    return



