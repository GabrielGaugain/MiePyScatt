import numpy as np
import matplotlib.pyplot as plt
import Mie_field_calculation as mie
import Gabriel_eq as diel
import scipy.constants as ctes
import coordinates as coor

def main():

    ## frequence(s) of interest
    f = 10**6
    lmda = ctes.c / f

    n_ext = 1                   # outer refrective index : 1 for the air
    k_ext = 2*np.pi/lmda * n_ext 


    ## medium of the sphere
    tissue = "Muscle"
    perm, cond = diel.get_diel_properties(tissue, freqs=f)

    n_int = np.sqrt(perm)   ## without loss => no cmplx part contrib
    k_int = 2*np.pi/lmda * n_int

    
    #####################################
    a = 0.1                             # 20 cm , radius of the sphere
    L_x, L_y, L_z = 0.4, 0.4, 0.4       # size of the grid on which to compute the fields
    N_x, N_y, N_z = 400,400,400         # number of points in each direction of the grid 

    ## grid (in cartesian coordinates)
    X = np.meshgrid( np.linspace(0.001,L_x,N_x) - L_x/2, \
                     np.linspace(0.001,L_y,N_y) - L_y/2, \
                     np.linspace(0.001,L_z,N_z) - L_z/2 )
    x,y,z = X[0].ravel(), X[1].ravel(), X[2].ravel()
   

    
    ## converting the grid into spherical coordinates for calculations
    (r, theta, phi) = coor.cart2sph(x,y,z)


    r_int, theta_int, phi_int = r[r<a], theta[r<a] , phi[r<a]    
    r_ext, theta_ext, phi_ext = r[r>=a], theta[r>=a] , phi[r>=a]

    print("sca field calc")
    #E_int = mie.E_internal(r_int,theta_int,phi_int, k_int,k_ext,a)
    E_s   = mie.E_scattered(r_ext,theta_ext,phi_ext, k_int,k_ext,a)

    #E_external = E_s + np.exp(2*np.pi * k_ext * r_ext)/r_ext 
    
    


    ##end main





if __name__ == "__main__":
    main()