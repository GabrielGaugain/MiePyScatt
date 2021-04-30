import numpy as np
import matplotlib.pyplot as plt
import Mie_field_calculation as mie
import Gabriel_eq as diel
import scipy.constants as ctes
import coordinates as coor

def main():

    ## frequency(ies) of interest
    f = 0.39035*10**9
    lmda = ctes.c / f

    n_ext = 1                   # outer refrective index : 1 for the air
    k_ext = 2*np.pi/lmda * n_ext 

    ## medium of the sphere
    tissue = "Muscle"
    perm, cond = diel.get_diel_properties(tissue, freqs=f)

    n_int = np.sqrt(perm)   ## without loss => no cmplx part contrib
    print("n_int : " ,n_int)
    
    k_int = 2*np.pi/lmda * n_int
    k_int = k_int[0,0]
    print("kint : ", k_int)

    #####################################
    a = 0.1                             # 20 cm , radius of the sphere
    L_x, L_y, L_z = 0.4, 0.4, 0.4       # size of the grid on which to compute the fields
    N_x, N_y, N_z = 150,150,150         # number of points in each direction of the grid 

    ## grid (in cartesian coordinates)
    ## 3D too computationnaly expensive 
    #X = np.meshgrid( np.linspace(0.001,L_x,N_x) - L_x/2, \
    #                 np.linspace(0.001,L_y,N_y) - L_y/2, \
    #                 np.linspace(0.001,L_z,N_z) - L_z/2 )
    #x,y,z = X[0].ravel(), X[1].ravel(), X[2].ravel()
   
    ## 2D, in the plane z=0
    X = np.meshgrid( np.linspace(0,L_x,N_x) - L_x/2, \
                     np.linspace(0,L_y,N_y) - L_y/2, \
                    )
    x,y = X[0].ravel(), X[1].ravel()
    z = np.zeros(x.shape)
    

    E = mie.calc_E_field([x,y,z], k_int, k_ext, a)    
    
    Ex = E[:,0].reshape(X[0].shape)
    Ey = E[:,1].reshape(X[0].shape)
    Ez = E[:,2].reshape(X[0].shape)
    
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True )
    axs[0].pcolormesh(X[0],X[1], np.abs(Ex))
    axs[1].pcolormesh(X[0],X[1], np.abs(Ey))
    axs[2].pcolormesh(X[0],X[1], np.abs(Ez))
    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")

    plt.show()    

    
    ##end main





if __name__ == "__main__":
    main()