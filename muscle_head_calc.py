import numpy as np
import matplotlib.pyplot as plt
import mie_field_calculation as mie
import Gabriel_eq as diel
import scipy.constants as ctes
import scipy.io as io
import coordinates as coor
import quasi_static as qs


def main():

    ## frequency(ies) of interest
    f = 10**6
    lmda = ctes.c / f

    n_ext = 1.0                   # outer refrective index : 1 for the air
    k_ext = 2*np.pi/lmda * n_ext 

    ## medium of the sphere
    tissue = "Muscle"
    perm, cond = diel.get_diel_properties(tissue, freqs=f)
    print(perm, cond)

    n_int = np.sqrt(perm)   ## without loss => no cmplx part contrib
    print("n_int : " ,n_int)
    
    k_int = 2*np.pi/lmda * n_int
    k_int = k_int[0,0]
    print("kint : ", k_int)

    e_int = perm
    e_ext = 1.0

    #####################################
    a = 0.1                             # 20 cm diameter sphere
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
    c1 = axs[0].pcolormesh(X[0],X[1], np.abs(Ex)**2)
    cbar1 = fig.colorbar(c1, ax=axs[0])
    c2 = axs[1].pcolormesh(X[0],X[1], np.abs(Ey)**2)
    cbar2 = fig.colorbar(c2, ax=axs[1])
    c3 = axs[2].pcolormesh(X[0],X[1], np.abs(Ez)**2)
    cbar3 =  fig.colorbar(c3, ax=axs[2])

    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")

    Eqs = qs.quasi_stat_field([x,y,z], e_int, e_ext, a)    

    Eqsx = Eqs[:,0].reshape(X[0].shape).T
    Eqsy = Eqs[:,1].reshape(X[0].shape)
    Eqsz = Eqs[:,2].reshape(X[0].shape)
    
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True )
    c1=axs[0].pcolormesh(X[0],X[1], np.abs(Eqsx)**2)
    cbar = fig.colorbar(c1, ax=axs[0])
    c2=axs[1].pcolormesh(X[0],X[1], np.abs(Eqsy)**2)
    cbar2 = fig.colorbar(c2, ax=axs[1])
    c3=axs[2].pcolormesh(X[0],X[1], np.abs(Eqsz)**2)
    cbar = fig.colorbar(c3, ax=axs[2])
    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")   
    plt.show()


       
    """
    matFile = io.loadmat('./res/resultatsMie_1_MHz.mat')
    E_mat = np.flipud( np.rot90(np.array(matFile['E']) ) )


    fig, axs = plt.subplots(1,3, sharex=True, sharey=True )
    axs[0].pcolormesh(X[0],X[1], np.abs(Ex -E_mat[:,:,0])**2)
    axs[1].pcolormesh(X[0],X[1], np.abs(Ey - E_mat[:,:,1])**2)
    axs[2].pcolormesh(X[0],X[1], np.abs(Ez - E_mat[:,:,2])**2)
    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")
    plt.show()

    print(np.nanmean(np.abs(Ex-E_mat[:,:,0])**2))
    """

    ##end main





if __name__ == "__main__":
    main()