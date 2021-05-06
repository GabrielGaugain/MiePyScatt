import numpy as np
import matplotlib.pyplot as plt
import mie_field_calculation as mie
import Gabriel_eq as diel
import scipy.constants as ctes
import scipy.io as io
import coordinates as coor
import quasi_static as qs


def internal_field_comparison(X, k_int, k_ext, a):


    (r, theta, phi) = coor.cart2sph(X[0],X[1],X[2])

    E = np.zeros((np.size(r),3), dtype=complex)

    r_int, theta_int, phi_int = r[r<a], theta[r<a] , phi[r<a]    
    r_ext, theta_ext, phi_ext = r[r>=a], theta[r>=a] , phi[r>=a]

    #print("Fields calculations ...")
    E_int = mie.E_internal(r_int,theta_int,phi_int, k_int,k_ext,a)
    #print("Fields computed !")

    for i in range(theta_int.size):
        R= coor.Mat_sph2cart(theta_int[i],phi_int[i])
        E_int[i,:] = R@E_int[i,:]

    E[r<a,:] = E_int
    
    m = (k_int/k_ext)**2

    Eqs = np.ones(E_int[:,0].shape)* 3.0/(m+2)

    rel_err_field = np.abs(np.abs(E_int[:,0])- Eqs)/Eqs
    
    plt.figure()
    plt.plot(np.abs(E_int[:,0]))
    plt.plot(Eqs)
    plt.show()

    return rel_err_field



def main():

    ## frequency(ies) of interest
    #f = 10**6
    freqs = np.logspace(1,8, 8)
    print(freqs)
    lmda = ctes.c / freqs

    n_ext = 1.0                   # outer refrective index : 1 for the air
    k_ext = 2*np.pi/lmda * n_ext 

    ## medium of the sphere
    tissue = "Muscle"
    perm, cond =  diel.get_diel_properties(tissue, freqs=freqs)

    n_int = np.sqrt(perm)   ## without loss => no cmplx part contrib    
    k_int = 2*np.pi/lmda * n_int
    k_int = k_int[0,:]
    print("kint : ", k_int)

    e_int = perm
    e_ext = 1.0

    #####################################
    a = 0.1                             # 20 cm diameter sphere
    L_x, L_y, L_z = 0.4, 0.4, 0.4       # size of the grid on which to compute the fields
    N_x, N_y, N_z = 150,150,150         # number of points in each direction of the grid 
   
    ## 2D, in the plane z=0
    X = np.meshgrid( np.linspace(0,L_x,N_x) - L_x/2, \
                     np.linspace(0,L_y,N_y) - L_y/2, \
                    )
    x,y = X[0].ravel(), X[1].ravel()
    z = np.zeros(x.shape)

    avg_E = np.zeros(freqs.shape)
    max_E = np.zeros(freqs.shape)

    for (i,f) in enumerate(freqs):

        rel_error = internal_field_comparison([x,y,z], k_int[i], k_ext[i], a)
    
        avg_E[i] = np.mean(rel_error)
        max_E[i] = np.max(rel_error)


    plt.figure()
    plt.plot(freqs,avg_E)
    plt.plot(freqs,max_E)
    plt.show()

    """    
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True )
    c1 = axs[0].pcolormesh(X[0],X[1], np.abs(Ex)**2)
    cbar1 = fig.colorbar(c1, ax=axs[0])
    c2 = axs[1].pcolormesh(X[0],X[1], np.abs(Ey)**2)
    cbar2 = fig.colorbar(c2, ax=axs[1])
    c3 = axs[2].pcolormesh(X[0],X[1], np.abs(Ez)**2)
    cbar3 =  fig.colorbar(c3, ax=axs[2])

    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")
    plt.show()
    """
    ##end main





if __name__ == "__main__":
    main()