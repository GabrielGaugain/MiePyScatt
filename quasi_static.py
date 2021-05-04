import numpy as np
import scipy.constants as ctes
import matplotlib.pyplot as plt
import coordinates as coor
import Gabriel_eq as diel

## See Jackson p180-184

def quasi_stat_field(X, e_int, e_ext, a, E_0=1):


    m = e_int/e_ext


    (r, theta, phi) = coor.cart2sph(X[0],X[1],X[2])

    E = np.zeros((np.size(r),3))  #, dtype=complex)
    
    r_int, theta_int, phi_int = r[r<a], theta[r<a] , phi[r<a]    
    r_ext, theta_ext, phi_ext = r[r>=a], theta[r>=a] , phi[r>=a]

    E_int = np.zeros((np.size(r_int),3))
    E_sca = np.zeros((np.size(r_ext),3))

    E_int[:,0] =   3.0/(2 + m) * E_0          #e_x
    #E_int[:,0] =   3.0/(2 + m) * E_0 * np.cos(theta_int)        #e_r
    #E_int[:,0] = - 3.0/(2 + m) * E_0 * np.sin(theta_int)        #e_theta
    #E_int[:,0] =            #e_phi

    E_sca[:,0] =  (m -1)/(m +2) * E_0   *a**3 /r_ext**3           #e_r
    E_sca[:,0] =  (m -1)/(m +2) * E_0   *a**3 /r_ext**3*np.sin(theta_ext)            #e_theta
    #E_sca[:,0] =            #e_phi

    ## transform back the fields into cartesian coordinates
    """
    for i in range(theta_int.size):
        R= coor.Mat_sph2cart(theta_int[i],phi_int[i])
        E_int[i,:] = R@E_int[i,:]
    """

    for i in range(theta_ext.size):
        R= coor.Mat_sph2cart(theta_ext[i],phi_ext[i])
        E_sca[i,:] = R@E_sca[i,:]

    E[r<a,:] = E_int
    E[r>=a,:] = E_sca 
    E[r>=a,0] =  E[r>=a,0] +  E_0


    return E

def main():

    ## frequency(ies) of interest
    f = 10**3
    e_ext = 1.0                   # outer relative permittivity : 1 for the air

    ## medium of the sphere
    tissue = "Muscle"
    perm, cond = diel.get_diel_properties(tissue, freqs=f)
    e_int = perm   ## without loss => no cmplx part contrib
    
    print("e_int : ", e_int)

    # GEOMETRY      #####################
    #####################################
    a = 0.1                             # 20 cm diameter sphere
    L_x, L_y, L_z = 0.4, 0.4, 0.4       # size of the grid on which to compute the fields
    N_x, N_y, N_z = 150,150,150         # number of points in each direction of the grid 

    X = np.meshgrid( np.linspace(0,L_x,N_x) - L_x/2, \
                     np.linspace(0,L_y,N_y) - L_y/2, \
                    )
    x,y = X[0].ravel(), X[1].ravel()
    z = np.zeros(x.shape)

    E = quasi_stat_field([x,y,z], e_int, e_ext, a)    

    Ex = E[:,0].reshape(X[0].shape)
    Ey = E[:,1].reshape(X[0].shape)
    Ez = E[:,2].reshape(X[0].shape)
    
    fig, axs = plt.subplots(1,3, sharex=True, sharey=True )
    c1=axs[0].pcolormesh(X[0],X[1], np.abs(Ex)**2)
    cbar = fig.colorbar(c1, ax=axs[0])
    c2=axs[1].pcolormesh(X[0],X[1], np.abs(Ey)**2)
    cbar2 = fig.colorbar(c2, ax=axs[1])
    c3=axs[2].pcolormesh(X[0],X[1], np.abs(Ez)**2)
    cbar3 = fig.colorbar(c3, ax=axs[2])
    axs[0].set_xlabel("x (cm)")
    axs[0].set_ylabel("y (cm)")   
    


    plt.show()


    
if __name__ == "__main__":

    main()