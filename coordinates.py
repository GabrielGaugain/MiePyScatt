import numpy as np
import time


def cart2sph(x,y,z):

    #tic = time.time()
   
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.arctan2(y,x)
    theta = np.arctan2( np.sqrt(x**2 + y**2), z  )

    #toc = time.time()
    #print("tooks "+ str(toc-tic) +" s to convert coorinates")
    return r, theta, phi




def Mat_sph2cart(theta,phi):
    """
    Function to compute the matrix to transform the basis vector inspherical coordinates
    to cartesian one or in other word to transform a vector in spherical coordinates to cartesian ones
    """

    assert np.size(theta)==np.size(phi), "Theta and phi must have the same size"

    M = np.zeros((np.size(theta),3,3))


    M[:,0,0]= np.sin(theta) * np.cos(phi)
    M[:,0,1]= np.cos(theta) * np.cos(phi)
    M[:,0,0]= - np.sin(phi)

    M[:,1,0]= np.sin(theta) * np.sin(phi)
    M[:,1,1]= np.cos(theta) * np.sin(phi)
    M[:,1,2]= np.cos(phi)

    M[:,2,0]= np.cos(theta) 
    M[:,2,1]= -np.sin(theta) 

    return M

def Mat_cart2sph(theta,phi):
    """
    Function to compute the matrix to transform the basis vector in cartesian coordinates
    to spherical ones or in other word to transform a vector in cartesian coordinates to spherical ones
    """

    return Mat_sph2cart(theta,phi).T