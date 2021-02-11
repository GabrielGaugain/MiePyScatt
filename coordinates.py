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




