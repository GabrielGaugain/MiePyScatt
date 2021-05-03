import numpy as np
import scipy.io as io

a = np.array([ [[1, 0],
               [0, 1]],
               [[2, 0],
               [0, 2]],
               [[3, 0],
               [0, 3]],
               [[3, 0],
               [0, 3]],               
            ])

b = np.array([[4, 1],
              [2, 2],
              [4, 1],
              [1, 2],])


print(a.shape)
print(b.shape)
print(np.shape(b.T))
#np.expand_dims(b, axis=0)

print( a.dot(b))



matFile = io.loadmat('./resultatsMie_0.001_MHz.mat')
E = np.array(matFile['E'])
