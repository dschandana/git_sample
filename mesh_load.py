import numpy as np

def mesh(degree,nelem,elemType):
  
  if elemType == 'Tri':
    if degree == 1:
      X = np.loadtxt("Tri_Lin_Mesh1_Nodes.dat")
      T = np.loadtxt("Tri_Lin_Mesh1_Elements.dat") 
  elif elemType == 'Quad':
      if degree == 1:
        X = np.loadtxt("Quad_Lin_Mesh1_Nodes.dat")
        T = np.loadtxt("Quad_Lin_Mesh1_Elements.dat")
  else:
    print('Error in loading mesh')
  return X,T
