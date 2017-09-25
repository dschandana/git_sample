#! /usr/bin/env python3



# Solves Laplace 2D equation using Finite element method
# Problem description: Square domain
# Specify:
# Discretization size
#         nx = ny = number of elements
#         per side in x and y directions
# Element type
#         'Tri'  for triangular elements
#         'Quad' for rectangular elements
# Degree   1 (linear)
#          2 (quadratic)

# Code by: Sai Chandana Divi (s.c.divi@tue.nl)
# Doctoral student, TU/e
# 01-09-2017

#from libraries import *
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

#Defined functions/objects
from structure import structtype
import mesh_gen
import setRef
import sysMatrix
import solution_u

def main(degree,nx,ny,dom,elemType):
  
  element = structtype()
  element.degree = degree
  element.elemType = elemType
  element.dom = dom

  #mesh_generation
  X,T,BCDir,BCNeu,element = mesh_gen.mesh(element,nx,ny)
  
  #mesh_loading
  #X,T = mesh_load.mesh(degree,n,elemType)
  
  element = setRef.setReferenceElement(element) 
  K,f = sysMatrix.getSystemMatrix(X,T,element)
  

  print('Coordinates')
  print(X)
  print('Connectivity')
  print(T)
  print('Basis functions')
  print(element.N)
  print('Derivatives')
  print(element.Nxi)
  print('second derivatives')
  print(element.Neta)
  print('K matrix')
  print(K.toarray()) 
  print('RHS')
  print(f)
  print('Boundary conditions')
  print(BCDir)
  print('Neumann conditions')
  print(BCNeu)
 
  u = solution_u.getSolution(K,f,BCDir,BCNeu,X,element)
  print('Solution')
  print(u)
 
  fig = plt.figure()
  ax = fig.gca(projection='3d')
  
  # Make data.
  x = np.array([row[0] for row in X])
  y = np.array([row[1] for row in X])
  z = np.array([row[0] for row in u])

  print(z)
  # Plot the surface.
  surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,linewidth=0, antialiased=False)


  # Customize the z axis.
  ax.set_zlim(0,2)
  ax.zaxis.set_major_locator(LinearLocator(10))
  ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

  # Add a color bar which maps values to colors.
  fig.colorbar(surf, shrink=0.5, aspect=5)
  
  plt.show()

  return

main(degree=1,nx=2,ny=2,dom=([0,1,0,1]),elemType='Quad')
