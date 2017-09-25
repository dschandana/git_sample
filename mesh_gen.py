# Generate mesh
# Part of laplace.py

import numpy as np
from math import ceil

def mesh(element,nx,ny):
  elemType = element.elemType
  degree = element.degree
  dom = element.dom
  element.nDim = 2
  x1 = dom[0]
  x2 = dom[1]
  y1 = dom[2]
  y2 = dom[3]

  npx = nx*degree+1
  npy = ny*degree+1
  npt = npx*npy
  hx = np.linspace(x1,x2,npx)
  hy = np.linspace(y1,y2,npy)

  xn,yn = np.meshgrid(hx,hy)
  x = np.reshape(xn,npt)
  y = np.reshape(yn,npt)
  
  X = np.zeros((len(x),element.nDim))
  for i in range(0,len(x)):
    X[i,:] = ([x[i],y[i]])

  if elemType == 'Quad':
    nElems = nx*ny
    if nx == ny:
      n = nx
    if degree == 1:
      nen = 4
      nNodes = (n+1)*(n+1)
      T = np.zeros((nElems,int(nen)))  #[[0] * nen for iT in range(0,nElems)]
      for n_elem in range(1,nElems+1):
        line = ceil((n_elem)/n)
        print(nen)
        T[n_elem-1,:] = ([n_elem + line - 1,n_elem+line,n_elem+line+n+1,n_elem+line+n])
     
    
    if degree == 2:
      T = np.zeros((nx*ny,9))
      for i in range(0,ny):
        for j in range(0,nx):
          ielem = (i-1)*nx + j
          inode = (i-1)*2*npx + 2*(j-1) + 1
          nodes_aux = ([inode+np.linspace(0,2,3), inode+npx+np.linspace(0,2,3), inode+2*npx+np.linspace(0,2,3)])
          nodes_auxx= [l for k in nodes_aux for l in k]
          nodes_auxx = np.array(nodes_auxx)
          T[ielem-1,:] = nodes_auxx[[0, 2, 8, 6, 1, 5, 7, 3, 4]]
 

      print(T)
      # Creating Boundary conditions
      # Problem: left inlet and right outlet with -1 and 1. Dirichlet at (0,0)
     # BCDir = np.array([[1,0]])

     # BCNeu= np.zeros((4*n,3))
     # BCNeu_left_n = 1
     # BCNeu_right_n = 1+n
     # BCNeu[0,0] = BCNeu_left_n
     # BCNeu[0,1] = BCNeu_left_n+(n+1)
     # BCNeu[n,0] = BCNeu_right_n + (n+1)
     # BCNeu[n,1] = BCNeu_right_n
     # BCNeu[2*n,0] = BCNeu_left_n
     # BCNeu[2*n,1] = BCNeu_left_n + 1
     # BCNeu[3*n,0] = (n+1)*(n+1)
     # BCNeu[3*n,1] = (n+1)*(n+1)-1
     # BCNeu[0,2] = -1
     # BCNeu[n,2] = 1
     # BCNeu[2*n,2] = 0
     # BCNeu[3*n,2] = 0

     # for i in range(1,n):
     #   BCNeu_left_n = BCNeu_left_n + (n+1)
     #   BCNeu_right_n = BCNeu_right_n + (n+1)
     #   BCNeu[i,0] = BCNeu_left_n
     #   BCNeu[i,1] = BCNeu_left_n+(n+1)
     #   BCNeu[i,2] = -1
     #   BCNeu[i+n,0] = BCNeu_right_n + (n+1)
     #   BCNeu[i+n,1] = BCNeu_right_n
     #   BCNeu[i+n,2] = 1
     #   BCNeu[i+2*n,0] = BCNeu[i+2*n-1,0] + 1
     #   BCNeu[i+2*n,1] = BCNeu[i+2*n-1,1] + 1
     #   BCNeu[i+2*n,2] = 0
     #   BCNeu[i+3*n,0] = BCNeu[i+3*n-1,0] -1 
     #   BCNeu[i+3*n,1] = BCNeu[i+3*n-1,1] - 1
     #   BCNeu[i+3*n,2] = 0


  # Problem:
  x = np.array([row[0] for row in X])
  y = np.array([row[1] for row in X])
  nodesDir1 = np.where(x == 0)[0] 
  nodesDir2 = np.where(x == 1)[0]
  nodesDir = ([nodesDir1,nodesDir2])
  BCDir_dof = np.unique(nodesDir)
  print(y)
  uD = np.zeros((len(X),1))
  for i in range(0,len(nodesDir1)):
    uD[nodesDir1[i]] = y[nodesDir1[i]]
      
  for i in range(0,len(nodesDir2)):
    uD[nodesDir2[i],:] = 1+y[nodesDir2[i]]
     
    BCDir = np.zeros((2*(n+1),2))
    for i in range(0,len(BCDir_dof)):
      BCDir[i,0] = BCDir_dof[i]
      BCDir[i,1] = uD[BCDir_dof[i]]

  nodesNeu1 = np.where(y == 0)[0]
  nodesNeu2 = np.where(y == 1)[0]
  nodesNeu = ([nodesNeu1,nodesNeu2])
  BCNeu_dof = np.unique(nodesNeu)
  face1 = BCNeu_dof[0:int((len(BCNeu_dof))/2)]
  face2 = BCNeu_dof[int((len(BCNeu_dof))/2):]
  print(BCNeu_dof[-1])
  print(face2)
  Neu_conn = np.zeros((2*n,2))
  Neu_conn[0:n,0] = face1[0:len(face1)-1]
  Neu_conn[0:n,1] = face1[1:len(face1)]
  Neu_conn[n:2*n,1] = face2[0:len(face2)-1]
  Neu_conn[n:2*n,0] = face2[1:len(face2)]
  print(Neu_conn) 

  du_n = np.zeros((len(Neu_conn),1))
  du_n[0:int((len(Neu_conn))/2),:] = -1
  du_n[int((len(Neu_conn))/2):] = 1
  BCNeu= np.zeros((2*n,3))
  BCNeu[:,0] = np.array([row[0] for row in Neu_conn])
  BCNeu[:,1] = np.array([row[1] for row in Neu_conn])
  BCNeu[:,2] = du_n[:,0]
  print(BCNeu)

  element.nen = nen
  element.nElems = nElems
  element.nNodes = nNodes
  return X,T,BCDir,BCNeu,element

