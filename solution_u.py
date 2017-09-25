# Computing solution = u
# Part of laplace.py

import numpy as np
import gauss1D

def getSolution(K,f,BCDir,BCNeu,X,element):

  degree = element.degree
  ndofT = len(X)
  nen = 2 
  nDim = element.nDim
  nIP = degree+1
  z,w = gauss1D.getGaussLegendre1D(nIP)
  N = gauss1D.getBasisFunction1D(element,z)

  Xe = np.zeros((nen,nDim))
  Te = np.zeros((1,nen))
  Neu = np.zeros((ndofT,1))
  for Neu_ie in range(0,len(BCNeu)):
    if degree == 1:
      Te = ([int(BCNeu[Neu_ie,0]),int(BCNeu[Neu_ie,1])])
      for i in range(0,nen):
        Xe[i,:] = X[Te[i],:]
      h = np.sqrt((Xe[-1,0]-Xe[0,0])**2+(Xe[-1,1]-Xe[0,1])**2)
      neu = np.zeros((1,len(N[0])))
      for ig in range(0,nIP):
        N_ig = np.array(N)[ig,:]
        w_ig = np.array(w)[ig]*h*0.5
        neu = neu + w_ig * N_ig*BCNeu[Neu_ie,2]
      for i in range(0,len(N[0])):
        Neu[Te[i],:] = Neu[Te[i]] + neu[0,i]
  kk = np.zeros(((ndofT+1)*(ndofT+1),1))
  for i in range(0,ndofT):
    for j in range(0,len(BCDir)):
      kk[i,:] = kk[i,:] + K[i,int(BCDir[j,0])]*BCDir[j,1] 
    f[i,:] = f[i,:] - kk[i,:] + Neu[i,:]
  udofs = np.setdiff1d(np.arange(0,ndofT,1),BCDir[:,0])
  Kred = np.zeros((len(udofs),len(udofs)))
  fred = np.zeros((len(udofs),1))
  for i in range(0,len(udofs)):
    for j in range(0,len(udofs)):
      Kred[i,j] = K[udofs[i],udofs[j]]
    fred[i] = f[udofs[i]]
  sol = np.zeros((len(udofs)))
  sol = np.linalg.solve(Kred,fred)
  print('Reduced K')
  print(fred) 
  u = np.zeros((ndofT,1))
  for i in range(0,len(BCDir)):
    u[int(BCDir[i,0]),0] = BCDir[i,1]
  
  for i in range(0,len(udofs)):
    u[udofs[i]] = sol[i]

  return u
