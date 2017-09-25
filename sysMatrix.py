#Computing system matrix
#Part of laplace.py

from scipy import sparse
import numpy as np

def getSystemMatrix(X,T,element):

  numGauss = element.numGauss
  w_gp = element.w_gp
  N = element.N
  Nxi = element.Nxi
  Neta = element.Neta
  nen = element.nen
  nElems = element.nElems
  ndofTe = element.nen
  ndofT = len(X)
  nNodes = element.nNodes
  nDim = element.nDim

  K = np.zeros((ndofT,ndofT))
  f = np.zeros((ndofT,1))
  Xe_n = np.zeros((nen,nDim))
  Te = np.zeros((1,nen))
  Te_dof = np.zeros((1,ndofTe))
  
  #Assembly (Sparse)
  ind = 0
  n = 2*((ndofTe**2)*nElems)
  ind_i = np.zeros((1,n))
  ind_j = np.zeros((1,n))
  coef_K = np.zeros((1,n))

  for ielem in range(1,nElems+1):
    Te = T[ielem-1,:].astype(int)
    Te_dof = Te
    # 2 dofs per node
    #for i in range(0,len(Te)):
      #Te_dof[0,i] = 2*Te[i]-1
      #Te_dof[0,i+nen] = 2*Te[i]
    for i in range(0,nen):
      Xe_n[i,:] = X[Te[i]-1,:]
    Xe = np.matrix(Xe_n)
    #Element matrix
    Ke,fe = getElementMatrix(Xe,N,Nxi,Neta,ndofTe,numGauss,w_gp,nDim)
    
    #assembly using connectivity
    #for i in range(0,nen):
    #  for j in range(0,nen):
    #    K[Te[i]-1,Te[j]-1] = K[Te[i]-1,Te[j]-1] + Ke[i][j]
    #  f[Te[i]-1,:] = f[Te[i]-1] + fe[i]
  
    #Assembly (Sparse)

    for irow in range(0,ndofTe):
      for icol in range(0,ndofTe):
        ind_i[0,ind] = Te_dof[irow]
        ind_j[0,ind] = Te_dof[icol]
        coef_K[0,ind] = Ke[irow,icol] 
        ind = ind+1
    for i in range(0,nen):
     f[Te[i]-1,:] = f[Te[i]-1] + fe[i]
  
  #Construction of the sparse matrix from the vectors of coefficients, row and column indices
  ind_ii = ind_i[0,0:ind]
  ind_jj = ind_j[0,0:ind]
  coef_Ke = coef_K[0,0:ind]
  K = sparse.coo_matrix((coef_Ke,(ind_ii-1,ind_jj-1)),shape=(ndofT,ndofT)).tocsr()

  return K,f


#---- Generates Elemental matrix -------------------------
def getElementMatrix(Xe,N,Nxi,Neta,ndofTe,numGauss,w_gp,nDim):

  Ke = np.zeros((ndofTe,ndofTe))
  fe = np.zeros((ndofTe,1))

  f_ig = 0
  nen = len(Xe)

  for ig in range(0,numGauss):
    #N_ig = np.array([item[ig] for item in N])
    N_ig = np.array(N)[ig,:]
    Nxi_ig = np.array(Nxi)[ig,:]
    Neta_ig = np.array(Neta)[ig,:]
    Jacob = ([[Nxi_ig],[Neta_ig]])*Xe
    J = np.linalg.det(Jacob)
    dvolu = np.array(w_gp)[ig]*J
    res = np.zeros((nDim,len(Nxi)))
    for i in range(0,len(Nxi)):
      res_n = np.linalg.solve(Jacob,([[Nxi_ig[i]],[Neta_ig[i]]])) 
      res[:,[i]] = res_n
    dNdx_ig = res[0,:]
    dNdy_ig = res[1,:]
    Ke = Ke + (dNdx_ig[:, np.newaxis]*dNdx_ig+dNdy_ig[:, np.newaxis]*dNdy_ig)*dvolu; 
    fe = fe + N_ig[:, np.newaxis]*f_ig*dvolu;
  return Ke,fe

