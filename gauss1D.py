# Gauss points, weights and basis functions for 1D
# Part of laplace.py

#import statements
import numpy as np 


def getGaussLegendre1D(nIP):
  
  if nIP == 1:
    z = 0
    w = 2
  elif nIP == 2:
    z = ([[-1/np.sqrt(3)],[1/np.sqrt(3)]])
    w = ([1,1]) 
  elif nIP == 3:
    z = ([[-np.sqrt(15)/5],[0],[np.sqrt(15)/5]])
    w = ([5/9,8/9,5/9]) 
  elif nIP == 4:
    z = ([[-0.8611363115940525],[-0.3399810435848563],[0.3399810435848563],[0.8611363115940525]]) 
    w = ([0.3478548451374541,0.6521451548625461,0.6521451548625461,0.3478548451374541])
  else:
    print('unavailable quadrature')

  return z,w

def getBasisFunction1D(element,zgp):

  xi = [item[0] for item in zgp] #zgp[:,0]
  nen = element.degree+1
  if nen == 2:
    N = ([np.array([1-x for x in xi])/2,np.array([1+x for x in xi])/2])
  elif nen == 3:
    N = ([xi*(np.array([x-1 for x in xi]))/2,(np.array([1+x for x in xi]))*(np.array([1-x for x in xi])),xi*(np.array([x+1 for x in xi]))/2])
  elif nen == 4:
    N = ([-1/16*(np.array([x-1 for x in xi]))*(3*(np.array([x-1 for x in xi])))*(3*(np.array([x+1 for x in xi]))),9/16*(np.array([x-1 for x in xi]))*(3*(np.array([x-1 for x in xi])))*(np.array([x+1 for x in xi])),-9/16*(np.array([x-1 for x in xi]))*(3*(np.array([x+1 for x in xi])))*(np.array([x+1 for x in xi])),1/16*(3*(np.array([x-1 for x in xi])))*(3*(np.array([x+1 for x in xi])))*(np.array([x+1 for x in xi]))])
  else: 
    print('Error in ShapeFunc_1D: wrong number of nodes')

  return N
