# Setting reference element
# Part of laplace.py

import numpy as np

def setReferenceElement(element):

  #Gauss points
  z_gp,w_gp,numGauss = getQuadrature(element)
  N,Nxi,Neta = getBasisfunctions(element,z_gp)
  
  element.z_gp = z_gp
  element.w_gp = w_gp
  element.numGauss = numGauss
  element.N = N
  element.Nxi = Nxi
  element.Neta = Neta

  return element

def getQuadrature(element):

  #2DQudrature
  elemType = element.elemType
  degree = element.degree
  if elemType == 'Tri':
  #Triangular
    numGauss_array = ([3,4,6])
    numGauss = numGauss_array[degree-1]
    if numGauss == 1:
    #Degree = 1
      z_gp = ([1/3,1/3])
      w_gp = 1
    elif numGauss == 3: 
    #Degree = 2
      z_gp = ([[2/3,1/6],[1/6,2/3],[1/6,1/6]])
      w_gp = ([1/3,1/3,1/3])
    elif numGauss == 4:
    #Degree = 3
      z_gp = ([[1/3,1/3],[0.6,0.2],[0.2,0.6],[0.2,0.2]])
      w_gp = ([-27.0/48.0,25.0/48.0,25.0/48.0,25.0/48.0])
    elif numGauss == 6:
      a = 0.659027622374092
      b = 0.231933368553031
      c = 0.109039009072877
      z_gp = ([[a,b],[b,a],[a,c],[c,a],[b,c],[c,b]])
      w_gp = ([1/6,1/6,1/6,1/6,1/6,1/6]) 
    elif numGauss == 12:
      a = 0.873821971016996
      b = 0.063089014491502
      c = 0.501426509658179
      d = 0.249286745170910
      e = 0.636502499121399
      f = 0.310352451033785
      g = 0.053145049844816
      z_gp = ([[a,b],[b,b],[b,a],[c,d],[d,d],[d,c],[e,f],[f,e],[e,g],[g,e],[f,g],[g,f]])
      w_gp =([0.050844906370207,0.050844906370207,0.050844906370207,0.116786275726379,0.116786275726379,0.116786275726379,0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374,0.082851075618374])
    else:
      print('Error in setting reference element')
  elif elemType == 'Quad':
  #Quadratic
    numGauss = degree*4 +(degree-1)**2
    if numGauss == 4:
    #Degree  = 1
      a = 1/np.sqrt(3)
      z_gp = ([[-a,-a],[a,-a],[a,a],[-a,a]])
      w_gp = ([1,1,1,1])
    elif numGauss == 9:
    #Degree = 2
      a = np.sqrt(3/5)
      z_gp = ([[-a,-a],[0,-a],[a,-a],[-a,0],[0,0],[a,0],[-a,a],[0,a],[a,a]])
      b = 5/9
      c = 8/9
      w_gp= [b*b,c*b,b*b,b*c,c*c,b*c,b*b,c*b,b*b];
    elif numGauss == 16:
      a = np.sqrt(525+70*np.sqrt(30))/35
      b = np.sqrt(525-70*np.sqrt(30))/35
      z_gp = ([[-b,-b],[-b,-a],[-b,a],[-b,b],[-a,-b],[-a,-a],[-a,a],[-a,b],[a,-b],[a,-a],[a,a],[a,b],[b,-b],[b,-a],[b,a],[b,b]])
      c = np.sqrt(30)*(3*np.sqrt(30)-5)/180
      d = np.sqrt(30)*(3*np.sqrt(30)+5)/180
      w_gp = ([d*d,d*c,d*c,d*d,c*d,c*c,c*c,c*d,c*d,c*c,c*c,c*d,d*d,d*c,d*c,d*d])
  else:
    print('Error in setting reference element')
  return z_gp,w_gp,numGauss

def getBasisfunctions(element,z_gp):
  #Basisfunction
  elemType = element.elemType
  degree = element.degree
  if elemType == 'Quad':
    if degree == 1:
      xi = [item[0] for item in z_gp] #z_gp[:,0]
      eta = [item[1] for item in z_gp] #z_gp[:,1]
      #nodesCoord = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0]
      onemxi = np.array([1-x for x in xi])
      onemeta = np.array([1-x for x in eta])
      onepxi = np.array([1+x for x in xi])
      onepeta = np.array([1+x for x in eta]) 
      N = ([(onemxi*onemeta)/4,(onepxi*onemeta)/4,(onepxi*onepeta)/4,(onemxi*onepeta)/4])
      Nxi = ([np.array([-x/4 for x in onemeta]),np.array([x/4 for x in onemeta]),np.array([x/4 for x in onepeta]),np.array([-x/4 for x in onepeta])])
      Nxi = np.transpose(Nxi)
      Neta = ([np.array([-x/4 for x in onemxi]),np.array([-x/4 for x in onepxi]),np.array([x/4 for x in onepxi]),np.array([x/4 for x in onemxi])])
      Neta = np.transpose(Neta)
  return N,Nxi,Neta
