# -*- coding: utf-8 -*-
# INT_GREEN3D_TRI integrates '1/r' and 'grad(1/r)' over a triangle
#
# USE:
# I1, Igrad = int_green3d_TRI(Pglo, V)
#
# INPUTS:
# 'Pglo': point(s) to calculate integral (in in global 3d coordinates)
# 'V': vertices of triangle (in global 3d coordinates)
#
# OUTPUTS:
# 'I1': value of the integral of '1/r'
# 'Igrad': value of the integral of 'grad(1/r)' in global 3d coordinates
# 
# DEMO:
# 
# If you run the following code:
#
# V = [[0,0,0],
#      [3,0,0],
#      [0,3,0]]
# Pglo = [[1,1,0],
#        [-3,3,0],
#        [2,0,2],
#        [-4,1,0],
#        [0,0,1]]
# I1, Igrad = int_green3d_tri(Pglo,V)
#
# You will get
#
# A 5*1 array I1:
# array([[ 7.22168977],
#       [ 1.03049665],
#       [ 1.78031383],
#       [ 0.90645893],
#       [ 2.48777229]])
# and a 5*3 array Igrad:
# array([[-0.24666258, -0.24666258,  0.        ],
#       [ 0.218763  , -0.10034869,  0.        ],
#       [-0.2156998 ,  0.22326257, -0.60312755],
#       [ 0.18338431,  0.00411424,  0.        ],
#       [ 0.66977517,  0.66977517, -0.95824159]])
#
# WARNING:
# Igrad will be diverge if the respective point lies on edges. We use a
# threshold to keep points away from edges. So if a input point is too
# close to the triangle edges, a huge error will be introduced.
#
# NOTE:
# See R. Graglia, "Numerical Integration of the linear shape functions
# times the 3-d Green's function or its gradient on a plane triangle", IEEE
# Transactions on Antennas and Propagation, vol. 41, no. 10, Oct 1993,
# pp. 1448--1455
#
# Adapted from MATLAB Version:
# Date: 03.03.2010
# Copyright(C) 2010-2014: Fabio Freschi (fabio.freschi@polito.it)
#
# This is a Python Version:
# Date: 09.07.2017
# by 郭伟轩 (cnmozzie@gmail.com)

import numpy as np

def versor(v):
    return v / np.linalg.norm(v)

def int_green3d_tri(Pglo, V):
    
    # number of field points
    nfp = int(np.size(Pglo)/3)
    
    # turn the input to array
    Pglo = np.array(Pglo)
    Pglo = Pglo.reshape(nfp, 3).transpose() # use different way to store P
    V = np.array(V)
    
    # local-to-global rotation matrix
    vec1 = V[1]-V[0]
    vec2 = V[2]-V[0]
    vec3 = V[2]-V[1]
    dx = versor(vec1)                      # local x-axis
    dz = versor(np.cross(vec1, vec2))      # local z-axis
    dy = versor(np.cross(dz, dx))          # local y-axis
    R = (np.mat([dx,dy,dz]).T).I  # the inverse of matrix [dx,dy,dz]
    
    # field points in local coordinates
    Ploc = Pglo-V[0].reshape(3,1)        # translation of field points Pglo
    Ploc = np.array(R*np.mat(Ploc))    # rotate the coordinate
    u0 = Ploc[0]                        # notation according to Graglia's paper
    v0 = Ploc[1]
    w0 = Ploc[2]
    
    # vertices in local coordinates
    Vloc = V-V[0]  # translation
    Vloc = np.array((R*np.mat(Vloc).T).T)      # rotate the coordinate
    l3 = Vloc[1][0]                     # notation according to Graglia's paper
    u3 = Vloc[2][0]
    v3 = Vloc[2][1]
    
    # edge lengths
    l1 = np.linalg.norm(vec3)
    l2 = np.linalg.norm(vec2)
    
    # threshold for small numbers
    threshold = 1e-6 * np.min([l1,l2,l3])
    
    # make sure Igrad behaves well near the plane
    w0[abs(w0) < threshold] = 0   
    
    # versors normal to edges Fig. 1(b)
    m = np.array([    # m = s x w, with w = [0 0 1]
        versor([v3,l3-u3,0]),  # m1
        versor([-v3,u3,0]),    # m2
        versor([0,-l3,0])      # m3
    ])     

    # useful quantities for integration
    sminus = np.array([-((l3-u3)*(l3-u0)+v3*v0)/l1, 
                        -(u3*(u3-u0)+v3*(v3-v0))/l2, -u0]) # eq. (3)
    
    splus = np.array([((u3-l3)*(u3-u0)+v3*(v3-v0))/l1, 
                                 (u3*u0+v3*v0)/l2, l3-u0]) # eq. (3)
    
    t0 = np.array([((u3-l3)*v0+v3*(l3-u0))/l1, (v3*u0-u3*v0)/l2, v0]) # eq. (4)
    
    tplus = np.array([np.sqrt((u3-u0)**2+(v3-v0)**2), 
                    np.sqrt(u0**2+v0**2), np.sqrt((l3-u0)**2+v0**2)]) # eq. (5)

    tminus = np.array([tplus[2],tplus[0], tplus[1]])  # eq. (5)
    
    R0 = np.sqrt(t0**2+w0**2)                # line 1 pp. 1450
    
    Rminus = np.sqrt(tminus**2+w0**2)        # line 2 pp. 1450
    
    Rplus = np.sqrt(tplus**2+w0**2)          # line 2 pp. 1450
    
    # initialize f2 and beta matrices
    f2 = np.zeros([3,nfp]);
    beta = np.zeros([3,nfp]);
    
    # field point not in the plane of triangle 
    id1 = np.where(abs(w0)  >= threshold)
    f2[:,id1] = np.log((Rplus+splus)/(Rminus+sminus))[:,id1]    # eq. (11)
    beta[:,id1] = (np.arctan(t0*splus/(R0**2+abs(w0)*Rplus))-
            np.arctan(t0*sminus/(R0**2+abs(w0)*Rminus)))[:,id1] # eq. (14)
    
    # field point in the plane of triangle
    id2 = np.array(np.where(abs(w0)  < threshold))
    id2 = id2.reshape(1,id2.size)
    
    for i in id2[0]:
        for j in range(3):
            # swap the the end points of the integral over 1/r for convenience
            # the absolute value will not be changed
            if abs(splus[j][i])<abs(sminus[j][i]):
                splus[j][i], sminus[j][i] = -sminus[j][i], -splus[j][i]
                tplus[j][i], tminus[j][i] = tminus[j][i], tplus[j][i]
             # field point lying at vertices
            if abs(tminus[j][i])  < 2*threshold:  
                f2[j][i] = np.log(tplus[j][i]/threshold) 
                beta[j][i] = 0
            # field point lying on edges
            elif abs(t0[j][i])  < threshold and sminus[j][i] < 0: 
                f2[j][i] = np.log(abs(2*sminus[j][i]*
                                    (tplus[j][i]+splus[j][i]))/(threshold**2)) 
                beta[j][i] = 0
            # field point not aligned with edges
            else:                                   
                f2[j][i] = np.log((tplus[j][i]+splus[j][i])/
                                        (tminus[j][i]+sminus[j][i])) # eq. (15)
                beta[j][i] = np.arctan(splus[j][i]/t0[j][i])-np.arctan(
                                            sminus[j][i]/t0[j][i])   # eq. (17)
                
    # integral value of '1/r'
    I1 = np.sum(t0*f2-abs(w0)*beta, axis=0).reshape(nfp,1) # eq. (19)

    # integral value of grad(1/r)
    Igradloc = np.array([                                  # eq. (34)
    -m[0][0]*f2[0]-m[1][0]*f2[1]-m[2][0]*f2[2],
    -m[0][1]*f2[0]-m[1][1]*f2[1]-m[2][1]*f2[2],
    -np.sign(w0)*np.sum(beta, axis=0)])     
    Igrad = np.array((np.mat([dx,dy,dz]).T*np.mat(Igradloc)).transpose())
    
    return I1, Igrad
