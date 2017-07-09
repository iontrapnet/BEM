import os
import numpy as np
from scipy.io import loadmat

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'

MSession = {}

def MEngine(*args):
    if not MSession.get('eng', None):
        import matlab.engine
        eng = matlab.engine.start_matlab(*args)
        eng.cd(os.path.join(os.path.dirname(__file__), '..', 'MATLAB'))
        MSession['eng'] = eng
    return MSession['eng']

def MFieldInit(path, xr, yr, zr, data=''):
    MSession['xr'] = xr
    MSession['yr'] = yr
    MSession['zr'] = zr
    if not data:
        import matlab.engine
        data = MEngine().DataHash(matlab.double(xr+yr+zr))
    matfile = os.path.join(os.path.dirname(__file__), path + '-' + data + '.mat')
    if not os.path.exists(matfile):
        import matlab.engine
        MEngine().FieldInit(path,matlab.double(xr),matlab.double(yr),matlab.double(zr), nargout = 0)
    pb = loadmat(matfile)['pb']
    if pb.ndim == 2:
        pb = pb.reshape(pb.shape[0],pb.shape[1],1)
    MSession['pb'] = pb

def MField(voltages, points, ratio=1):
    xr = MSession['xr']
    yr = MSession['yr']
    zr = MSession['zr']
    pb = MSession['pb']
    nx = xr[2] + 1
    ny = yr[2] + 1
    nz = zr[2] + 1
    r = []
    for point in points:
        dx = (nx-1)*(point[0]-xr[0])/(xr[1]-xr[0])
        dy = (ny-1)*(point[1]-yr[0])/(yr[1]-yr[0])
        dz = (nz-1)*(point[2]-zr[0])/(zr[1]-zr[0])
        ix = int(dx)
        iy = int(dy)
        iz = int(dz)
        dx -= ix
        dy -= iy
        dz -= iz
        coe = [(1-dx)*(1-dy)*(1-dz),(1-dx)*(1-dy)*dz,(1-dx)*dy*(1-dz),(1-dx)*dy*dz,dx*(1-dy)*(1-dz),dx*(1-dy)*dz,dx*dy*(1-dz),dx*dy*dz]
        idx = [iz+iy*nz+ix*ny*nz,1+iz+iy*nz+ix*ny*nz,iz+(1+iy)*nz+ix*ny*nz,1+iz+(1+iy)*nz+ix*ny*nz,iz+iy*nz+(1+ix)*ny*nz,1+iz+iy*nz+(1+ix)*ny*nz,iz+(1+iy)*nz+(1+ix)*ny*nz,1+iz+(1+iy)*nz+(1+ix)*ny*nz]
        r.append(ratio*np.dot(np.dot(pb[1:, idx, :], voltages), coe))
    return r

import time

if __name__ == '__main__':
    #points = [[0,0,0],[0,-1,-2],[1,1,1]]
    #MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40])
    #print MField([1],points)
    #MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40],'7a5fb13c63c4d748307ca9824e34bdd7')
    #print MField([1],points)
    MFieldInit(PATH_4ROD,[-0.005,0.005,100], [-0.005,0.005,100], [2.095,2.105,100], '5d49bb9d6704e98ea598bb82a952b646')
    print MField([1,0,0,0,0,1],[[0,0,2.0999],[0,0,2.0101]])
    #t = time.time()
    #for i in xrange(100000):
    #    MField([1,0,0,0,0,1],[[0,0,2.0999],[0,0,2.0101]])
    #print time.time() - t