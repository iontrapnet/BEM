import matlab.engine
import itertools, os
import numpy as np
from scipy.io import loadmat

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'

MSession = {}

def MEngine(*args):
    if not MSession.get('eng', None):
        eng = matlab.engine.start_matlab(*args)
        eng.cd('../MATLAB')
        MSession['eng'] = eng
    return MSession['eng']

def MFieldInit(path, xr=[], yr=[], zr=[], data=''):
    if xr and yr and zr:
        MSession['xr'] = xr
        MSession['yr'] = yr
        MSession['zr'] = zr
        if not data:
            data = MEngine().DataHash(matlab.double(xr+yr+zr))
        matfile = path + '-' + data + '.mat'
        if not os.path.exists(matfile):
            MEngine().FieldInit(path,matlab.double(xr),matlab.double(yr),matlab.double(zr))
        pb = loadmat(matfile)['pb']
        if pb.ndim == 2:
            pb = pb.reshape(pb.shape[0],1,pb.shape[1])
    else:
        pb = MEngine().FieldInit(path)
    MSession['pb'] = pb

def flatten(seq):
    return list(itertools.chain.from_iterable(seq))

def EqualRangeInterpolation(point, xr, yr, zr, *args):
    '''
    # test
    r = range(-10,11)
    x = flatten([[i]*(21*21) for i in r])
    y = flatten([[i]*21 for i in r])*21
    z = r*(21*21)
    print EqualRangeInterpolation([3.7,-8.4,9.9],[-10,10,20],[-10,10,20],[-10,10,20],x,y,z)
    > [3.6999999999999993, -8.4, 9.899999999999999]
    '''
    nx = xr[2] + 1
    ny = yr[2] + 1
    nz = zr[2] + 1
    dx = 0
    dy = 0
    dz = 0
    if nx > 1:
       dx = (nx-1)*(point[0]-xr[0])/(xr[1]-xr[0])
    if ny > 1:
       dy = (ny-1)*(point[1]-yr[0])/(yr[1]-yr[0])
    if nz > 1:
        dz = (nz-1)*(point[2]-zr[0])/(zr[1]-zr[0])
    ix = int(dx);
    iy = int(dy);
    iz = int(dz);
    dx -= ix;
    dy -= iy;
    dz -= iz;
    coe = [(1-dx)*(1-dy)*(1-dz),(1-dx)*(1-dy)*dz,(1-dx)*dy*(1-dz),(1-dx)*dy*dz,dx*(1-dy)*(1-dz),dx*(1-dy)*dz,dx*dy*(1-dz),dx*dy*dz]
    idx = [iz+iy*nz+ix*ny*nz,1+iz+iy*nz+ix*ny*nz,iz+(1+iy)*nz+ix*ny*nz,1+iz+(1+iy)*nz+ix*ny*nz,iz+iy*nz+(1+ix)*ny*nz,1+iz+iy*nz+(1+ix)*ny*nz,iz+(1+iy)*nz+(1+ix)*ny*nz,1+iz+(1+iy)*nz+(1+ix)*ny*nz]
    r = []
    for arg in args:
        r.append(sum([coe[i]*arg[idx[i]] for i in range(8)]))
    if len(r) == 1:
        return r[0]
    return r

def MField(voltages, points):
    pb = MSession['pb']
    if np.ndim(pb) == 0:
        return MEngine().Field(matlab.double(voltages), matlab.double(points))
    else:
        pb = np.dot(voltages,pb)
        xr = MSession['xr']
        yr = MSession['yr']
        zr = MSession['zr']
        return [EqualRangeInterpolation(point, xr, yr, zr, pb[1,:], pb[2,:], pb[3,:]) for point in points]

if __name__ == '__main__':
    points = [[0,0,0],[0,-1,-2],[1,1,1]]
    MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40],'7a5fb13c63c4d748307ca9824e34bdd7')
    print MField([1],points)
    MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40])
    print MField([1],points)
    MFieldInit(PATH_SPHERE)
    print MField([1],points)