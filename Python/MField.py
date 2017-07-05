import matlab.engine

# load BemSolver (32 bit dll)
#from ctypes import *
#import os
#os.chdir(os.path.dirname(__file__))
#os.chdir('../Lua')
#print windll.BEM

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'

MSession = {}

def MEngine(*args):
    if not MSession.get('eng', None):
        eng = matlab.engine.start_matlab(*args)
        eng.cd('../MATLAB')
        MSession['eng'] = eng
    return MSession['eng']

def MFieldInit(path, xr=[], yr=[], zr=[]):
    eng = MEngine()
    if xr and yr and zr:
        pb = eng.FieldInit(path,matlab.double(xr),matlab.double(yr),matlab.double(zr))
    else:
        pb = eng.FieldInit(path)
    MSession['pb'] = pb

def EqualRangeInterpolation(point, xr, yr, zr, *args):
    '''
# test
> print EqualRangeInterpolation([0.1,-0.3,0.9],[-1,1,10],[-1,1,10],[-1,1,10],range(1000))
> 713.5
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
        dz = (nz-1)*(point[2]-zr[0])/(zr[1]-zr[0]);
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
    eng = MEngine()
    #pb = MSession['pb']
    return eng.Field(matlab.double(voltages), matlab.double(points))
   	
if __name__ == '__main__':
    points = [[0,0,0],[0,-1,-2],[1,1,1]]
    MFieldInit(PATH_SPHERE)
    print MField([1],points)
    MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40])
    print MField([1],points)