import matlab.engine

# load BemSolver (32 bit dll)
#from ctypes import *
#import os
#os.chdir(os.path.dirname(__file__))
#os.chdir('../Lua')
#print windll.BEM

PATH_SPHERE = '..\\Model\\sphere\\922019645250603132'
PATH_4ROD = '..\\Model\\4rod\\167634622912717531'

def MStart():
	global eng
	eng = matlab.engine.start_matlab()
	eng.cd('../MATLAB')

def MFieldInit(path, xr=[], yr=[], zr=[]):
	global eng
	if xr and yr and zr:
		eng.FieldInit(path,matlab.double(xr),matlab.double(yr),matlab.double(zr),nargout=0)
	else:
		eng.FieldInit(path,nargout=0)

def MField(voltages, points, method='Calc'):
	global eng
	if method == 'Calc':
		func = eng.FieldCalc
	else:
		func = eng.FieldInterp
	return func(matlab.double(voltages), matlab.double(points))

if __name__ == '__main__':
	MStart()
	points = [[0,0,0],[0,-1,-2],[1,1,1]]
	MFieldInit(PATH_SPHERE)
	print MField([1],points)
	MFieldInit(PATH_SPHERE,[-2,2,40],[-2,2,40],[-2,2,40])
	print MField([1],points,method='Interp')