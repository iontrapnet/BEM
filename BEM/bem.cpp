//todo get surface charges per electrode

//this is already done:

//to include fortran libraries read especially "Intel Fortran/C Mixed-Language Programs Overview" in Intel Visual Fortran User Guide Vol I
//but it also works without including libifcore.lib
//you just have to set in the ctest Project property tab at linker/input/additional dependencies to ..\Lib2\Debug\Lib2.lib or  ..\Lib2\Release\Lib2.lib
//and the linker/common/additional library directory to D:\Programme\Intel\Fortran\compiler80\IA32\LIB
//rightclick on Projectmappe to set dependencies and define the file to be debugged and started

//if linker does not find ifconsol.lib then change setting in properties of whole project/linker/allgemein/zusaetzliche bibliotheksverzeichnisse to the proper directory of the intel fortran compiler

// testlapack.cpp : Definiert den Einstiegspunkt für die Konsolenanwendung.
//in the fortran project:
//Remove "/libdir:noauto". In the project property pages, this is Fortran..Libraries..Disable default library search rules - change to "No".
//This setting, which is "Yes" by default in static library projects, tells the compiler to not include in the object file directives to pull in the run-time libraries. This is correct when the library is used by another Fortran project, but not when the main program is not Fortran.
//In order to avoid library collisions also Librarian..General..Ignore all default libraries to yes
//Take care that the project type in Fortran..Libraries is set to multithreaded debug/ multithreaded release respectively (if the current c project has the same settings)
//for old style function definition to work under intel c++ 9.1 use c/language/enable c99 option (Yes (/Qc99))
//for debug mode to work set /MTd or /MT option and ignore specific library libcmtd.lib
//time vc 50000 seconds
//time intel floating point precision optimized: 75587.2

//todo check that if spaceunit is changed all fields are rescaled even when there exists a cache

#ifdef WIN32
#include <windows.h>
#include <aclapi.h>
//#define CHECKLICENSE
#else
//Licence only works in Windows (fz)
#undef CHECKLICENSE
#endif

//#pragma comment(lib, "libCore.lib");
//#pragma comment(lib, "libGeom.lib");

#undef StrDup;

#define _WIN32_WINNT 0x0501
#include <stdio.h>
#include <iostream>


using namespace std;

#include "bem.h"
#include "md5.hpp"
#include <strstream>
#include <iostream>
#include <stdio.h>
#include <complex>
using namespace std;
#include <vector>
#include <list>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>

#ifndef _TIME_			/* if not on a Sun4 */
#include <time.h>
#endif

/* Types of Problems. */
#define FIELD 0                 /* fastlap should only do a field computation.
This is equivalent to computing a RHS.*/
#define GREEN 1                 /* fastlap should solve a linear system with
both right- and left-hand side matrices
as in a Green formulation.*/
#define INDIRECT 2              /* fastlap should solve a linear system with
a left-hand side matrix and a right-hand
side vector as in a single- or double-
layer formulation.*/

/* Shapes of singularities. */
#define POINT 1                 /* Indicates snglrty geometry is a point. */
#define TRIANGLE 3              /* Indicates snglrty geometry is a triangle. */
#define QUADRILAT 4             /* Indicates snglrty geometry is a 
quadrilateral. */

/* Types of singularities. */
#define POINT_SOURCE 1          /* Indicates snglrty is a point source. */
#define CONSTANT_SOURCE 11      /* Indicates snglrty is a constant source. */
#define CONSTANT_DIPOLE 12      /* Indicates snglrty is a constant dipole. */
#define LINEAR_SOURCE 21        /* Indicates snglrty is a linear source. */
#define LINEAR_DIPOLE 22        /* Indicates snglrty is a linear dipole. */

#define XI 0
#define YI 1
#define ZI 2
#define DIRICHLET 0
#define NEUMANN 1
#define DIMEN 3
#define VERTS 4
#define ONE3 0.3333333333333
#define Dot_Product(V1,V2) V1[XI]*V2[XI]+V1[YI]*V2[YI]+V1[ZI]*V2[ZI]
#define DotP_Product(V1,R,S,T) (V1[XI])*(R)+(V1[YI])*(S)+(V1[ZI])*(T)

extern "C" void Cross_Product(double vector1[],double vector2[],double result_vector[]);
extern "C" double normalize(double vector[3]);
extern "C" int fastlap(int *plhsSize,int *prhsSize,int *pnumSing,double *px,int *pshape,int *pdtype,
			int *plhsType,int *prhsType,int *plhsIndex,int *prhsIndex,
			double *plhsVect,double *prhsVect,double *pxf,double *pxnrm,int *pnumLev,int *pnumMom,
			int *pmaxItr,double *ptol,int *pjob);

typedef long int integer;
typedef char *address;
typedef short int shortint;
typedef float real_;
typedef double doublereal;
typedef struct { real_ r, i; } complex_;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
typedef short int shortlogical;
typedef char logical1;
typedef char integer1;
extern "C" int __cdecl DGESV(integer *n, integer *nrhs, doublereal *a, integer 
							 *lda, integer *ipiv, doublereal *b, integer *ldb, integer *info);

extern "C" void freeunusedmem();
//void freeunusedmem() {};

double *dfdnbackup;
/*
 inline void * __cdecl operator new(unsigned int size)
{
  void *ptr = (void *)kcalloc2(size);
//	void *ptr = (void *)malloc(size);
  if(!ptr){
	  fprintf(stderr,"new allocation error. Size: %d\n",size);
	  exit(0);
  }
  return(ptr);
};
inline void __cdecl operator delete(void *p)
{
//	free(p);
  reusemem(p);
  //freeunusedmem();
};
 inline void * __cdecl operator new[](unsigned int size)
{
  void *ptr = (void *)kcalloc2(size);
//	void *ptr = (void *)malloc(size);
  if(!ptr){
	  fprintf(stderr,"new allocation error. Size: %d\n",size);
	  exit(0);
  }
  return(ptr);
};
 inline void __cdecl operator delete[](void *p)
{
//	free(p);
  reusemem(p);
  //freeunusedmem();
}; 
*/

  
double testfn(double x,double y,double z){
	double *feld1=new double[100000000];
	double *feld2=new double[100];
	double *feld3=new double[100];
	delete[] feld2;
	delete[] feld3;
	delete[] feld1;
	return sqrt(x+y);
}

void q(){
#ifdef WIN32
	TerminateProcess(GetCurrentProcess(),0);
#endif
}

double asin_safe(double val){
	if (val>1.) return M_PI/2.;
	else if(val<-1.) return -M_PI/2.;
	else return asin(val);
}

double acos_safe(double val){
	if (val>1.) return 0.;
	else if(val<-1.) return M_PI;
	else return acos(val);
}

calcD3triangle::calcD3triangle():numquad(3),EPS(3.0e-11),n(numquad*numquad){
	xi=new double[numquad*numquad];
	eta=new double[numquad*numquad];
	w=new double[numquad*numquad];
	gauslegsimplextriangle(xi,eta,w,numquad);

};
calcD3triangle::~calcD3triangle(){
	delete[] xi;
	delete[] eta;
	delete[] w;
};



double calcD3triangle::GetIntL1divR3(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double rp){
double a,b,c,asq,bsq,csq,psi,ax,ay,az,bx,by,bz,cx,cy,cz,cospsi,sinpsi;
	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	asq=sqr(ax)+sqr(ay)+sqr(az);
	a=sqrt(asq);
	bsq=sqr(bx)+sqr(by)+sqr(bz);
	b=sqrt(bsq);
	csq=sqr(cx)+sqr(cy)+sqr(cz);
	c=sqrt(csq);

	//a=abs(x3-x1);b=abs(x2-x1);c=abs(x3-x2);

	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
	if(1.-abs(cospsi)<1e-13) return 0;
	if(abs(a)<1e-13) return 0;
	if(abs(b)<1e-13) return 0;
	if(abs(c)<1e-13) return 0;
	psi=acos_safe(cospsi);
	sinpsi=sqrt(1.-sqr(cospsi));

	double s;//psi entspricht phi1p2
	s=asin_safe(rp*(b*cospsi-a)/sqrt(sqr(a*b*sinpsi)+sqr(c*rp))) -
	  asin_safe(rp*(b-a*cospsi)/sqrt(sqr(a*b*sinpsi)+sqr(c*rp)));

	return (s+psi)/(rp);

};

double calcD3triangle::GetIntL1(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3){
	double a,b,c,psi,ax,ay,az,bx,by,bz,cx,cy,cz,cospsi,sinpsi;
	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	a=sqrt(sqr(ax)+sqr(ay)+sqr(az));
	b=sqrt(sqr(bx)+sqr(by)+sqr(bz));
	c=sqrt(sqr(cx)+sqr(cy)+sqr(cz));

	//a=abs(x3-x1);b=abs(x2-x1);c=abs(x3-x2);
	if(abs(a)<1e-13) return 0;
	if(abs(b)<1e-13) return 0;
	if(abs(c)<1e-13) return 0;
	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
    if(1.-abs(cospsi)<1e-13) return 0;
	
	sinpsi=sqrt(1.-sqr(cospsi));


	return  2.*(a*b*sinpsi/4./c)*(log((c-b*cospsi+a)/(c+b*cospsi-a)/(c-b+a*cospsi)*(c+b-a*cospsi)));

};

double calcD3triangle::GetIntL1OutOfPlane(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double rp){
double a,b,c,asq,bsq,csq,psi,ax,ay,az,bx,by,bz,cx,cy,cz,cospsi,sinpsi;
	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	asq=sqr(ax)+sqr(ay)+sqr(az);
	a=sqrt(asq);
	bsq=sqr(bx)+sqr(by)+sqr(bz);
	b=sqrt(bsq);
	csq=sqr(cx)+sqr(cy)+sqr(cz);
	c=sqrt(csq);

	//a=abs(x3-x1);b=abs(x2-x1);c=abs(x3-x2);
	if(abs(a)<1e-13) return 0;
	if(abs(b)<1e-13) return 0;
	if(abs(c)<1e-13) return 0;
	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
    if(1.-abs(cospsi)<1e-13) return 0;

	psi=acos_safe(cospsi);
	sinpsi=sqrt(1.-sqr(cospsi));

	double f4,l2,l3,d,rpsq;//psi entspricht phi1p2
	rpsq=sqr(rp);

	d=acos_safe((rp*(b*cospsi-a))/sqrt(sqr(a*b*sinpsi)+sqr(c*rp)))-
	  acos_safe((rp*(b-a*cospsi))/sqrt(sqr(a*b*sinpsi)+sqr(c*rp)));
	f4=(a*b*sinpsi)/c;
	l2=log((c/a*sqrt(asq+rpsq)-b*cospsi+a)/(c/a*sqrt(asq+rpsq)+b*cospsi-a ));
	l3=log((c/b*sqrt(bsq+rpsq)-b+a*cospsi)/(c/b*sqrt(bsq+rpsq)+b-a*cospsi));

	return (f4/2.*(l2-l3)+rp*(d-psi));

}

//triangleint(&testfn,0,0,0, 1,0,0, 0,1,0);

double calcD3triangle::triangleint(trianglefn *fn,double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3){
	double hs=sqrt(sqr(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) + sqr(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + sqr(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3));
	double sum=0;
	int i;
	for(i=0;i<n;i++){
		double zeta,eta2,xi2;
			
 		eta2=eta[i];
		xi2=xi[i];
		zeta=1-xi2-eta2;
		sum+=fn(xp,yp,zp,x1*zeta+x2*xi2+x3*eta2,y1*zeta+y2*xi2+y3*eta2,z1*zeta+z2*xi2+z3*eta2)*w[i];
	}
	return sum*hs;
	//return sum*hs/2.;//use this if you use trianglequad
};
/*
void calcD3triangle::triangleint_pot_exyz(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm,double &pot,double &ex,double &ey,double &ez){
	pot=0;
	ex=0;
	ey=0;
	ez=0;
	double hs=sqrt(sqr(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) + sqr(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + sqr(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3));
	double sum=0;
	int i;
	
	// b x a
	double nx,ny,nz,nn;
	double ax,ay,az,bx,by,bz;
	ax=(x2-x1);
	ay=(y2-y1);
	az=(z2-z1);
	bx=(x3-x1);
	by=(y3-y1);
	bz=(z3-z1);
	nx=az*by - ay*bz;
	ny=ax*bz - az*bx;
	nz=ay*bx - ax*by;
	nn=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
	if(!inversenorm) nn=-nn;

	nx/=nn;ny/=nn;nz/=nn;

	for(i=0;i<n;i++){
		double zeta,eta2,xi2;
			
 		eta2=eta[i];
		xi2=xi[i];
		zeta=1-xi2-eta2;

		double xq=x1*zeta+x2*xi2+x3*eta2-x0;
		double yq=y1*zeta+y2*xi2+y3*eta2-y0;
		double zq=z1*zeta+z2*xi2+z3*eta2-z0;

		double sqr_xq=sqr(xq);
		double sqr_yq=sqr(yq);
		double sqr_zq=sqr(zq);
		
		double oneover_sqrt_all=1./sqrt(sqr_xq+sqr_yq+sqr_zq);
		double oneover_sqrt_all_3=oneover_sqrt_all*oneover_sqrt_all*oneover_sqrt_all*w[i];
		oneover_sqrt_all*=w[i];

		pot+=oneover_sqrt_all;//speedup idea put pot in array of sizeof(w) and ausklammere w[i]

		ex+=(xq*nx)*oneover_sqrt_all_3;
		ey+=(yq*ny)*oneover_sqrt_all_3;
		ez+=(zq*nz)*oneover_sqrt_all_3;

	}
	pot*=hs;
	ex*=hs;
	ey*=hs;
	ez*=hs;
}
*/
void calcD3triangle::triangleint_pot_exyz(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm,double &pot,double &ex,double &ey,double &ez){
	pot=0;
	ex=0;
	ey=0;
	ez=0;
	double hs=sqrt(sqr(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) + sqr(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + sqr(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3));
	double sum=0;
	int i;
	for(i=0;i<n;i++){
		double zeta,eta2,xi2;
			
 		eta2=eta[i];
		xi2=xi[i];
		zeta=1-xi2-eta2;

		double xq=x1*zeta+x2*xi2+x3*eta2-x0;
		double yq=y1*zeta+y2*xi2+y3*eta2-y0;
		double zq=z1*zeta+z2*xi2+z3*eta2-z0;

		double sqr_xq=sqr(xq);
		double sqr_yq=sqr(yq);
		double sqr_zq=sqr(zq);
		
		double oneover_sqrt_all=1./sqrt(sqr_xq+sqr_yq+sqr_zq);
		double oneover_sqrt_all_3=oneover_sqrt_all*oneover_sqrt_all*oneover_sqrt_all*w[i];
		oneover_sqrt_all*=w[i];

		pot+=oneover_sqrt_all;//speedup idea put pot in array of sizeof(w) and ausklammere w[i]

		ex+=(xq)*oneover_sqrt_all_3;
		ey+=(yq)*oneover_sqrt_all_3;
		ez+=(zq)*oneover_sqrt_all_3;

	}
	pot*=hs;
	ex*=hs;
	ey*=hs;
	ez*=hs;
//	if(inversenorm){
//		ex=-ex;
//		ey=-ey;
//		ez=-ez;
//	}
}

/*
double calcD3triangle::triangleintall(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double &G3D,double &dG3Ddx,double &dG3Ddy,double &dG3Ddz,double &dG3Ddxdx,double &dG3Ddxdy,double &dG3Ddxdz,double &dG3Ddydy,double &dG3Ddydz,double &dG3Ddzd){
	double hs=sqrt(sqr(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) + sqr(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + sqr(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3));
	double sum=0;
	int i;
	for(i=0;i<n;i++){
		double zeta,eta2,xi2;
			
 		eta2=eta[i];
		xi2=xi[i];
		zeta=1-xi2-eta2;
		sum+=fn(xp,yp,zp,x1*zeta+x2*xi2+x3*eta2,y1*zeta+y2*xi2+y3*eta2,z1*zeta+z2*xi2+z3*eta2)*w[i];
	}
	return sum*hs;
	//return sum*hs/2.;//use this if you use trianglequad
};
*/
//x0 is the point where the charge is located
//x1 is the point of field calculation



double calcD3triangle::G3D(double x1,double y1,double z1,double x0,double y0,double z0){
	return 1./(sqrt(sqr(x1-x0)+sqr(y1-y0)+sqr(z1-z0)));
};
double calcD3triangle::dG3Ddx(double x1,double y1,double z1,double x0,double y0,double z0){
	return (x1-x0)/(pow(sqr(x1-x0)+sqr(y1-y0)+sqr(z1-z0),3./2.));
}; 
double calcD3triangle::dG3Ddy(double x1,double y1,double z1,double x0,double y0,double z0){
	return (y1-y0)/(pow(sqr(x1-x0)+sqr(y1-y0)+sqr(z1-z0),3./2.));
};
double calcD3triangle::dG3Ddz(double x1,double y1,double z1,double x0,double y0,double z0){
	return (z1-z0)/(pow(sqr(x1-x0)+sqr(y1-y0)+sqr(z1-z0),3./2.));
};

void calcD3triangle::gauleg(double x1, double x2, double x[], double w[], int n) //from nr c p152 c4-5.pdf
//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
//Legendre n-point quadrature formula.
{
	x=x-1;
	w=w-1;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1; //High precision is a good idea for this routine.
	m=(n+1)/2; //The roots are symmetric in the interval, so
				//we only have to find half of them. 
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {// Loop over the desired roots.
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		//Starting with the above approximation to the ith root, we enter the main loop of
		//refinement by Newton's method.
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {// Loop up the recurrence relation to get the
								//Legendre polynomial evaluated at z. 
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			//p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			//by a standard relation involving also p2, the polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp; //Newton's method.
		}while(fabs(z-z1) > EPS);
		x[i]=xm-xl*z; //Scale the root to the desired interval,
		x[n+1-i]=xm+xl*z;// and put in its symmetric counterpart.
		w[i]=2.0*xl/((1.0-z*z)*pp*pp); //Compute the weight
		w[n+1-i]=w[i]; //and its symmetric counterpart.
	}
};

void calcD3triangle::gauslegsimplextriangle(double xi[],double eta[],double c[],int n){
	double *x=new double[n];
	double *w=new double[n];

	gauleg(-1.,1., x , w, n );
	int i,j,k=0;
	for(j=0;j<n;j++){
		for(i=0;i<n;i++){
			c[k]=(1.-x[i])/8.*w[i]*w[j];
			xi[k]=(1.+x[i])/2.;
			eta[k]=(1.-x[i])*(1+x[j])/4.;
			k++;
		}
	}
	delete[] x;
	delete[] w;
};

double calcD3triangle::G3DdnAnalytic(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm){
	//first find out if the point of evaluation is in the plane of the triangle or outside the plane
	//therefore calculate norm of triangle

	
	double px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,nx,ny,nz,n;
	double a,b,c,psi,cospsi,sinpsi;

	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	a=sqrt(sqr(ax)+sqr(ay)+sqr(az));
	b=sqrt(sqr(bx)+sqr(by)+sqr(bz));
	c=sqrt(sqr(cx)+sqr(cy)+sqr(cz));
	if(abs(a)<1e-13) return 0;
	if(abs(b)<1e-13) return 0;
	if(abs(c)<1e-13) return 0;
	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
	if(1.-abs(cospsi)<1e-13) return 0;
	sinpsi=sqrt(1.-sqr(cospsi));
	psi=acos_safe(cospsi);
	px=-x1+xp;py=-y1+yp;pz=-z1+zp;

	// b x a
	nx=az*by - ay*bz;
	ny=ax*bz - az*bx;
	nz=ay*bx - ax*by;
	n=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
	//if(!inversenorm) n=-n;

	nx/=n;ny/=n;nz/=n;

	

	double scalarprod=px*nx+py*ny+pz*nz;

	bool onplane=(abs(scalarprod)<1e-13);
	if(onplane) return 0;
	
		double rp=0;
	if(!onplane){
		rp=scalarprod;
		xp=xp-rp*nx;
		yp=yp-rp*ny;
		zp=zp-rp*nz;
		rp=abs(rp);
	}
	double Trianglei12;
	double Trianglei23;
	double Trianglei13;
	
	double ri1x=xp-x1,ri1y=yp-y1,ri1z=zp-z1;//definition different to paper Appl.Math.Moedel 1989_vol13_450
	double ri1=sqrt(sqr(x1-xp)+sqr(y1-yp)+sqr(z1-zp));
	if(abs(ri1)<1e-13) return (inversenorm?-1:1)*scalarprod*GetIntL1divR3(xp,yp,zp,x2,y2,z2,x3,y3,z3,rp);
	double cospsi_i=(ri1x*bx+ri1y*by+ri1z*bz)/(ri1*b);
	if(cospsi_i>1.) cospsi_i=1;
	if(cospsi_i<-1.) cospsi_i=-1;

	double sinpsi_i=sqrt(1.-sqr(cospsi_i));
	double psi_i=acos_safe(cospsi_i);
	//negate psi_i if (b x a) . (  b x ri1)<=0
	//n=b x a
	double ri1xb_x=-(-(ri1z*by) + ri1y*bz);
	double ri1xb_y=-(ri1z*bx - ri1x*bz);
	double ri1xb_z=-(-(ri1y*bx) + ri1x*by);

	double scalarprod2=ri1xb_x*nx+ri1xb_y*ny+ri1xb_z*nz;
	if(scalarprod2<0.) {
		psi_i=-psi_i;
		sinpsi_i=-sinpsi_i;
	}



	double Ui=ri1*sin(psi-psi_i)/b/sinpsi;
	double Vi=ri1*sinpsi_i/a/sinpsi;
	double Wi=1.-Vi-Ui;
	double n1,n2,n3;
	if(abs(Wi)<1e-13) {
		n1=0;
		Trianglei23=0;
	}
	else{
		n1=(Wi>=0)?1:-1;
		Trianglei23=GetIntL1divR3(xp,yp,zp,x2,y2,z2,x3,y3,z3,rp);
	}
	
	if(abs(Ui)<1e-13) {
		n2=0;
		Trianglei13=0;
	}
	else{
		n2=(Ui>=0)?1:-1;
		Trianglei13=GetIntL1divR3(xp,yp,zp,x1,y1,z1,x3,y3,z3,rp);
	}

	if(abs(Vi)<1e-13) {
		n3=0;
		Trianglei12=0;
	}
	else{
		n3=(Vi>=0)?1:-1;
		Trianglei12=GetIntL1divR3(xp,yp,zp,x1,y1,z1,x2,y2,z2,rp);
	}
	return (inversenorm?-1:1)*scalarprod*(n1*Trianglei23+n2*Trianglei13+n3*Trianglei12);
}
double calcD3triangle::G3Danalytic(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm){
	//first find out if the point of evaluation is in the plane of the triangle or outside the plane
	//therefore calculate norm of triangle
//undokument if numerical sollution should be enforced
	return calc.triangleint(calcD3triangle::G3D,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3);
	
	double px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,nx,ny,nz,n;
	double a,b,c,psi,cospsi,sinpsi;

	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	a=sqrt(sqr(ax)+sqr(ay)+sqr(az));
	b=sqrt(sqr(bx)+sqr(by)+sqr(bz));
	c=sqrt(sqr(cx)+sqr(cy)+sqr(cz));
	if(abs(a)<1e-13) return 0;
	if(abs(b)<1e-13) return 0;
	if(abs(c)<1e-13) return 0;
	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
    if(1.-abs(cospsi)<1e-13) return 0;
    sinpsi=sqrt(1.-sqr(cospsi));
	psi=acos_safe(cospsi);
	px=-x1+xp;py=-y1+yp;pz=-z1+zp;

	// b x a
	nx=az*by - ay*bz;
	ny=ax*bz - az*bx;
	nz=ay*bx - ax*by;
	n=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
	//if(!inversenorm) n=-n;

	nx/=n;ny/=n;nz/=n;

	

	double scalarprod=px*nx+py*ny+pz*nz;
	bool onplane=(abs(scalarprod)<1e-13);
	//if(!onplane) return calc.triangleint(calcD3triangle::G3D,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3);
	double rp=0;
	if(!onplane){
		rp=scalarprod;
		xp=xp-rp*nx;
		yp=yp-rp*ny;
		zp=zp-rp*nz;
		rp=abs(rp);
	}
	double Trianglei12;
	double Trianglei23;
	double Trianglei13;
	
	double ri1x=xp-x1,ri1y=yp-y1,ri1z=zp-z1;//definition different to paper Appl.Math.Moedel 1989_vol13_450
	double ri1=sqrt(sqr(x1-xp)+sqr(y1-yp)+sqr(z1-zp));
	if(abs(ri1)<1e-13) return onplane?GetIntL1(xp,yp,zp,x2,y2,z2,x3,y3,z3):GetIntL1OutOfPlane(xp,yp,zp,x2,y2,z2,x3,y3,z3,rp);
	double cospsi_i=(ri1x*bx+ri1y*by+ri1z*bz)/(ri1*b);
	if(cospsi_i>1) cospsi_i=1;
	if(cospsi_i<-1) cospsi_i=-1;

	double sinpsi_i=sqrt(1.-sqr(cospsi_i));
	double psi_i=acos_safe(cospsi_i);
	//negate psi_i if (b x a) . (  b x ri1)<=0
	//n=b x a
	double ri1xb_x=-(-(ri1z*by) + ri1y*bz);
	double ri1xb_y=-(ri1z*bx - ri1x*bz);
	double ri1xb_z=-(-(ri1y*bx) + ri1x*by);

	double scalarprod2=ri1xb_x*nx+ri1xb_y*ny+ri1xb_z*nz;
	if(scalarprod2<0.) {
		psi_i=-psi_i;
		sinpsi_i=-sinpsi_i;
	}



	double Ui=ri1*sin(psi-psi_i)/b/sinpsi;
	double Vi=ri1*sinpsi_i/a/sinpsi;
	double Wi=1.-Vi-Ui;
	double n1,n2,n3;
	if(abs(Wi)<1e-13) {
		n1=0;
		Trianglei23=0;
	}
	else{
		n1=(Wi>=0)?1:-1;
		Trianglei23=onplane?GetIntL1(xp,yp,zp,x2,y2,z2,x3,y3,z3):GetIntL1OutOfPlane(xp,yp,zp,x2,y2,z2,x3,y3,z3,rp);
	}
	
	if(abs(Ui)<1e-13) {
		n2=0;
		Trianglei13=0;
	}
	else{
		n2=(Ui>=0)?1:-1;
		Trianglei13=onplane?GetIntL1(xp,yp,zp,x1,y1,z1,x3,y3,z3):GetIntL1OutOfPlane(xp,yp,zp,x1,y1,z1,x3,y3,z3,rp);
	}

	if(abs(Vi)<1e-13) {
		n3=0;
		Trianglei12=0;
	}
	else{
		n3=(Vi>=0)?1:-1;
		Trianglei12=onplane?GetIntL1(xp,yp,zp,x1,y1,z1,x2,y2,z2):GetIntL1OutOfPlane(xp,yp,zp,x1,y1,z1,x2,y2,z2,rp);
	}
	return (n1*Trianglei23+n2*Trianglei13+n3*Trianglei12);

	
}



//D3element needs list of sollutions for each 0 ...0 1 0...0 electrode voltage configuration
//new electrode container class
//new world container class
//MatlabDisplay needs other depth than real element depth
void D3element::correctNorm(double x0,double y0,double z0){
	int	n=GetAmountOfSubelements();
	PD3element *el=new PD3element[n];
	int cnt=0;
	GetListOfBaseElements(el,cnt);
	int row,col,col2;
	int elcnt;
	double a[3],b[3],c[3],d[3],color[3],xdir,ydir,zdir,mul;	
	
	for(col2=0;col2<n;col2++){//geht alle Flächenelemente durch
		el[col2]->GetReferencePoint(xdir,ydir,zdir);
		bool towardsx0y0z0=true;
		for(col=0;col<n;col++){//geht alle Flächenelemente durch
			if(col!=col2 && el[col]->IntersectWithRay(x0,y0,z0,xdir-x0,ydir-y0,zdir-z0,mul))
				if((mul<1.) && (mul>0.)) 
					towardsx0y0z0=!towardsx0y0z0;
		}
		el[col2]->SetNormTowards(x0,y0,z0,towardsx0y0z0);
	}
	delete[] el;
}

double D3element::GetSelfPotential(){return 0;};
double D3element::GetPotentialAt(double x,double y,double z){return 0;};
void D3element::GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez){pot=0;ex=0;ey=0;ez=0;}
double D3element::GetSelfDoubleLayerPotential(){return 0;};
double D3element::GetDoubleLayerPotentialAt(double x,double y,double z){return 0;};
D3element::D3element(	double x_,double y_,double z_,
			double phi_,double theta_,double psi_,
			bool inversenorm_,bool refineable_,D3element* parent_,double epsilon,string name):
				isBaseElement(false),x(x_),y(y_),z(z_),
				phi(-phi_/180.*pi),theta(-theta_/180.*pi),psi(-psi_/180*pi),
				inversenorm(inversenorm_),
				refineable(refineable_),parent(parent_),epsilon(epsilon),Color(kBlue),name(name){};
D3element::~D3element(){
	if (!isBaseElement){
		deleteSubelements();
	}
};

string D3element::getName(){
	return name;
}


double D3element::GetArea(){
	return 0;
}

void D3element::GetCenter(double &x,double &y,double &z){
	x=0;
	y=0;
	z=0;
}

bool D3element::GetInversenorm(){
	return inversenorm;
}

void D3element::setName(string name2){
	name=name2;
}

bool D3element::IntersectWithRay(double x0,double y0,double z0,double xdir,double ydir,double zdir,double &mul){return true;}

void D3element::SetNormTowards(double x,double y,double z,bool towards){}

void D3element::Add(PD3element el){subelement.push_back(el);};
void D3element::GetReferencePoint(double &x,double &y, double &z){};
void D3element::init(bool ignorefirst){
	PD3element el=this;
	/*TGeoRotation r;
	TGeoTranslation t;
	TGeoHMatrix id;
	h=id;*/
	if(ignorefirst) {//if element is a Rectangle generated by Constructor with string as first parameter
		/*r.SetAngles(-el->psi*180./pi,-el->theta*180./pi,-el->phi*180./pi);
		t.SetTranslation(el->x,el->y,el->z);
		h*=r;
		h*=t;*/
		el=el->parent;
	}
	while(el!=NULL){
		/*r.SetAngles(-el->psi*180./pi,-el->theta*180./pi,-el->phi*180./pi);
		t.SetTranslation(el->x,el->y,el->z);
		h*=r;
		h*=t;*/

		rotate(el->phi,el->theta,el->psi);
		shift(el->x,el->y,el->z);
		el=el->parent;
	};
}

void D3element::refine(double length){//TODO change refine to refine parts not in power of 2 but for any number
	if(refineable){
		if(!isBaseElement) {//rollback refinement first
			deleteSubelements();
		}
		if(length>0){
			//Do refine here
			createNewSubelements(length);
		}
	}
	else{//!refineable
		if(isBaseElement) return;//This baseElement cannot be refined
		else{//!isBaseElement then refine Subelements
			PD3elementList::iterator i;
			for(i=subelement.begin();i!=subelement.end();++i) 
				(*i)->refine(length);
		}
	}
}

void D3element::createNewSubelements(double length){}

void D3element::refine(double length,int num){//This is needed to skip refinement under world of first num electrodes. eg for bounding box
	PD3elementList::iterator i;
	i=subelement.begin();
	for(;num>0;num--) i++;
	(*i)->refine(length);
};
int D3element::GetAmountOfSubelements(){
	if (isBaseElement) return 1;
	int cnt=0;
	PD3elementList::iterator i;
	for(i=subelement.begin();i!=subelement.end();++i){
		cnt+=(*i)->GetAmountOfSubelements();
	}
	return cnt;
};

int D3element::PrintAmountOfSubelements(){
	if (isBaseElement) {
		cout <<"1\n";
		return 1;
	}
	int cnt=0;
	PD3elementList::iterator i;
	int elcnt;
	for(i=subelement.begin();i!=subelement.end();++i){
		PD3element p=*i;
		elcnt=p->GetAmountOfSubelements();
		cnt+=elcnt;
		cout <<"class: "<<typeid(*(p)).name()<<" name:"<<(*(p)).getName()<<" "<<elcnt<<"\n";
	}
	return cnt;
}

void D3element::GetListOfBaseElements(PD3element *el,int &cnt){
	if (isBaseElement) {
		el[cnt++]=this;
		return;
	}
	PD3elementList::iterator i;
	for(i=subelement.begin();i!=subelement.end();++i){
		(*i)->GetListOfBaseElements(el,cnt);
	}
};
void D3element::rotate(double phi_,double theta_,double psi_){
int i=1;
};
void D3element::shift(double xs,double ys,double zs){
int iii=1;
};
void D3element::rotate2(double phi_,double theta_,double psi_,double x,double y,double z,double &x2,double &y2,double &z2){
	double cosphi=cos(phi_);
	double sinphi=sin(phi_);
	double costheta=cos(theta_);
	double sintheta=sin(theta_);
	double cospsi=cos(psi_);
	double sinpsi=sin(psi_);

	x2=(cospsi*cosphi-costheta*sinphi*sinpsi)*x		+	(cospsi*sinphi+costheta*cosphi*sinpsi)*y	+	sinpsi*sintheta*z;
	y2=(-sinpsi*cosphi-costheta*sinphi*cospsi)*x	+	(-sinpsi*sinphi+costheta*cosphi*cospsi)*y	+	cospsi*sintheta*z;
	z2=(sintheta*sinphi)*x							+	(-sintheta*cosphi)*y						+	costheta*z;
	

	//x2=(cospsi*cosphi-costheta*sinphi*sinpsi)*x		+	(-sinpsi*cosphi-costheta*sinphi*cospsi)*y		+ (sintheta*sinphi)*z;
	//y2=	(cospsi*sinphi+costheta*cosphi*sinpsi)*x	+	(-sinpsi*sinphi+costheta*cosphi*cospsi)*y		+(-sintheta*cosphi)*z;
	//z2=	sinpsi*sintheta*x							+   cospsi*sintheta*y								+costheta*z;

	//als
	//x2=(cospsi*cosphi-costheta*sinphi*sinpsi)*x		+	(-sinpsi*cosphi-costheta*sinphi*cospsi)*y	+	sintheta*sinphi*z;
	//y2=(cospsi*sinphi+costheta*cosphi*sinpsi)*x		+	(-sinpsi*sinphi+costheta*cosphi*cospsi)*y	-	sintheta*cosphi*z;
	//z2=(sintheta*sinpsi)*x							+	(sintheta*cospsi)*y							+	costheta*z;

};

void D3element::rotateinv(double phi_,double theta_,double psi_,double x,double y,double z,double &x2,double &y2,double &z2){
	double cosphi=cos(phi_);
	double sinphi=sin(phi_);
	double costheta=cos(theta_);
	double sintheta=sin(theta_);
	double cospsi=cos(psi_);
	double sinpsi=sin(psi_);
	//Attention constructor of D3element inverts sign of all angles

	x2=(cospsi*cosphi-costheta*sinphi*sinpsi)*x		+	(-sinpsi*cosphi-costheta*sinphi*cospsi)*y		+ (sintheta*sinphi)*z;
	y2=	(cospsi*sinphi+costheta*cosphi*sinpsi)*x	+	(-sinpsi*sinphi+costheta*cosphi*cospsi)*y		+(-sintheta*cosphi)*z;
	z2=	sinpsi*sintheta*x							+   cospsi*sintheta*y								+costheta*z;

	//als
	//x2=(cospsi*cosphi-costheta*sinphi*sinpsi)*x		+	(-sinpsi*cosphi-costheta*sinphi*cospsi)*y	+	sintheta*sinphi*z;
	//y2=(cospsi*sinphi+costheta*cosphi*sinpsi)*x		+	(-sinpsi*sinphi+costheta*cosphi*cospsi)*y	-	sintheta*cosphi*z;
	//z2=(sintheta*sinpsi)*x							+	(sintheta*cospsi)*y							+	costheta*z;

};

//int D3element::Get3DMatlab( mxArray *plhs[],int cnt){
//	if(cnt==-1){
//		int subel=GetAmountOfSubelements();
//		plhs[0]=mxCreateDoubleMatrix(3,subel, mxREAL);		
//		plhs[1]=mxCreateDoubleMatrix(3,subel, mxREAL);
//		plhs[2]=mxCreateDoubleMatrix(3,subel, mxREAL);
//		plhs[3]=mxCreateDoubleMatrix(3,subel, mxREAL);
//		cnt=0;
//	}
//	if(!isBaseElement){
//		PD3elementList::iterator i;
//		for(i=subelement.begin();i!=subelement.end();++i){
//			cnt=(*i)->Get3DMatlab(plhs,cnt);
//		}
//		return cnt;
//	}
//	else{
//		double *A=mxGetPr(plhs[0])+cnt*3;
//		double *B=mxGetPr(plhs[1])+cnt*3;
//		double *C=mxGetPr(plhs[2])+cnt*3;
//		double *COL=mxGetPr(plhs[3])+cnt*3;
//		GetTriangle(A,B,C,COL);
//		cnt++;
//		return cnt;
//	}
//};

void D3element::GetTriangle(double *A,double *B,double *C,double *COL){};
void D3element::GetRectangle(double *A,double *B,double *C,double *D,double *COL){};
void D3element::deleteSubelements(){ 
	PD3elementList::iterator i;
	for(i=subelement.begin();i!=subelement.end();++i) delete *i;
	subelement.clear();
	isBaseElement=true;
};


/////////////////////////////////////////////////////////////////////////////////////////////
D3electrode::D3electrode():voltage(0),totalrot(-1),D3element(0,0,0,0,0,0,false,false,NULL){
	
};

//D3electrode::D3electrode(const D3electrode &cp):voltage(0),
//				D3element(0,0,0,0,0,0,false,false,NULL){}
typedef D3electrode* PD3electrode;
void D3electrode::SetSymmetryWith(bool xsym,bool ysym,bool zsym,bool diagmirror,D3electrode *symel2,int rotsym,int totalrot2){
	if(totalrot<0) {
		totalrot=totalrot2;
		int siz=(totalrot<<4)+15+1;
		symel=new PD3electrode[siz];
		symel[0]=this;
		for(int i=1;i<siz;++i) 
			symel[i]=NULL;
	}
	symel[(xsym?1:0)+(ysym?2:0)+(zsym?4:0)+(diagmirror?8:0)+(rotsym<<4)]=symel2;
}
void D3electrode::SetSymmetry(bool xsym,bool ysym,bool zsym,bool diagmirror,int rotsym,int totalrot2){
	SetSymmetryWith(xsym,ysym,zsym,diagmirror,this,rotsym,totalrot2);
}
D3electrode *D3electrode::GetSymmetryWith(int num){
	return symel[num];
}


void D3electrode::SetCardinalNumber(int num){
	cardinalnum=num;
}
int D3electrode::GetCardinalNumber(){
	return cardinalnum;
}
D3electrode::~D3electrode(){

};

void D3electrode::insert(PD3element el){
	el->parent=this;
	subelement.push_back(el);//triangle at the left bottom corner
	if(el->getName()!=""){;
		name=name+el->getName()+";";
	}
};

void D3electrode::SetVoltage(double voltage_){voltage=voltage_;};
double D3electrode::GetVoltage(){return voltage;};
int D3electrode::GetTotalrot(){return totalrot;};

//////////////////////////////////////////////////////////////////////////////////////////////////

//                
//               /\
//              /  \
//             /    \
//         ^  /      \
//        x3 /________\
//          / \      / \
//         /   \    /   \
//        /     \  /     \
//       /_______\/_______\
//      xyz     x2>
//
//phi,theta,psi are euler coordinates. rotation point is (x,y,z)
calcD3triangle calc;
void D3triangle::GetReferencePoint(double &x,double &y, double &z){x=(x1+x2+x3)/3.,y=(y1+y2+y3)/3.;z=(z1+z2+z3)/3.;};
bool D3triangle::IntersectWithRay(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double dirx,double diry,double dirz,double &mul){
	// a=xyz2-xyz1, b=xyz3-xyz1;
	// c=xyz0-xyz1;
	// n=a x b
	// Eqation to solve: c+mul*dir == amul*a+bmul*b
	// |dir . n| * a = sign(dir . n) * dir . (c x b)
	// |dir . n| * b = sign(dir . n) * dir . (a x c)
	// |dir . n| * mul = -sign(dir . n) * (c . n)
	
	
	// a x b
	double nx,ny,nz;
	nx=-az*by + ay*bz;

	ny=-ax*bz + az*bx;
	nz=-ay*bx + ax*by;

	double dirDOTn=dirx*nx+diry*ny+dirz*nz;
	double signDirDOTn;

	if(dirDOTn>0.) signDirDOTn=1.;
	else if(dirDOTn<0.) {
		signDirDOTn=-1.;
		dirDOTn=-dirDOTn;
	}
	else return false; //Ray lies fully on plane
	
	// c x b
	double cXb_x=-cz*by + cy*bz;
	double cXb_y=-cx*bz + cz*bx;
	double cXb_z=-cy*bx + cx*by;
	
	double dirDOTcXb=signDirDOTn*(dirx*cXb_x+diry*cXb_y+dirz*cXb_z);

	if (dirDOTcXb>0.0){
		// a x c
		double aXc_x=-az*cy + ay*cz;
		double aXc_y=-ax*cz + az*cx;
		double aXc_z=-ay*cx + ax*cy;
		
		double dirDOTaXc=signDirDOTn*(dirx*aXc_x+diry*aXc_y+dirz*aXc_z);
		if(dirDOTaXc>0.0){
			if(dirDOTcXb+dirDOTaXc<=dirDOTn){
				double cDOTn=-signDirDOTn*(cx*nx+cy*ny+cz*nz);
				if(cDOTn>=0.){
					mul=cDOTn/dirDOTn;
					return true;
				}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
			}
		}
	}
	mul=-1;
	return false;
}
void D3triangle::createNewSubelements(double length){
	if(stripeit) Stripeit(length);
	else Newtriangles(length);
}

void D3triangle::Newtriangles(double length){
	double la=sqrt(xa*xa+ya*ya);
	double lb=sqrt(xb*xb+yb*yb);
	double lc=sqrt(sqr(xa-xb)+sqr(ya-yb));


	
	int na=la/length;
	int nb=lb/length;
	int nc=lb/length;
	int n;
	if     (na>nb && na>nc) n=na;
	else if(nb>na && nb>nc) n=nb;
	else n=nc;
	if(n<=1) {
		isBaseElement=true;
		return;
	}
		
	isBaseElement=false;

	double xadiv=xa/double(n);
	double yadiv=ya/double(n);
	double xbdiv=xb/double(n);
	double ybdiv=yb/double(n);

	int i,j;
	for(j=0;j<n;++j){
		for(i=0;i<(n-j);++i){
			subelement.push_back(        new D3triangle(xadiv,yadiv,xbdiv,ybdiv,xadiv*double(i)+xbdiv*double(j),yadiv*double(i)+ybdiv*double(j),0,0,0,0,inversenorm,this));
			if(j>0) subelement.push_back(new D3triangle(-xadiv,-yadiv,-xbdiv,-ybdiv,xadiv*double(i+1)+xbdiv*double(j),yadiv*double(i+1)+ybdiv*double(j),0,0,0,0,inversenorm,this));
		}
	}



}

double D3triangle::GetArea(){
	double nz;
	nz=-ya*xb + xa*yb;
	if (nz<0) nz=-nz;
	return 0.5*nz;
}

void D3triangle::GetCenter(double &xx,double &yy,double &zz){
	xx=(x1+x2+x3)/3.;
	yy=(y1+y2+y3)/3.;
	zz=(z1+z2+z3)/3.;
}



void D3triangle::Stripeit(double length){
	double la=sqrt(xa*xa+ya*ya);
	double lb=sqrt(xb*xb+yb*yb);
	double lc=sqrt(sqr(xa-xb)+sqr(ya-yb));
	if(length>=la && length>=lb && length>=lc){
		isBaseElement=true;	
		return; //do no refinement
	}
	double x0_,y0_,xa_,ya_,xb_,yb_;
	//First find smallest side and take coordinates in such a way that xy0_ is at the sharpest edge
	//we want to stripe the triangle starting from the sharpest edge
	if(la>=lc && lb>=lc){
		x0_=0;y0_=0;
		xa_=xa;ya_=ya;
		xb_=xb;yb_=yb;
	}
	else if(la>=lb && lc>=lb){
		x0_=xa;y0_=ya;
		xb_=-xa;yb_=-ya;
		xa_=xb-xa;ya_=yb-ya;
	}
	else{//lb>la && lc>la
		x0_=xb;y0_=yb;
		xa_=-xb;ya_=-yb;
		xb_=xa-xb;yb_=ya-yb;		
	}
	double la_=sqrt(xa_*xa_+ya_*ya_);
	double lb_=sqrt(xb_*xb_+yb_*yb_);
	double lc_=sqrt(sqr(xa_-xb_)+sqr(ya_-yb_));
	//the triangle has a very great angle on one of the sides so first subdivide it in two parts
	if(lc_/min(la_,lb_)>0.7 && min(la_,lb_)/max(la_,lb_)<0.6){
		isBaseElement=false;
		if(la_>lb_){
			D3element* pD3triangle1=new D3triangle(xa_*lb_/la_,ya_*lb_/la_,xb_,yb_, x0_,y0_,0, 0,0,0,inversenorm,this);
			pD3triangle1->createNewSubelements(length);
			subelement.push_back(pD3triangle1);
			D3element* pD3triangle2=new D3triangle(xa_*(1.-lb_/la_),ya_*(1.-lb_/la_),xb_-xa_*lb_/la_,yb_-ya_*lb_/la_, x0_+xa_*lb_/la_,y0_+ya_*lb_/la_,0, 0,0,0,inversenorm,this);
			pD3triangle2->createNewSubelements(length);
			subelement.push_back(pD3triangle2);
		}
		else{
			D3element* pD3triangle1=new D3triangle(xa_,ya_,xb_*la_/lb_,yb_*la_/lb_, x0_,y0_,0, 0,0,0,inversenorm,this);
			pD3triangle1->createNewSubelements(length);
			subelement.push_back(pD3triangle1);
			D3element* pD3triangle2=new D3triangle(xa_-xb_*la_/lb_,ya_-yb_*la_/lb_,xb_*(1.-la_/lb_),yb_*(1.-la_/lb_), x0_+xb_*la_/lb_,y0_+yb_*la_/lb_,0, 0,0,0,inversenorm,this);
			pD3triangle2->createNewSubelements(length);
			subelement.push_back(pD3triangle2);
		}
		return;
	}

	int div1=ceil(max(la_,lb_)/length);
	int n;
	//This is the triangle at the sharpest edge
	subelement.push_back(new D3triangle(xa_/div1,ya_/div1,xb_/div1,yb_/div1, x0_,y0_,0, 0,0,0,inversenorm,this));
	for(n=1;n<div1;n++){
		//make rectangular stripes and refine them accordingly
		D3element* pD3rectangle=new D3rectangle(xa_/div1, ya_/div1,
			(xb_/div1)*(n+1)-(xa_/div1)*n,(yb_/div1)*(n+1)-(ya_/div1)*n,
			(xb_/div1)*n-(xa_/div1)*n,(yb_/div1)*n-(ya_/div1)*n,
			x0_+ (xa_/div1)*n,y0_+(ya_/div1)*n,0, 0,0,0,inversenorm,this);
		pD3rectangle->createNewSubelements(length);
		subelement.push_back(pD3rectangle);
	}
	isBaseElement=false;
	return;
}

bool D3triangle::IntersectWithRay(double x0,double y0,double z0,double dirx,double diry,double dirz,double &mul){
	double ax,ay,az,bx,by,bz,cx,cy,cz;
	ax=x2-x1;
	ay=y2-y1;
	az=z2-z1;
	
	bx=x3-x1;
	by=y3-y1;
	bz=z3-z1;

	cx=x0-x1;
	cy=y0-y1;
	cz=z0-z1;
	return IntersectWithRay(ax,ay,az,bx,by,bz,cx,cy,cz,dirx,diry,dirz,mul);
}
              
double D3triangle::GetSelfPotential(){
	double xp,yp,zp,L11,L12,L13;
	GetReferencePoint(xp,yp,zp);
	L11=calc.GetIntL1(xp,yp,zp,x1,y1,z1,x2,y2,z2);
	L12=calc.GetIntL1(xp,yp,zp,x2,y2,z2,x3,y3,z3);
	L13=calc.GetIntL1(xp,yp,zp,x3,y3,z3,x1,y1,z1);
	return (L11+L12+L13);
};
double D3triangle::GetPotentialAt(double x,double y,double z){
	return calc.G3Danalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm);
};

void D3triangle::GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez){
	calc.triangleint_pot_exyz(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm,pot,ex,ey,ez);
}

double D3triangle::GetSelfDoubleLayerPotential(){
	return 0;
};
double D3triangle::GetDoubleLayerPotentialAt(double x,double y,double z){
	//return calc.G3DdnAnalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm);//undocument if nonanalytical sollution is wanted

	double ax,ay,az,bx,by,bz,nx,ny,nz,nasq,nbsq,n,dG3Ddx,dG3Ddy,dG3Ddz;
		ax=x2-x1;ay=y2-y1;az=z2-z1;
		bx=x3-x1;by=y3-y1;bz=z3-z1;
	
		
		
		nx=(-(az*by) + ay*bz);
		ny=(az*bx - ax*bz);
		nz=(-(ay*bx) + ax*by);
		n=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
		
		if(inversenorm) n=-n;

		dG3Ddx=calc.triangleint(calcD3triangle::dG3Ddx,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);
		dG3Ddy=calc.triangleint(calcD3triangle::dG3Ddy,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);
		dG3Ddz=calc.triangleint(calcD3triangle::dG3Ddz,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);

		return (nx*dG3Ddx+ny*dG3Ddy+nz*dG3Ddz)/n;


};

D3triangle::D3triangle(	double xa_,double ya_,double xb_,double yb_,
			double x,double y,double z,
			double phi_,double theta_,double psi_,bool inversenorm,PD3element parent,bool stripeit):
				xa(xa_),ya(ya_),xb(xb_),yb(yb_),
				x1(0) ,y1(0) ,z1(0),
				x2(xa_),y2(ya_),z2(0),
				x3(xb_),y3(yb_),z3(0),stripeit(stripeit),
				D3element(x,y,z,phi_,theta_,psi_,inversenorm,true,parent,(abs(xa_)+abs(ya_))/100.){
					isBaseElement=true;					
					init();
				};

D3triangle::D3triangle(double x1_,double y1_,double z1_,double x2_,double y2_,double z2_,double x3_,double y3_,double z3_,
			bool inversenorm,PD3element parent,bool stripeit):x1(x1_),y1(y1_),z1(z1_),x2(x2_),y2(y2_),z2(z2_),x3(x3_),y3(y3_),z3(z3_),stripeit(stripeit),
					D3element(x1_,y1_,z1_,0,0,0,inversenorm,false,parent,0){
	isBaseElement=true;		
	refineable=true;
	double nx,ny,nz,nnorm,nxynorm,za,zb,zc,xa2,ya2,za2,xb2,yb2,zb2;
	xa2=x2-x1;xb2=x3-x1;
	ya2=y2-y1;yb2=y3-y1;
	za2=z2-z1;zb2=z3-z1;
	epsilon=(abs(xa2)+abs(ya2))/100.;
	//n= a x b
	nx=-za2*yb2 + ya2*zb2;
	ny=-xa2*zb2 + za2*xb2;
	nz=-ya2*xb2 + xa2*yb2;
	
	nnorm=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
	nx/=nnorm;
	ny/=nnorm;
	nz/=nnorm;
	nxynorm=sqrt(sqr(nx)+sqr(ny));
	theta=acos(nz);
	phi=0;
	if (nxynorm!=0){
		nx/=nxynorm;
		ny/=nxynorm;
		psi=acos(-ny);
		if (nx<0) psi=-psi;
	}
	else psi=0;
	psi=-psi;
	theta=-theta;

	rotateinv(phi,theta,psi,xa2,ya2,za2,xa,ya,za);
	rotateinv(phi,theta,psi,xb2,yb2,zb2,xb,yb,zb);
	
	init(true);
}

					
					
					
void D3triangle::SetNormTowards(double x0,double y0,double z0,bool towards){
	double nx,ny,nz,xa2,ya2,za2,xb2,yb2,zb2,x_,y_,z_;
	xa2=x2-x1;xb2=x3-x1;
	ya2=y2-y1;yb2=y3-y1;
	za2=z2-z1;zb2=z3-z1;
	
	//n= a x b
	nx=-za2*yb2 + ya2*zb2;
	ny=-xa2*zb2 + za2*xb2;
	nz=-ya2*xb2 + xa2*yb2;
	GetReferencePoint(x_,y_,z_);
	double nDOTxMINUSx0=nx*(x_-x0)+ny*(y_-y0)+nz*(z_-z0);
	inversenorm =((nDOTxMINUSx0<0) ^ towards);
}

D3triangle::~D3triangle(){};	
void D3triangle::rotate(double phi_,double theta_,double psi_){
	rotate2(phi_,theta_,psi_,x1,y1,z1,x1,y1,z1);
	rotate2(phi_,theta_,psi_,x2,y2,z2,x2,y2,z2);
	rotate2(phi_,theta_,psi_,x3,y3,z3,x3,y3,z3);
};
void D3triangle::shift(double xs,double ys,double zs){
	x1+=xs;y1+=ys;z1+=zs;
	x2+=xs;y2+=ys;z2+=zs;
	x3+=xs;y3+=ys;z3+=zs;
};

void D3triangle::GetTriangle(double *A,double *B,double *C,double *COL){

	if(inversenorm){
		A[0]=x3;A[1]=y3;A[2]=z3;
		B[0]=x2;B[1]=y2;B[2]=z2;
		C[0]=x1;C[1]=y1;C[2]=z1;
		COL[0]=0;COL[1]=0.5;COL[2]=0.7;
	}
	else{
		A[0]=x1;A[1]=y1;A[2]=z1;
		B[0]=x2;B[1]=y2;B[2]=z2;
		C[0]=x3;C[1]=y3;C[2]=z3;
		COL[0]=0;COL[2]=0.5;COL[1]=0.7;
	}
}


////////////////////////////////////////////////////////
#define TRIANGULATE_D3RECTANGLE
#undef TRIANGULATE_D3RECTANGLE
D3rectangle::D3rectangle(double xa_,double ya_,double xc_,double yc_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm,PD3element parent,bool subdivideInSquares):
					xa(xa_),ya(ya_),xb(xa_+xc_),yb(ya_+yc_),xc(xc_),yc(yc_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent,(abs(xa_)+abs(ya_))/100.){
		if(subdivideInSquares){
			//try to use squares and not long rectangles
			double lena=sqrt(sqr(xa)+sqr(ya));
			double lenc=sqrt(sqr(xc)+sqr(yc));
			if(lena<lenc){
				int subdivc=int(lenc/lena);
				int i;
				for(i=0;i<=subdivc-1;i++){
#ifdef TRIANGULATE_D3RECTANGLE
					subelement.push_back(new D3triangle(xa,ya,xa+xc/double(subdivc),ya+yc/double(subdivc),xc/double(subdivc)*double(i),yc/double(subdivc)*double(i),0,0,0,0,inversenorm,this));//triangle at the left bottom corner
					subelement.push_back(new D3triangle(xa+xc/double(subdivc),ya+yc/double(subdivc),xc/double(subdivc),yc/double(subdivc),xc/double(subdivc)*double(i),yc/double(subdivc)*double(i),0,0,0,0,inversenorm,this));//triangle at the right bottom corner
#else
					subelement.push_back(new D3rectangle(xa,ya,xc/double(subdivc),yc/double(subdivc),xc/double(subdivc)*double(i),yc/double(subdivc)*double(i),0,0,0,0,inversenorm,this,false));	
#endif
				}
			}
			else{
				int subdiva=int(lena/lenc);
				int i;
				for(i=0;i<=subdiva-1;i++){
#ifdef TRIANGULATE_D3RECTANGLE
					subelement.push_back(new D3triangle(xa/double(subdiva),ya/double(subdiva),xa/double(subdiva)+xc,ya/double(subdiva)+yc,xa/double(subdiva)*double(i),ya/double(subdiva)*double(i),0,0,0,0,inversenorm,this));//triangle at the left bottom corner
					subelement.push_back(new D3triangle(xa/double(subdiva)+xc,ya/double(subdiva)+yc,xc,yc,xa/double(subdiva)*double(i),ya/double(subdiva)*double(i),0,0,0,0,inversenorm,this));//triangle at the right bottom corner
#else
					subelement.push_back(new D3rectangle(xa/double(subdiva),ya/double(subdiva),xc,yc,xa/double(subdiva)*double(i),ya/double(subdiva)*double(i),0,0,0,0,inversenorm,this,false));
#endif
				}
			}
		}
		else{
#ifdef TRIANGULATE_D3RECTANGLE
			subelement.push_back(new D3triangle(xa,ya,xb,yb,0,0,0,0,0,0,inversenorm,this));//triangle at the left bottom corner
			subelement.push_back(new D3triangle(xb,yb,xc,yc,0,0,0,0,0,0,inversenorm,this));//triangle at the right bottom corner
#else
			isBaseElement=true;		
			refineable=true;
			x1=0;y1=0;z1=0;
			x2=xa;y2=ya;z2=0;
			x3=xc+xa;y3=yc+ya;z3=0;
			x4=xc;y4=yc;z4=0;
			init();
#endif
		}
}

double D3rectangle::GetArea(){
	double nz;
	nz=-ya*xb + xa*yb;
	if (nz<0) nz=-nz;
	double nz2;
	nz2=-yb*xc + xb*yc;
	if (nz2<0) nz2=-nz2;
	return 0.5*nz+0.5*nz2;
}


void D3rectangle::GetCenter(double &xx,double &yy,double &zz){
	xx=(x1+x2+x3+x4)/4.;
	yy=(y1+y2+y3+y4)/4.;
	zz=(z1+z2+z3+z4)/4.;
}

void D3rectangle::createNewSubelements(double length){
#ifdef TRIANGULATE_D3RECTANGLE
	isBaseElement=false;
	D3triangle *pD3triangle1=new D3triangle(xa,ya,xb,yb,0,0,0,0,0,0,inversenorm,this);
	pD3triangle1->createNewSubelements(length);
	subelement.push_back(pD3triangle1);//triangle at the left bottom corner
	D3triangle *pD3triangle2=new D3triangle(xb,yb,xc,yc,0,0,0,0,0,0,inversenorm,this);//triangle at the right bottom corner
	pD3triangle2->createNewSubelements(length);
	subelement.push_back(pD3triangle2);//triangle at the left bottom corner

#else
	double la=sqrt(sqr(xa)+sqr(ya));
	double lb=sqrt(sqr(xb-xa)+sqr(yb-ya));
	double lc=sqrt(sqr(xc-xb)+sqr(yc-yb));
	double ld=sqrt(sqr(xc)+sqr(yc));
	if((la<=length) && (lb <=length) && (lc <=length) && (ld<=length)){
		isBaseElement=true;
		return;
	}

	isBaseElement=false;

	int div1=ceil(max(la,lc)/length);
	int div2=ceil(max(lb,ld)/length);

	int n1,n2;
	double x0sumstart=0,y0sumstart=0;
	for(n2=0;n2<div2;n2++){
		double x0sum=xc/div2*n2,y0sum=yc/div2*n2;		
		for(n1=0;n1<div1;n1++){
			subelement.push_back(new D3rectangle(
				xa/div1*(div2-n2)/div2+(xb-xc)/div1*n2/div2,
				ya/div1*(div2-n2)/div2+(yb-yc)/div1*n2/div2,
				
				xc/div2*(div1-n1-1)/div1+(xb-xa)/div2*(n1+1)/div1+  xa/div1*(div2-n2)/div2+(xb-xc)/div1*n2/div2,
				yc/div2*(div1-n1-1)/div1+(yb-ya)/div2*(n1+1)/div1+  ya/div1*(div2-n2)/div2+(yb-yc)/div1*n2/div2,


				xc/div2*(div1-n1)/div1+(xb-xa)/div2*n1/div1,
				yc/div2*(div1-n1)/div1+(yb-ya)/div2*n1/div1,


			
				x0sum,y0sum,
				0,0,0,0,inversenorm,this));
			x0sum+=xa/div1*(div2-n2)/div2+(xb-xc)/div1*n2/div2;
			y0sum+=ya/div1*(div2-n2)/div2+(yb-yc)/div1*n2/div2;
		}
		
	}
#endif
}


D3rectangle::D3rectangle(double xa_,double ya_,double xb_,double yb_,double xc_,double yc_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm,PD3element parent):
					xa(xa_),ya(ya_),xb(xb_),yb(yb_),xc(xc_),yc(yc_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent,(abs(xa_)+abs(ya_))/100.){
#ifdef TRIANGULATE_D3RECTANGLE
	subelement.push_back(new D3triangle(xa,ya,xb,yb,0,0,0,0,0,0,inversenorm,this));//triangle at the left bottom corner
	subelement.push_back(new D3triangle(xb,yb,xc,yc,0,0,0,0,0,0,inversenorm,this));//triangle at the right bottom corner
#else
	isBaseElement=true;		
	refineable=true;
	x1=0;y1=0;z1=0;
	x2=xa;y2=ya;z2=0;
	x3=xb;y3=yb;z3=0;
	x4=xc;y4=yc;z4=0;
	init();
#endif

}

D3rectangle::D3rectangle(char *str,double x1_,double y1_,double z1_,double x2_,double y2_,double z2_,double x3_,double y3_,double z3_,double x4_,double y4_,double z4_,
			bool inversenorm,PD3element parent):x1(x1_),y1(y1_),z1(z1_),x2(x2_),y2(y2_),z2(z2_),x3(x3_),y3(y3_),z3(z3_),x4(x4_),y4(y4_),z4(z4_),
					D3element(x1_,y1_,z1_,0,0,0,inversenorm,false,parent,0){
#ifdef TRIANGULATE_D3RECTANGLE
	x=0;y=0;z=0;
	subelement.push_back(new D3triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm,this));
	subelement.push_back(new D3triangle(x1,y1,z1,x3,y3,z3,x4,y4,z4,inversenorm,this));
	
#else
	isBaseElement=true;		
	refineable=true;
	double nx,ny,nz,nnorm,nxynorm,za,zb,zc,xa2,ya2,za2,xb2,yb2,zb2,xc2,yc2,zc2;
	xa2=x2-x1;xb2=x3-x1;xc2=x4-x1;
	ya2=y2-y1;yb2=y3-y1;yc2=y4-y1;
	za2=z2-z1;zb2=z3-z1;zc2=z4-z1;
	epsilon=(abs(xa2)+abs(ya2)+abs(za2))/100.;
	//n= a x b
	nx=-za2*yb2 + ya2*zb2;
	ny=-xa2*zb2 + za2*xb2;
	nz=-ya2*xb2 + xa2*yb2;
	
	nnorm=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
	nx/=nnorm;
	ny/=nnorm;
	nz/=nnorm;
	nxynorm=sqrt(sqr(nx)+sqr(ny));
	theta=acos(nz);
	phi=0;
	if (nxynorm!=0){
		nx/=nxynorm;
		ny/=nxynorm;
		psi=acos(-ny);
		if (nx<0) psi=-psi;
	}
	else psi=0;
	psi=-psi;
	theta=-theta;

	rotateinv(phi,theta,psi,xa2,ya2,za2,xa,ya,za);
	rotateinv(phi,theta,psi,xb2,yb2,zb2,xb,yb,zb);
	rotateinv(phi,theta,psi,xc2,yc2,zc2,xc,yc,zc);
	
	init(true);
#endif
}

void D3rectangle::SetNormTowards(double x0,double y0,double z0,bool towards){
	double nx,ny,nz,xa2,ya2,za2,xb2,yb2,zb2,x_,y_,z_;
	xa2=x2-x1;xb2=x3-x1;
	ya2=y2-y1;yb2=y3-y1;
	za2=z2-z1;zb2=z3-z1;
	
	//n= a x b
	nx=-za2*yb2 + ya2*zb2;
	ny=-xa2*zb2 + za2*xb2;
	nz=-ya2*xb2 + xa2*yb2;

	GetReferencePoint(x_,y_,z_);
	double nDOTxMINUSx0=nx*(x_-x0)+ny*(y_-y0)+nz*(z_-z0);
	
	inversenorm= ((nDOTxMINUSx0<0) ^ towards);
	
}
bool D3rectangle::IntersectWithRay(double x0,double y0,double z0,double dirx,double diry,double dirz,double &mul){
	// a=xyz2-xyz1, b=xyz3-xyz1;
	// c=xyz0-xyz1;
	// n=a x b
	// Eqation to solve: c+mul*dir == amul*a+bmul*b
	// |dir . n| * a = sign(dir . n) * dir . (c x b)
	// |dir . n| * b = sign(dir . n) * dir . (a x c)
	// |dir . n| * mul = -sign(dir . n) * (c . n)
	
	double ax,ay,az,bx,by,bz,cx,cy,cz;
	ax=x2-x1;
	ay=y2-y1;
	az=z2-z1;
	
	bx=x3-x1;
	by=y3-y1;
	bz=z3-z1;

	cx=x0-x1;
	cy=y0-y1;
	cz=z0-z1;
	if (D3triangle::IntersectWithRay(ax,ay,az,bx,by,bz,cx,cy,cz,dirx,diry,dirz,mul)) return true;
	ax=x3-x1;
	ay=y3-y1;
	az=z3-z1;
	
	bx=x4-x1;
	by=y4-y1;
	bz=z4-z1;

	cx=x0-x1;
	cy=y0-y1;
	cz=z0-z1;

	return D3triangle::IntersectWithRay(ax,ay,az,bx,by,bz,cx,cy,cz,dirx,diry,dirz,mul);
}

void D3rectangle::GetReferencePoint(double &x,double &y, double &z){x=(x1+x2+x3+x4)/4.,y=(y1+y2+y3+y4)/4.;z=(z1+z2+z3+z4)/4.;}
double D3rectangle::GetSelfPotential(){
	double xp,yp,zp,L11,L12,L13,L14;
	GetReferencePoint(xp,yp,zp);
	L11=calc.GetIntL1(xp,yp,zp,x1,y1,z1,x2,y2,z2);
	L12=calc.GetIntL1(xp,yp,zp,x2,y2,z2,x3,y3,z3);
	L13=calc.GetIntL1(xp,yp,zp,x3,y3,z3,x4,y4,z4);
	L14=calc.GetIntL1(xp,yp,zp,x4,y4,z4,x1,y1,z1);
	return (L11+L12+L13+L14);
}
double D3rectangle::GetPotentialAt(double x,double y,double z){
	return	calc.G3Danalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm)+
			calc.G3Danalytic(x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1,inversenorm);
}
double D3rectangle::GetSelfDoubleLayerPotential(){return 0;}
double D3rectangle::GetDoubleLayerPotentialAt(double x,double y,double z){
	return calc.G3DdnAnalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm)+//undocument if nonanalytical sollution is wanted
		   calc.G3DdnAnalytic(x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1,inversenorm);
	double ax,ay,az,bx,by,bz,nx,ny,nz,nasq,nbsq,n,dG3Ddx,dG3Ddy,dG3Ddz,dG3Ddx2,dG3Ddy2,dG3Ddz2;
		ax=x2-x1;ay=y2-y1;az=z2-z1;
		bx=x3-x1;by=y3-y1;bz=z3-z1;
	
		
		
		nx=(-(az*by) + ay*bz);
		ny=(az*bx - ax*bz);
		nz=(-(ay*bx) + ax*by);
		n=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
		
		if(inversenorm) n=-n;

		dG3Ddx=calc.triangleint(calcD3triangle::dG3Ddx,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);
		dG3Ddy=calc.triangleint(calcD3triangle::dG3Ddy,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);
		dG3Ddz=calc.triangleint(calcD3triangle::dG3Ddz,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3);

		dG3Ddx2=calc.triangleint(calcD3triangle::dG3Ddx,x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1);
		dG3Ddy2=calc.triangleint(calcD3triangle::dG3Ddy,x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1);
		dG3Ddz2=calc.triangleint(calcD3triangle::dG3Ddz,x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1);

		return (nx*(dG3Ddx+dG3Ddx2)+ny*(dG3Ddy+dG3Ddy2)+nz*(dG3Ddz+dG3Ddz2))/n;

}

void D3rectangle::GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez){
	double pot2,ex2,ey2,ez2;
	calc.triangleint_pot_exyz(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm,pot,ex,ey,ez);
	calc.triangleint_pot_exyz(x,y,z,x3,y3,z3,x4,y4,z4,x1,y1,z1,inversenorm,pot2,ex2,ey2,ez2);
	pot+=pot2;
	ex+=ex2;
	ey+=ey2;
	ez+=ez2;
}
void D3rectangle::GetRectangle(double *A,double *B,double *C,double *D,double *COL){

	if(inversenorm){
		A[0]=x4;A[1]=y4;A[2]=z4;
		B[0]=x3;B[1]=y3;B[2]=z3;
		C[0]=x2;C[1]=y2;C[2]=z2;
		D[0]=x1;D[1]=y1;D[2]=z1;
		COL[0]=0;COL[1]=0.5;COL[2]=0.7;
	}
	else{
		A[0]=x1;A[1]=y1;A[2]=z1;
		B[0]=x2;B[1]=y2;B[2]=z2;
		C[0]=x3;C[1]=y3;C[2]=z3;
		D[0]=x4;D[1]=y4;D[2]=z4;
		COL[0]=0;COL[2]=0.5;COL[1]=0.7;
	}
}

void D3rectangle::rotate(double phi_,double theta_,double psi_){
	rotate2(phi_,theta_,psi_,x1,y1,z1,x1,y1,z1);
	rotate2(phi_,theta_,psi_,x2,y2,z2,x2,y2,z2);
	rotate2(phi_,theta_,psi_,x3,y3,z3,x3,y3,z3);
	rotate2(phi_,theta_,psi_,x4,y4,z4,x4,y4,z4);
}
void D3rectangle::shift(double xs,double ys,double zs){
	x1+=xs;y1+=ys;z1+=zs;
	x2+=xs;y2+=ys;z2+=zs;
	x3+=xs;y3+=ys;z3+=zs;
	x4+=xs;y4+=ys;z4+=zs;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////D3world::
D3world::D3world(unsigned long num1,unsigned long num2,unsigned long num3,unsigned long num4,unsigned long d,unsigned long m,unsigned long y):D3element(0,0,0,0,0,0,false,false,NULL),refreshchecksum(false){//register in registry
	//num is registration key, d is day, m is month, y is year
	{
#ifdef CHECKLICENSE
#include "license1.cpp"
#else
}
#endif
};

int D3world::cut(int n,int max){
		rangeerror|=(n>max);
		return (n<max)?n:max-1;
}


D3world::D3world(const char *_cachefilename,double tol,int maxit,int numMom,int numLev,double spaceunit,int segmentation):D3element(0,0,0,0,0,0,false,false,NULL),f(NULL),dfdnAll(NULL),el(NULL),electrodeIndexLimit(NULL),max_diri(0),ave_diri(0),
tol(tol),maxit(maxit),numMom(numMom),numLev(numLev),max_neum(0),ave_neum(0),checksum(0),docache(false),spaceunit(spaceunit),segmentation(segmentation),refreshchecksum(false),rangeerror(false),xscale(1),yscale(1),zscale(1){
	cout <<"bemsolver Version 2.00 15.08.2009\nbug report to: kilian.singer@uni-ulm.de\n";
	strcpy(cachefilename,_cachefilename);
#ifdef CHECKLICENSE
#include "license2.cpp"
#endif
};

D3world::~D3world(){
	if (docache){
		delete[] potcache;
		delete[] feldxcache;
		delete[] feldycache;
		delete[] feldzcache;
	}
	delete[] f;
	delete[] dfdnAll;
	delete[] electrodeIndexLimit;
	delete[] el;
	delete[] electrodes;

	
};
void D3world::insert(D3electrode *el){
	el->parent=this;
	subelement.push_back(el);//triangle at the left bottom corner
};

typedef D3electrode *PD3electrode;
void D3world::GetListOfBaseElements(PD3element *el,int *electrodeIndexLimit,int &cnt){
	if (isBaseElement) {
		el[cnt++]=this;
		return;
	}
	PD3elementList::iterator i;
	int electrodecnt;
	electrodes=new PD3electrode[amountOfElectrodes];
	for(i=subelement.begin(),electrodecnt=0;i!=subelement.end();++i,electrodecnt++){
		electrodeIndexLimit[electrodecnt]=cnt;
		(*i)->GetListOfBaseElements(el,cnt);
		electrodes[electrodecnt]=(D3electrode*)(*i);
		((D3electrode*)(*i))->SetCardinalNumber(electrodecnt);
		
	}
	electrodeIndexLimit[amountOfElectrodes]=cnt;

};


/**
 * Default constructor.
 */
D3ImportedElectrodes::D3ImportedElectrodes():importingPolyline(false),ignore_Model_Space_Handle(true) {}

/**
 * Destructor
 */
D3ImportedElectrodes::~D3ImportedElectrodes() { }

bool D3ImportedElectrodes::Import(const char* file,bool ignore3DFace_,bool ignorePolyline_) {
    // Load DXF file into memory:
	ignore3DFace=ignore3DFace_;
	ignorePolyline=ignorePolyline_;
    std::cout << "Reading file " << file << "...\n";
   
    DL_Dxf* dxf = new DL_Dxf();
    if (!dxf->in(file, this)) { // if file open failed
        std::cerr << file << " could not be opened.\n";
        return false;
    }
    delete dxf;
	return true;

}

void D3ImportedElectrodes::ListElectrodes(){
   map<string,D3electrode>::iterator elit; // iterates over all electrodes
   for (elit = electrode.begin(); elit != electrode.end(); ++elit){
	   cout <<"Electrode :"<<(*elit).first<<"\n";
	}
}

D3electrode & D3ImportedElectrodes::FindElectrode(const char* name){
	map<string,D3electrode>::iterator it;
	it=electrode.find(name);
	if(it!=electrode.end()) return (*it).second;
	else{
		cerr <<"ERROR: Cannot find electrode: "<<name<<endl;
		exit(-1);
	}
}

/**
 * Sample implementation of the method which handles layers.
 */
void D3ImportedElectrodes::addLayer(const DL_LayerData& data) {
    printf("LAYER: %s flags: %d\n", data.name.c_str(), data.flags);
    printAttributes();
	/*D3electrode dummy;//this code needs one more contructor/destructor call than the following one
	pair<string,D3electrode> pair(data.name,dummy);
	electrode.insert(pair);
	*/
	electrode[data.name]; //ctor cp cp destr destr
	(electrode[data.name]).name=data.name;
}


/**
 * Sample implementation of the method which handles point entities.
 */
void D3ImportedElectrodes::addPoint(const DL_PointData& data) {
    //printf("POINT    (%6.3f, %6.3f, %6.3f)\n",
    //       data.x, data.y, data.z);
    //printAttributes();
}

/**
 * Sample implementation of the method which handles line entities.
 */
void D3ImportedElectrodes::addLine(const DL_LineData& data) {
    //printf("LINE     (%6.3f, %6.3f, %6.3f) (%6.3f, %6.3f, %6.3f)\n",
    //       data.x1, data.y1, data.z1, data.x2, data.y2, data.z2);
    //printAttributes();
}

/**
 * Sample implementation of the method which handles arc entities.
 */
void D3ImportedElectrodes::addArc(const DL_ArcData& data) {
    //printf("ARC      (%6.3f, %6.3f, %6.3f) %6.3f, %6.3f, %6.3f\n",
    //       data.cx, data.cy, data.cz,
    //      data.radius, data.angle1, data.angle2);
    //printAttributes();
}

/**
 * Sample implementation of the method which handles circle entities.
 */
void D3ImportedElectrodes::addCircle(const DL_CircleData& data) {
    //printf("CIRCLE   (%6.3f, %6.3f, %6.3f) %6.3f\n",
    //       data.cx, data.cy, data.cz,
    //       data.radius);
    //printAttributes();
}


/**
 * Sample implementation of the method which handles polyline entities.
 */
void D3ImportedElectrodes::addPolyline(const DL_PolylineData& data) {
    //printf("POLYLINE \n");
    //printf("flags: %d\n", (int)data.flags);
    //printAttributes();
	if(!ignore_Model_Space_Handle)
		if(data.owner!=Model_Space_Handle) return;

	if(ignorePolyline) return;
	importingPolyline=true;
	vertexcnt=0;
}
 
/**
 * Called when a SEQEND occurs (when a POLYLINE or ATTRIB is done)
 */
void D3ImportedElectrodes::endSequence(){
	if(ignorePolyline) return;
	if(vertexcnt==3) electrode[attributes.getLayer()].insert(new D3triangle(x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2],true));
	else if(vertexcnt==4) electrode[attributes.getLayer()].insert(new D3rectangle("",x[0],y[0],z[0],x[1],y[1],z[1],x[2],y[2],z[2],x[3],y[3],z[3],true));
	importingPolyline=false;
}
/**
 * Sample implementation of the method which handles vertices.
 */
void D3ImportedElectrodes::addVertex(const DL_VertexData& data) {
	if(ignorePolyline) return;
	if(importingPolyline){
		if (vertexcnt<=3){ 
			x[vertexcnt]=data.x;y[vertexcnt]=data.y;z[vertexcnt]=data.z;
			vertexcnt++;
		}
	}
    //printf("VERTEX   (%6.3f, %6.3f, %6.3f) %6.3f\n",
    //       data.x, data.y, data.z,
    //       data.bulge);
    //printAttributes();
}

void D3ImportedElectrodes::add3DFace(const DL_3DFaceData &data){
	if(!ignore_Model_Space_Handle)
		if(data.owner!=Model_Space_Handle) return;
	if(ignore3DFace) return;
	if (data.cx==data.dx && data.cy==data.dy && data.cz==data.dz) electrode[attributes.getLayer()].insert(new D3triangle(data.ax,data.ay,data.az,data.bx,data.by,data.bz,data.cx,data.cy,data.cz,true));//dreieck
	else electrode[attributes.getLayer()].insert(new D3rectangle("",data.ax,data.ay,data.az,data.bx,data.by,data.bz,data.cx,data.cy,data.cz,data.dx,data.dy,data.dz,true));//viereck
	importingPolyline=false;


	/* printf("3DFace   (%6.3f, %6.3f, %6.3f,   %6.3f, %6.3f, %6.3f   %6.3f, %6.3f, %6.3f   %6.3f, %6.3f, %6.3f)\n",
           data.ax, data.ay, data.az,
		   data.bx, data.by, data.bz,
		   data.cx, data.cy, data.cz,
		   data.dx, data.dy, data.dz
           );
    printAttributes();*/
}

void D3ImportedElectrodes::addBlockRecord(const DL_BlockRecordData &data){
	if(data.name=="*Model_Space") {
		Model_Space_Handle=data.handle;
		ignore_Model_Space_Handle=false;
	}
	printf("BlockRecord %d %s\n",data.handle,data.name.c_str());
}

void D3ImportedElectrodes::printAttributes() {
    printf("  Attributes: Layer: %s, ", attributes.getLayer().c_str());
    printf(" Color: ");
    if (attributes.getColor()==256)	{
        printf("BYLAYER");
    } else if (attributes.getColor()==0) {
        printf("BYBLOCK");
    } else {
        printf("%d", attributes.getColor());
    }
    printf(" Width: ");
    if (attributes.getWidth()==-1) {
        printf("BYLAYER");
    } else if (attributes.getWidth()==-2) {
        printf("BYBLOCK");
    } else if (attributes.getWidth()==-3) {
        printf("DEFAULT");
    } else {
        printf("%d", attributes.getWidth());
    }
    printf(" Type: %s\n", attributes.getLineType().c_str());
}

void D3world::Dcentroid(int shape,double * pc,double* xcout)

{
  double corner[4][3], X[3], Y[3], Z[3], vertex1[3], vertex3[3];
  double sum, delta, dl, x1, y1, x2, x3, y3, xc, yc;
  int i, j;
  //double normalize();
  /* Load the corners. */
  for(i=0; i<4; i++) { 
      for(j=0; j<3; j++) { 
	  corner[i][j] = *(pc++);
      }
  }

  /* Use vertex 0 as the origin and get diags and lengths. */
  for(sum=0, i=0; i<3; i++) {
    X[i] = delta = corner[2][i] - corner[0][i];
    sum += delta * delta;
    vertex1[i] = corner[1][i] - corner[0][i];
    if(shape == QUADRILAT) {
      vertex3[i] = corner[3][i] - corner[0][i];
      Y[i] = corner[1][i] - corner[3][i];
    }
    else if(shape == TRIANGLE) {
      vertex3[i] = corner[2][i] - corner[0][i];
      Y[i] = corner[1][i] - corner[0][i];
    }
    else {
      printf("Dcentroid FE: Shape indicator is neither triangle nor quadrilateral");
      exit(0);
    }
  }
  x2 = sqrt(sum);

  /* Z-axis is normal to two diags. */
  Cross_Product(X, Y, Z);
  normalize(X);
  normalize(Z);

  /* Real Y-axis is normal to X and Z. */
  Cross_Product(Z, X, Y);

  /* Project into the panel axes. */
  y1 = Dot_Product(vertex1, Y);
  y3 = Dot_Product(vertex3, Y);
  x1 = Dot_Product(vertex1, X);
  x3 = Dot_Product(vertex3, X);

  yc = ONE3 * (y1 + y3);
  xc = ONE3 * (x2 + ((x1 * y1 - x3 * y3)/(y1 - y3)));

  *(xcout+0) = corner[0][XI] + xc * X[XI] + yc * Y[XI];
  *(xcout+1) = corner[0][YI] + xc * X[YI] + yc * Y[YI];
  *(xcout+2) = corner[0][ZI] + xc * X[ZI] + yc * Y[ZI];
}

void D3world::save(char *fname){
	ofstream ofs(fname,ios::out|ios::binary);
	ofs<<n<<"\n";
	ofs<<amountOfElectrodes<<"\n";
	ofs<< tol<<"\n";
	ofs<< maxit<<"\n";
	ofs<< numMom<<"\n";
	ofs<< numLev<<"\n";
	ofs<<checksum<<"\n";
	ofs.write((char *)dfdnAll,n*amountOfElectrodes*sizeof(double));
	ofs.close();
}

void D3world::savecalc(char *fname){
	ofstream ofs(fname,ios::out|ios::binary);
	ofs.precision(20);
	ofs<<n<<"\n";
	ofs<<amountOfElectrodes<<"\n";
	ofs<< tol<<"\n";
	ofs<< maxit<<"\n";
	ofs<< numMom<<"\n";
	ofs<< numLev<<"\n";
	ofs<<checksum<<"\n";
	ofs<<xmin<<"\n";
	ofs<<xmax<<"\n";
	ofs<<nx<<"\n";
	ofs<<ymin<<"\n";
	ofs<<ymax<<"\n";
	ofs<<ny<<"\n";
	ofs<<zmin<<"\n";
	ofs<<zmax<<"\n";
	ofs<<nz<<"\n";

	ofs.write((char *)potcache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	int ii;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldxcache[ii]*=spaceunit;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldycache[ii]*=spaceunit;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldzcache[ii]*=spaceunit;
	ofs.write((char *)feldxcache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	ofs.write((char *)feldycache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	ofs.write((char *)feldzcache,nx*ny*nz*amountOfElectrodes*sizeof(double));

	ofs.close();
}

//UP020506///////////////////////////////////////////////////////////////////////////////

void D3world::AssignColors(){
	
	double *SurfaceCharge;
	double MinSurfaceCharge,MaxSurfaceCharge;
	int elcnt,row,cnt;
	int elcnt2;
	double V;
	PD3elementList::iterator i;
	
	SurfaceCharge=new double[n];
	for(row=0;row<n;row++){

		SurfaceCharge[row]=0;
		
		for(elcnt=0;elcnt<amountOfElectrodes;elcnt++){
			for(i=subelement.begin(),elcnt2=0;i!=subelement.end();++i,++elcnt2){ 
				if(elcnt==elcnt2) V=((D3electrode*)(*i))->GetVoltage();
			}
			SurfaceCharge[row]+=V*dfdnAll[row+elcnt*n];	
		}
	}

	MinSurfaceCharge=1e10;
	MaxSurfaceCharge=-1e10;
	for(row=0;row<n;row++){
		//printf("%i %lf %lf %lf\n",row,SurfaceCharge[row],MaxSurfaceCharge,MinSurfaceCharge);
		if(SurfaceCharge[row]>MaxSurfaceCharge) MaxSurfaceCharge=SurfaceCharge[row];
		if(SurfaceCharge[row]<MinSurfaceCharge) MinSurfaceCharge=SurfaceCharge[row];
	}

	for(row=0;row<n;row++){
		if(SurfaceCharge[row]>0)((D3triangle*)el[row])->Color=100+100.0*SurfaceCharge[row]/MaxSurfaceCharge;
		else ((D3triangle*)el[row])->Color=100-100.0*SurfaceCharge[row]/MinSurfaceCharge;
		
		//if(row%5==0) printf("%i %i %lf\n",row,((D3triangle*)el[row])->Color,SurfaceCharge[row]);
	}
	system("PAUSE");

	double r,g,b;
	for(cnt=1;cnt<=200;cnt++) {
		if(cnt<100) {
			r=0;
			g=(double)cnt/100.0;
			b=1.0-(double)cnt/100.0;
		} else {
			r=(double)(cnt-100)/100.0;
			g=1.0-(double)(cnt-100)/100.0;
			b=0;
		}		
	}
}
bool D3world::loadcalc(char *fname){
	checksum=	update_adler32double(BEMREVISION,x,n*VERTS*DIMEN);
	checksum=	update_adler32(checksum,(unsigned char*)shape,n*sizeof(int));
	unsigned long checksum2;
	double tol2,xmin2,xmax2,ymin2,ymax2,zmin2,zmax2;
	int n2,amountOfElectrodes2,maxit2,numMom2,numLev2,nx2,ny2,nz2;
	ifstream ifs(fname,ios::in|ios::binary);
	if(ifs.fail()) return false;
	
	ifs>>n2;
	ifs>>amountOfElectrodes2;
	ifs>>tol2;
	ifs>>maxit2;
	ifs>>numMom2;
	ifs>>numLev2;
	ifs>>checksum2;



	ifs>>xmin2;
	ifs>>xmax2;
	ifs>>nx2;
	ifs>>ymin2;
	ifs>>ymax2;
	ifs>>ny2;
	ifs>>zmin2;
	ifs>>zmax2;
	ifs>>nz2;

	




	if(n2!=n || amountOfElectrodes!=amountOfElectrodes2 ||
		tol2!=tol || maxit2!=maxit||numMom2!=numMom||numLev2!=numLev||
		checksum2!=checksum||
		xmin2!=xmin||xmax2!=xmax||nx2!=nx||
		ymin2!=ymin||ymax2!=ymax||ny2!=ny||
		zmin2!=zmin||zmax2!=zmax||nz2!=nz) {
		cout <<"Checksum or settings for calc different! checksum:"<<checksum<<" <> "<<checksum2<<"\n";
		if(refreshchecksum){
			cout <<"refreshing checksum!"<<endl;
		}
		else{
			ifs.close();
			return false;
		}
	}
	cout <<"Loading calc cache!\n";

	
	char slashn;
	ifs.read(&slashn,1);

	ifs.read((char *)potcache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	ifs.read((char *)feldxcache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	ifs.read((char *)feldycache,nx*ny*nz*amountOfElectrodes*sizeof(double));
	ifs.read((char *)feldzcache,nx*ny*nz*amountOfElectrodes*sizeof(double));
//	int ii;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldxcache[ii]/=spaceunit;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldycache[ii]/=spaceunit;
	//for(ii=0;ii<nx*ny*nz*amountOfElectrodes;ii++) feldzcache[ii]/=spaceunit;
	ifs.close();
	if(refreshchecksum)	savecalc(fname);
	return true;
}

void D3world::exportGeometry(const char *fname){
	ofstream ofs(fname,ios::out|ios::binary);
	ofs <<n<<"\n";
	ofs.write((char *)x,n*VERTS*DIMEN*sizeof(double));
	ofs.write((char *)shape,n*sizeof(int));
	ofs.close();
}
bool D3world::load(char *fname){

	checksum=	update_adler32double(BEMREVISION,x,n*VERTS*DIMEN);
	checksum=	update_adler32(checksum,(unsigned char*)shape,n*sizeof(int));


	unsigned long checksum2;
	double tol2;
	int n2,amountOfElectrodes2,maxit2,numMom2,numLev2;
	ifstream ifs(fname,ios::in|ios::binary);
	if(ifs.fail()) return false;
	
	ifs>>n2;
	ifs>>amountOfElectrodes2;
	ifs>>tol2;
	ifs>>maxit2;
	ifs>>numMom2;
	ifs>>numLev2;
	ifs>>checksum2;

	if(n2!=n || amountOfElectrodes!=amountOfElectrodes2 ||
		tol2!=tol || maxit2!=maxit||numMom2!=numMom||numLev2!=numLev||
		checksum2!=checksum) {
		cout <<"Checksum or settings for solve different! checksum:"<<checksum<<" <> "<<checksum2<<"\n";
		
		if(refreshchecksum){
			cout <<"refreshing checksum!"<<endl;
		}
		else{
			return false;
			ifs.close();
		}
	}
	cout <<"Loading solve cache!\n";

	
	char slashn;
	ifs.read(&slashn,1);
	ifs.read((char *)dfdnAll,n*amountOfElectrodes*sizeof(double));
	ifs.close();
	if(refreshchecksum) save(fname);
	
	return true;
}

void D3world::RefreshChecksum(){
	refreshchecksum=true;
}

unsigned long D3world::update_adler32(unsigned long old, unsigned char *buf, unsigned long len)
{
	#define BASE 65521                     /* largest prime smaller than 65536 */
	unsigned long s1 = old & 0xffff;
	unsigned long s2 = (old >> 16) & 0xffff;
	int n;
	for (n = 0; n < len; n++)
	{
		s1 = (s1 + buf[n]) % BASE;
		s2 = (s2 + s1) % BASE;
	}
	return( (s2 << 16) + s1 );
} // update_adler32() 

unsigned long D3world::update_adler32double(unsigned long old, double *buf, unsigned long len,int ignorebytes)
{
	#define BASE 65521                     /* largest prime smaller than 65536 */
	
	//strstream s;
	//s.precision(13);
	//s <<scientific;

	len*=sizeof(double);

	#define BASE 65521                     /* largest prime smaller than 65536 */
	unsigned long s1 = old & 0xffff;
	unsigned long s2 = (old >> 16) & 0xffff;
	int n;
	for (n = 0; n < len; n++)
	{
		if(n%sizeof(double)<=ignorebytes) continue;
		s1 = (s1 + ((char*)buf)[n]) % BASE;
		s2 = (s2 + s1) % BASE;
	}
	return( (s2 << 16) + s1 );
} // update_adler32() 


void D3world::solve(){
	//set cardinalnumbers;
	//main3(0,NULL);
	//return;
 	n=GetAmountOfSubelements();
	
	amountOfElectrodes=subelement.size();

	shape=new int[n];
	x=new double[n*VERTS*DIMEN];
	f=new double[n];
	dfdnAll=new double[n*amountOfElectrodes];
	double *dfdn;
	type=new int[n];
	xcoll = new double[n*DIMEN];
	xnrm = new double[n*DIMEN];
	dtype = new int[n];
	lhsvect = new double[n];
	rhsvect = new double[n];
	rhstype = new int[n]; 
	lhstype = new int[n];
	rhsindex = new int[n*VERTS];
	lhsindex = new int[n*VERTS];


	//delete[] electrodeIndexLimit;
	//delete[] el;








dfdn=dfdnAll;
i=0;
 int xx;
 int yy;

 
   i=0;
/*
  
   
   for(xx=0;xx<=19;xx++){
		for(yy=0;yy<=18;yy++){
			shape[i]=TRIANGLE;//QUADRILAT;//TRIANGLE
			x[i*VERTS*DIMEN]=xx;
			x[i*VERTS*DIMEN+1]=yy;
			x[i*VERTS*DIMEN+2]=0;
			x[i*VERTS*DIMEN+3]=xx+1;
			x[i*VERTS*DIMEN+4]=yy;
			x[i*VERTS*DIMEN+5]=0;
			x[i*VERTS*DIMEN+6]=xx+1;
			x[i*VERTS*DIMEN+7]=yy+1;
			x[i*VERTS*DIMEN+8]=0;
			//x[i*VERTS*DIMEN+9]=xx;
			//x[i*VERTS*DIMEN+10]=yy+1;
			//x[i*VERTS*DIMEN+11]=0;
			f[i]=(xx>=9 && xx<=11 && yy>=9 && yy<=11)?1:0;
			type[i]=DIRICHLET;
			i++;
		}
	}

	size=i;
*/ 
	electrodeIndexLimit=new int[amountOfElectrodes+1];
	el=new PD3element[n];
	int cnt=0;
	GetListOfBaseElements(el,electrodeIndexLimit,cnt);
	
	int row,col;
	int elcnt;
	double a[3],b[3],c[3],d[3],color[3];
	for(col=0;col<n;col++){//geht alle Flächenelemente durch
		type[col]=DIRICHLET;
		f[col]=1;
		dfdn[col]=0;
		if(dynamic_cast<D3triangle *>(el[col])) {
			el[col]->GetTriangle(a,b,c,color);
			shape[col]=TRIANGLE;
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];
			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];
			x[col*VERTS*DIMEN+9]=0;
			x[col*VERTS*DIMEN+10]=0;
			x[col*VERTS*DIMEN+11]=0;
//			cout <<"T ";
//			int iii;
//			for(iii=0;iii<=8;iii++)cout <<x[col*VERTS*DIMEN+iii]<<" ";
//			cout <<"1.0 0.0 0\n";*/
			
		}
		else if(dynamic_cast<D3rectangle*>(el[col])) {
			el[col]->GetRectangle(a,b,c,d,color);
			shape[col]=QUADRILAT;//TRIANGLE
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];
			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];	
			x[col*VERTS*DIMEN+9]=d[0];
			x[col*VERTS*DIMEN+10]=d[1];
			x[col*VERTS*DIMEN+11]=d[2];
//			int iii;
//			cout <<"Q ";
//			for(iii=0;iii<=11;iii++)cout <<x[col*VERTS*DIMEN+iii]<<" ";
//			cout <<"1.0 0.0 0\n";*/
		}
		else{
			exit(0);
		}
	}
	int elcnt2;

	if(load(cachefilename)) {
		return;
	}


	


	for(elcnt=0;elcnt<amountOfElectrodes;elcnt++){
		cout <<"*******************************************************\n";
		cout <<"Solving for electrode: "<<elcnt<<"/"<<amountOfElectrodes-1<<"\n";
		cout <<"*******************************************************\n";
	
		PD3elementList::iterator i;
		
		//for(i=subelement.begin(),elcnt2=0;i!=subelement.end();++i,++elcnt2){ 
		//	if(elcnt==elcnt2) ((D3electrode*)(*i))->SetVoltage(1);
		//	else ((D3electrode*)(*i))->SetVoltage(0);
		//}
		
		for(elcnt2=0;elcnt2<amountOfElectrodes;elcnt2++)
			for(col=electrodeIndexLimit[elcnt2];col<electrodeIndexLimit[elcnt2+1];col++) f[col]=(elcnt2==elcnt)?1:0;
		dfdn=dfdnAll+n*elcnt;
		for(col=0;col<n;col++) dfdn[col]=0;

		/* This is a Single Layer formulation. */
		fljob = 2;
		/* Set up for the fastlap call and save the exact solution for comparison with the 
			computed solution.  Note that recovery of the correct signs for Green's Thm. are 
			obtained by kidding fastlap about the signs on the lhs and rhs vectors. */
		int cnt;
		for(cnt=0; cnt<n; cnt++) {
			xnrm[col*DIMEN]=1;
			xnrm[col*DIMEN+1]=1;
			xnrm[col*DIMEN+2]=1;
			if(type[cnt] == DIRICHLET) {
				rhstype[cnt] = CONSTANT_DIPOLE;
				lhstype[cnt] = CONSTANT_SOURCE;
				rhsvect[cnt] = f[cnt];
				rhsindex[cnt*VERTS] = cnt;
				lhsindex[cnt*VERTS] = cnt;
				
				Dcentroid(shape[cnt], &x[cnt*VERTS*DIMEN], &xcoll[cnt*DIMEN]);
				dtype[cnt]=0;
				// fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",cnt, xcoll[cnt*DIMEN],xcoll[cnt*DIMEN+1],xcoll[cnt*DIMEN+2]); 

			}
			else if(type[cnt] == NEUMANN) {
				rhstype[cnt] = CONSTANT_SOURCE;
				lhstype[cnt] = CONSTANT_DIPOLE;
				rhsvect[cnt] = dfdn[cnt];
				rhsindex[cnt*VERTS] = cnt;
				lhsindex[cnt*VERTS] = cnt; 
				dtype[cnt]=0;
				Dcentroid(shape[cnt], &x[cnt*VERTS*DIMEN], &xcoll[cnt*DIMEN]);
				//fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",cnt, xcoll[cnt*DIMEN],xcoll[cnt*DIMEN+1],xcoll[cnt*DIMEN+2]); 
			}
			else {
				//printf("driverc FE: You're missing a boundary condition type");
				exit(0);
			}
		}
		double tol2=tol;
		numit = fastlap(&n,&n,&n,x,shape,dtype,lhstype,rhstype,lhsindex,rhsindex,dfdn,rhsvect,xcoll,xnrm,&numLev,&numMom,&maxit,&tol2,&fljob);				
	}
	
	save(cachefilename);

	

};



void D3world::calc(int xxxnum,double *xxx,double *result){
	/*electrodeIndexLimit=new int[amountOfElectrodes+1];
	el=new PD3element[n];
	int cnt=0;
	GetListOfBaseElements(el,electrodeIndexLimit,cnt);
	*/
	int row,col;
	int elcnt;
	double a[3],b[3],c[3],d[3],color[3];
	 double *dfdn=dfdnAll;
	for(col=0;col<n;col++){//geht alle Flächenelemente durch
		type[col]=NEUMANN;
		f[col]=0;
		if(dynamic_cast<D3triangle *>(el[col])) {
			el[col]->GetTriangle(a,b,c,color);
			shape[col]=TRIANGLE;
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];
			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];
		}
		else if(dynamic_cast<D3rectangle*>(el[col])) {
			el[col]->GetRectangle(a,b,c,d,color);
			shape[col]=QUADRILAT;//TRIANGLE
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];
			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];	
			x[col*VERTS*DIMEN+9]=d[0];
			x[col*VERTS*DIMEN+10]=d[1];
			x[col*VERTS*DIMEN+11]=d[2];
		}
		else{
			int hh=10;
			
		}
	}
	


	 /* This is a Fieldcalculation */
  fljob = 0;
  /* Set up for the fastlap call and save the exact solution for comparison with the 
     computed solution.  Note that recovery of the correct signs for Green's Thm. are 
     obtained by kidding fastlap about the signs on the lhs and rhs vectors. */
	
	

	for(i=0; i<n; i++) {
    
      rhstype[i] = CONSTANT_SOURCE;
	  lhstype[i] = CONSTANT_SOURCE;
     
      rhsvect[i] =0; 
	  PD3elementList::iterator it;
	  for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
		  rhsvect[i]+=dfdnAll[elcnt*n+i]* ((D3electrode*)(*it))->GetVoltage();
	  }

      rhsindex[i*VERTS] = i;
      lhsindex[i*VERTS] = i; 
	 // Dcentroid(shape[i], &x[i*VERTS*DIMEN], &xcoll[i*DIMEN]);
      /* fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",i, xcoll[i*DIMEN],xcoll[i*DIMEN+1],xcoll[i*DIMEN+2]); */
    
  }
  nlhs=xxxnum;
	
  int *dtype2=new int[xxxnum];
  for(i=0;i<xxxnum;i++) dtype2[i]=0;
 

  double tol2=tol;
   numit = fastlap(&nlhs,&n,&n,x,shape,dtype2,lhstype,rhstype,lhsindex,rhsindex,result,rhsvect,xxx,xnrm,&numLev,&numMom,&maxit,&tol2,&fljob);


}

bool D3world::IsEqualSurfaceElement(int amountOfVertices,double eps,double *x,double *X){
#define x1 x[0]
#define y1 x[1]
#define z1 x[2]
#define x2 x[3]
#define y2 x[4]
#define z2 x[5]
#define x3 x[6]
#define y3 x[7]
#define z3 x[8]
#define x4 x[9]
#define y4 x[10]
#define z4 x[11]

#define X1 X[0]
#define Y1 X[1]
#define Z1 X[2]
#define X2 X[3]
#define Y2 X[4]
#define Z2 X[5]
#define X3 X[6]
#define Y3 X[7]
#define Z3 X[8]
#define X4 X[9]
#define Y4 X[10]
#define Z4 X[11]
	if(amountOfVertices==3){
		if(         abs(x1-X1)<eps && abs(y1-Y1)<eps && abs(z1-Z1)<eps){
			if(     abs(x2-X2)<eps && abs(y2-Y2)<eps && abs(z2-Z2)<eps){
				if( abs(x3-X3)<eps && abs(y3-Y3)<eps && abs(z3-Z3)<eps) 
					return true;
			}
			else if(abs(x2-X3)<eps && abs(y2-Y3)<eps && abs(z2-Z3)<eps){
				if( abs(x3-X2)<eps && abs(y3-Y2)<eps && abs(z3-Z2)<eps) 
					return true;
			}
		}
		
		if(         abs(x1-X2)<eps && abs(y1-Y2)<eps && abs(z1-Z2)<eps){
			if(     abs(x2-X3)<eps && abs(y2-Y3)<eps && abs(z2-Z3)<eps){
				if( abs(x3-X1)<eps && abs(y3-Y1)<eps && abs(z3-Z1)<eps) 
					return true;
			}
			else if(abs(x2-X1)<eps && abs(y2-Y1)<eps && abs(z2-Z1)<eps){
				if( abs(x3-X3)<eps && abs(y3-Y3)<eps && abs(z3-Z3)<eps) 
					return true;
			}
		}
		
		
		if(         abs(x1-X3)<eps && abs(y1-Y3)<eps && abs(z1-Z3)<eps){
			if(     abs(x2-X1)<eps && abs(y2-Y1)<eps && abs(z2-Z1)<eps){
				if( abs(x3-X2)<eps && abs(y3-Y2)<eps && abs(z3-Z2)<eps) 
					return true;
			}
			else if(abs(x2-X2)<eps && abs(y2-Y2)<eps && abs(z2-Z2)<eps){
				if( abs(x3-X1)<eps && abs(y3-Y1)<eps && abs(z3-Z1)<eps) 
					return true;
			}
		}
		return false;
	}
	else{
	
		if(            abs(x1-X1)<eps && abs(y1-Y1)<eps && abs(z1-Z1)<eps){
			if(        abs(x2-X2)<eps && abs(y2-Y2)<eps && abs(z2-Z2)<eps){
				if(    abs(x3-X3)<eps && abs(y3-Y3)<eps && abs(z3-Z3)<eps) 
					if(abs(x4-X4)<eps && abs(y4-Y4)<eps && abs(z4-Z4)<eps)
						return true;
			}
			else if(   abs(x2-X4)<eps && abs(y2-Y4)<eps && abs(z2-Z4)<eps){
				if(    abs(x3-X3)<eps && abs(y3-Y3)<eps && abs(z3-Z3)<eps) 
					if(abs(x4-X2)<eps && abs(y4-Y2)<eps && abs(z4-Z2)<eps)
						return true;
			}
		}
		
		if(            abs(x1-X2)<eps && abs(y1-Y2)<eps && abs(z1-Z2)<eps){
			if(        abs(x2-X3)<eps && abs(y2-Y3)<eps && abs(z2-Z3)<eps){
				if(    abs(x3-X4)<eps && abs(y3-Y4)<eps && abs(z3-Z4)<eps) 
					if(abs(x4-X1)<eps && abs(y4-Y1)<eps && abs(z4-Z1)<eps)
						return true;
			}
			else if(   abs(x2-X1)<eps && abs(y2-Y1)<eps && abs(z2-Z1)<eps){
				if(    abs(x3-X4)<eps && abs(y3-Y4)<eps && abs(z3-Z4)<eps) 
					if(abs(x4-X3)<eps && abs(y4-Y3)<eps && abs(z4-Z3)<eps)
						return true;
			}
		}

		if(            abs(x1-X3)<eps && abs(y1-Y3)<eps && abs(z1-Z3)<eps){
			if(        abs(x2-X4)<eps && abs(y2-Y4)<eps && abs(z2-Z4)<eps){
				if(    abs(x3-X1)<eps && abs(y3-Y1)<eps && abs(z3-Z1)<eps) 
					if(abs(x4-X2)<eps && abs(y4-Y2)<eps && abs(z4-Z2)<eps)
						return true;
			}
			else if(   abs(x2-X2)<eps && abs(y2-Y2)<eps && abs(z2-Z1)<eps){
				if(    abs(x3-X1)<eps && abs(y3-Y1)<eps && abs(z3-Z1)<eps) 
					if(abs(x4-X4)<eps && abs(y4-Y4)<eps && abs(z4-Z4)<eps)
						return true;
			}
		}

		if(            abs(x1-X4)<eps && abs(y1-Y4)<eps && abs(z1-Z4)<eps){
			if(        abs(x2-X1)<eps && abs(y2-Y1)<eps && abs(z2-Z1)<eps){
				if(    abs(x3-X2)<eps && abs(y3-Y2)<eps && abs(z3-Z2)<eps) 
					if(abs(x4-X3)<eps && abs(y4-Y3)<eps && abs(z4-Z3)<eps)
						return true;
			}
			else if(   abs(x2-X3)<eps && abs(y2-Y3)<eps && abs(z2-Z3)<eps){
				if(    abs(x3-X2)<eps && abs(y3-Y2)<eps && abs(z3-Z2)<eps) 
					if(abs(x4-X1)<eps && abs(y4-Y1)<eps && abs(z4-Z1)<eps)
						return true;
			}
		}
		return false;

	}
#undef x1
#undef y1
#undef z1
#undef x2
#undef y2
#undef z2
#undef x3
#undef y3
#undef z3
#undef x4
#undef y4
#undef z4

#undef X1
#undef Y1
#undef Z1
#undef X2
#undef Y2
#undef Z2
#undef X3
#undef Y3
#undef Z3
#undef X4
#undef Y4
#undef Z4	
																					
}


enum XYZ {
	X0,Y0,Z0,
	X1,Y1,Z1,
	X2,Y2,Z2,
	X3,Y3,Z3
};


void D3world::exch(double &a,double &b){
	double t;
	t=b;
	b=a;
	a=t;
}


void D3world::RotMirrorSurfaceElement(int axis,int symm,double *xx){//erst rotieren dann spiegeln
	if((symm>>4)!=0){ 
		double phi=double(symm>>4)/double(totalrot)*2.*pi;
		double cosp=cos(phi);
		double sinp=sin(phi);
		double xxx[12];
		switch(axis){//aktive drehung punkt wird gedreht mathematisch positiv
			case 1://x
				//  yz         y            z
				xxx[Y0]=cosp*xx[Y0] - sinp*xx[Z0];
				xxx[Z0]=sinp*xx[Y0] + cosp*xx[Z0];
				xxx[X0]=xx[X0];

				xxx[Y1]=cosp*xx[Y1] - sinp*xx[Z1];
				xxx[Z1]=sinp*xx[Y1] + cosp*xx[Z1];
				xxx[X1]=xx[X1];

				xxx[Y2]=cosp*xx[Y2] - sinp*xx[Z2];
				xxx[Z2]=sinp*xx[Y2] + cosp*xx[Z2];
				xxx[X2]=xx[X2];

				xxx[Y3]=cosp*xx[Y3] - sinp*xx[Z3];
				xxx[Z3]=sinp*xx[Y3] + cosp*xx[Z3];
				xxx[X3]=xx[X3];
				break;
			case 2://y
				// zx          z            x
				xxx[Z0]=cosp*xx[Z0] - sinp*xx[X0];
				xxx[X0]=sinp*xx[Z0] + cosp*xx[X0];
				xxx[Y0]=xx[Y0];

				xxx[Z1]=cosp*xx[Z1] - sinp*xx[X1];
				xxx[X1]=sinp*xx[Z1] + cosp*xx[X1];
				xxx[Y1]=xx[Y1];

				xxx[Z2]=cosp*xx[Z2] - sinp*xx[X2];
				xxx[X2]=sinp*xx[Z2] + cosp*xx[X2];
				xxx[Y2]=xx[Y2];

				xxx[Z3]=cosp*xx[Z3] - sinp*xx[X3];
				xxx[X3]=sinp*xx[Z3] + cosp*xx[X3];
				xxx[Y3]=xx[Y3];
				break;
			case 3://z
				//  xy          x             y
				xxx[X0]=cosp*xx[X0] - sinp*xx[Y0];
				xxx[Y0]=sinp*xx[X0] + cosp*xx[Y0];
				xxx[Z0]=xx[Z0];

				xxx[X1]=cosp*xx[X1] - sinp*xx[Y1];
				xxx[Y1]=sinp*xx[X1] + cosp*xx[Y1];
				xxx[Z1]=xx[Z1];

				xxx[X2]=cosp*xx[X2] - sinp*xx[Y2];
				xxx[Y2]=sinp*xx[X2] + cosp*xx[Y2];
				xxx[Z2]=xx[Z2];

				xxx[X3]=cosp*xx[X3] - sinp*xx[Y3];
				xxx[Y3]=sinp*xx[X3] + cosp*xx[Y3];
				xxx[Z3]=xx[Z3];
				break;
		}
		for(int i=0;i<12;i++)
			xx[i]=xxx[i];
	}

	if(symm&1) {//xsym
		
		xx[X0]=-xx[X0];
		xx[X1]=-xx[X1];
		xx[X2]=-xx[X2];
		xx[X3]=-xx[X3];
	}
	if(symm&2) {//ysym
		
		xx[Y0]=-xx[Y0];
		xx[Y1]=-xx[Y1];
		xx[Y2]=-xx[Y2];
		xx[Y3]=-xx[Y3];
	}
	if(symm&4) {//zsym
		
		xx[Z0]=-xx[Z0];
		xx[Z1]=-xx[Z1];
		xx[Z2]=-xx[Z2];
		xx[Z3]=-xx[Z3];
	}
	if(symm&8) {//spiegelung an der winkelhalb
		switch(axis){
			case 1://x
				exch(xx[Y0],xx[Z0]);
				exch(xx[Y1],xx[Z1]);
				exch(xx[Y2],xx[Z2]);
				exch(xx[Y3],xx[Z3]);
				break;
			case 2://y
				exch(xx[X0],xx[Z0]);
				exch(xx[X1],xx[Z1]);
				exch(xx[X2],xx[Z2]);
				exch(xx[X3],xx[Z3]);
				break;
			case 3://z
				exch(xx[X0],xx[Y0]);
				exch(xx[X1],xx[Y1]);
				exch(xx[X2],xx[Y2]);
				exch(xx[X3],xx[Y3]);
				break;
		}
	}
	
}

void D3world::SymmetrizeCharges(int axis,double epsilon,bool ignoremirror){
	int symmstart;
	int symminc;
	if(ignoremirror){
		symmstart=16;
		symminc=16;
	}
	else{
		symmstart=1;
		symminc=1;
	}
	int elcnt,symm,totalsym;
	totalrot=-1;
	for(elcnt=0;elcnt<amountOfElectrodes;++elcnt){
		if(totalrot<0) {
			totalrot=electrodes[elcnt]->GetTotalrot();
			totalsym=(totalrot<<4)+15+1;
		}
		if(electrodes[elcnt]->GetTotalrot()!=totalrot){
			cerr <<"The total amount of rotational symmetries totalrot must be equal for all electrodes!"<<endl;
			return;
		}
		for(symm=1;symm<totalsym;++symm){
			D3electrode *symel=electrodes[elcnt]->GetSymmetryWith(symm);
			if(symel==electrodes[elcnt]) continue;
			if(symel==NULL) continue;
			bool xsym,ysym,zsym,diagmirror;
			int rotsym;
			xsym=symm&1;ysym=symm&2;zsym=symm&4;diagmirror=symm&8;rotsym=symm>>4;
			int mirrorcnt=0;
			if(xsym && (axis != 1) ) mirrorcnt++;
			if(ysym && (axis != 2) ) mirrorcnt++;
			if(zsym && (axis != 3) ) mirrorcnt++;
			if(diagmirror) mirrorcnt++;
			int backrot;
			if(mirrorcnt&1) backrot=rotsym;//bei ungeraden spiegelungen ist die rückdrehung die gleiche drehung
			else backrot=totalrot-rotsym;


			if(!diagmirror) symel->SetSymmetryWith(xsym,ysym,zsym,diagmirror,electrodes[elcnt],backrot,totalrot);
			else{
				switch(axis){
				case 1://x
					if(ysym&&zsym) symel->SetSymmetryWith(xsym,ysym,zsym,diagmirror,electrodes[elcnt],backrot,totalrot);
					else symel->SetSymmetryWith(xsym,zsym,ysym,diagmirror,electrodes[elcnt],backrot,totalrot);
					break;
				case 2://y
					if(xsym&&zsym) symel->SetSymmetryWith(xsym,ysym,zsym,diagmirror,electrodes[elcnt],backrot,totalrot);
					else symel->SetSymmetryWith(zsym,ysym,xsym,diagmirror,electrodes[elcnt],backrot,totalrot);
					break;
				case 3://z
					if(xsym&&ysym) symel->SetSymmetryWith(xsym,ysym,zsym,diagmirror,electrodes[elcnt],backrot,totalrot);
					else symel->SetSymmetryWith(ysym,xsym,zsym,diagmirror,electrodes[elcnt],backrot,totalrot);
					break;

				}
			}
		}
	}

	D3sortedelement *sortel=new D3sortedelement[n];
	int col;
	double xc,yc,zc;
	for(elcnt=0;elcnt<amountOfElectrodes;++elcnt){
		for(col=electrodeIndexLimit[elcnt];col<electrodeIndexLimit[elcnt+1];col++){
			el[col]->GetCenter(xc,yc,zc);
			switch(axis){
				case 1: 
					sortel[col].axis=xc;
					sortel[col].radius=sqrt(sqr(yc)+sqr(zc));
					break;
				case 2:
					sortel[col].axis=yc;
					sortel[col].radius=sqrt(sqr(xc)+sqr(zc));
					break;
				case 3:
					sortel[col].axis=zc;
					sortel[col].radius=sqrt(sqr(xc)+sqr(yc));
					break;
			}
			sortel[col].index=col;
			sortel[col].electrodeptr=electrodes[elcnt];
		}
	}
	D3sorter sort;
	sort.sort(sortel,n,epsilon);
	int symcnt=0;
	int *foundindex=new int[totalsym+1];
	int *foundsymmetry=new int[totalsym+1];
	int foundanz;
	int i;

	for(col=0;col<n;++col){
		foundanz=0;
		if(!sortel[col].done){
			foundindex[foundanz]=col;
			foundsymmetry[foundanz++]=0;//foundsymmetry[0]=0, so that the first element is the electrode itself, as electrode->symel[0] is this
			double xx[12];
			if(sortel[col].done) continue;
			sortel[col].done=true;
			for(symm=symmstart;symm<totalsym;symm+=symminc){
				for(i=0;i<12;i++) xx[i]=x[sortel[col].index*DIMEN*VERTS+i];
				RotMirrorSurfaceElement(axis,symm,xx);
				for(int col2=col+1;col2<n;++col2){
					
					if(sortel[col2].done) continue;
					if(shape[sortel[col2].index]!=shape[sortel[col].index]) continue;
					if(fabs(sortel[col2].radius-sortel[col].radius)>epsilon) break;
					if(fabs(sortel[col2].axis-sortel[col].axis)>epsilon) break;
					if(IsEqualSurfaceElement(shape[sortel[col].index],epsilon,xx,&x[sortel[col2].index*DIMEN*VERTS])){
						sortel[col2].done=true;//this works now because of second loop
						//old comment: does not work as then Symmetries with respect to other electrodes are not treated
						//so we accept multiple selection
						//ASSERT(foundanz<totalsym);
						foundindex[foundanz]=col2;
						foundsymmetry[foundanz++]=symm;
						symcnt++;
					}
				}
			}
		}
		cout.precision(40);

		if(foundanz<=1) continue;
		symcnt++;//to count also first element
		//now symmetrize charges
		int foundcnt;
		int foundcntstart;
		for(foundcntstart=0;foundcntstart<foundanz-1;++foundcntstart){//we have to iterate over the start of foundcnt in order to treat symmetries with which are ignored if we only start with the first found surface element
			//symmetries of a non first surface element with respect to an other which is adressed by an symmetry with of an electrode
			for(elcnt=0;elcnt<amountOfElectrodes;++elcnt){
				double chargesum=0;
				int chargeanz=0;
				
				for(foundcnt=foundcntstart;foundcnt<foundanz;++foundcnt){
					D3electrode *sym=electrodes[elcnt]->GetSymmetryWith(foundsymmetry[foundcnt]);
					if(sym==NULL) continue;
					chargesum+=dfdnAll[n*sym->GetCardinalNumber()+sortel[foundindex[foundcnt]].index];
					chargeanz++;
				}

				for(foundcnt=foundcntstart;foundcnt<foundanz;++foundcnt){
					D3electrode *sym=electrodes[elcnt]->GetSymmetryWith(foundsymmetry[foundcnt]);
					if(sym==NULL) continue;
					dfdnAll[n*sym->GetCardinalNumber()+sortel[foundindex[foundcnt]].index]=chargesum/double(chargeanz);
				}
				
			}
		}
	}
	cout <<symcnt<<" of "<<n<<" surface elements are symmetric!"<<endl;
	delete[] foundsymmetry;
	delete[] foundindex;
	delete[] sortel;
}

double D3world::calc(double xxx,double y,double z){//PostCalcScaled
	double pot,feldx,feldy,feldz;
	calc(xxx,y,z,pot,feldx,feldy,feldz);
	return pot;	
}

void D3world::calc(double xxx,double y,double z,double &pot,double &feldx,double &feldy,double &feldz){//PostCalcScaled
	xxx/=xscale;
	y/=yscale;
	z/=zscale;
	if(docache){
		int ix,iy,iz;
		double dx,dy,dz,XA,XB,YA,YB,ZA,ZB;
		if(nx>1) dx=(xxx-xmin)/(xmax-xmin)*double(nx-1);
		else dx=0;
		if(ny>1) dy=(y-ymin)/(ymax-ymin)*double(ny-1);
		else dy=0;
		if(nz>1) dz=(z-zmin)/(zmax-zmin)*double(nz-1);
		else dz=0;

		ix=int(dx);
		iy=int(dy);
		iz=int(dz);
		dx-=double(ix);
		dy-=double(iy);
		dz-=double(iz);

		XB=dx;XA=1.-XB;
		YB=dy;YA=1.-YB;
		ZB=dz;ZA=1.-ZB;

		cout.flush();
		PD3elementList::iterator it;
		int elcnt;
		pot=0;feldx=0;feldy=0;feldz=0;
		for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
			pot+=(	 XA*YA*ZA*potcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YA*ZB*potcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZA*potcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZB*potcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZA*potcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZB*potcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZA*potcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZB*potcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt])
						*((D3electrode*)(*it))->GetVoltage();
			feldx+=( XA*YA*ZA*feldxcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YA*ZB*feldxcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZA*feldxcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZB*feldxcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZA*feldxcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZB*feldxcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZA*feldxcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZB*feldxcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt])
						*((D3electrode*)(*it))->GetVoltage();
			feldy+=	(XA*YA*ZA*feldycache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YA*ZB*feldycache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZA*feldycache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZB*feldycache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZA*feldycache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZB*feldycache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZA*feldycache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZB*feldycache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt])
						*((D3electrode*)(*it))->GetVoltage();
			feldz+=	(XA*YA*ZA*feldzcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YA*ZB*feldzcache[cut(ix  ,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZA*feldzcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XA*YB*ZB*feldzcache[cut(ix  ,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZA*feldzcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YA*ZB*feldzcache[cut(ix+1,nx)*nz*ny+cut(iy  ,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZA*feldzcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz  ,nz)+nx*ny*nz*elcnt]
					+XB*YB*ZB*feldzcache[cut(ix+1,nx)*nz*ny+cut(iy+1,ny)*nz+cut(iz+1,nz)+nx*ny*nz*elcnt])
						*((D3electrode*)(*it))->GetVoltage();
				
			
		}

	}
	else{
		double xxxarray[3]={xxx,y,z};
		double potarr[1],feldxarr[1],feldyarr[1],feldzarr[1];
		calc_slow(1,xxxarray,potarr,feldxarr,feldyarr,feldzarr);
		feldx=feldxarr[0];
		feldy=feldyarr[0];
		feldz=feldzarr[0];
		pot= potarr[0];
	}
	feldx/=xscale;
	feldy/=yscale;
	feldz/=zscale;

};

void D3world::calc(double xmin2,double xmax2,int nx2,double ymin2,double ymax2,int ny2,double zmin2,double zmax2,int nz2){
	freeunusedmem();
	xmin=xmin2;
	xmax=xmax2;
	nx=nx2;
	ymin=ymin2;
	ymax=ymax2;
	ny=ny2;
	zmin=zmin2;
	zmax=zmax2;
	nz=nz2;

	docache=true;
	int num=nx*ny*nz;
	potcache=new double[num*amountOfElectrodes];
	feldxcache=new double[num*amountOfElectrodes];
	feldycache=new double[num*amountOfElectrodes];
	feldzcache=new double[num*amountOfElectrodes];
	
	//potcache=new double[segmentation];
	//feldxcache=new double[segmentation];
	//feldycache=new double[segmentation];
	//feldzcache=new double[segmentation];
	
	
	
	char fname[2000];
	strcpy(fname,cachefilename);
	strcat(fname,"calc");
	//if(loadcalc(fname))return;
	if(loadcalc(fname))return;

	double *xxx=new double[num*3];
	
	int ix,iy,iz;
	double x,y,z;

	for(ix=0;ix<nx;ix++){
		if(nx>1) x=xmin+double(ix)*(xmax-xmin)/double(nx-1);
		else x=xmin;
		for(iy=0;iy<ny;iy++){
			if(ny>1) y=ymin+double(iy)*(ymax-ymin)/double(ny-1);
			else y=ymin;
			for(iz=0;iz<nz;iz++){
				if(nz>1) z=zmin+double(iz)*(zmax-zmin)/double(nz-1);
				else z=zmin;
				xxx[(ix*nz*ny+iy*nz+iz)*3]=x;
				xxx[(ix*nz*ny+iy*nz+iz)*3+1]=y;
				xxx[(ix*nz*ny+iy*nz+iz)*3+2]=z;
			}
		}
	}
	
	PD3elementList::iterator it;
	int elcnt,elcnt2;
	double *voltages=new double[amountOfElectrodes];

	for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
	  voltages[elcnt]=((D3electrode*)(*it))->GetVoltage();
	}

	for(elcnt2=0;elcnt2<amountOfElectrodes;elcnt2++){
		for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
			((D3electrode*)(*it))->SetVoltage(elcnt==elcnt2?1:0);
		}
		cout <<"*******************************************************\n";
		cout <<"Calculating electrode: "<<elcnt2<<"/"<<amountOfElectrodes-1<<"\n";
		cout <<"*******************************************************\n";
		currentel=elcnt2;
		calc(num,xxx,potcache+num*elcnt2,feldxcache+num*elcnt2,feldycache+num*elcnt2,feldzcache+num*elcnt2);
		//calc(num,xxx,potcache,feldxcache,feldycache,feldzcache);
		freeunusedmem();
	}
	
	for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
		((D3electrode*)(*it))->SetVoltage(voltages[elcnt]);
	}
	
	savecalc(fname);
	delete[] xxx;
	delete[] voltages;

}
void D3world::calc_slow(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ){///<Calculating potentials and fields.
	
	unsigned long count=xxxnum/segmentation;
	if(xxxnum%segmentation!=0) ++count;
	
	unsigned long pos=0;
	for(int i=0;i<count;i++){
		//cout <<"*******************************************************\n";
		//cout <<"Calculating partition: "<<i<<"/"<<(count-1)<<"\n";
		//cout <<"*******************************************************\n";

		if(i<count-1) calc_slow2(segmentation,xxx+pos*3,Potential+pos,FieldX+pos,FieldY+pos,FieldZ+pos);
		else calc_slow2(xxxnum-pos,xxx+pos*3,Potential+pos,FieldX+pos,FieldY+pos,FieldZ+pos);

		pos+=segmentation;

	}
}

void D3world::calc_slow2(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ){
	n=GetAmountOfSubelements();
	amountOfElectrodes=subelement.size();

	freeunusedmem();
	
	docache=false;
	
	
	
	
	double x,y,z;
	int cnt=0;
	int col;
	double charge,pot,ex,ey,ez;
	int index;


	for(index=0;index<xxxnum;index++){
		Potential[index]=0;
		FieldX[index]=0;
		FieldY[index]=0;
		FieldZ[index]=0;

		int elcnt;
	
		PD3elementList::iterator it;
		for(col=0;col<n;++col){	
			//cout <<col<<" x: "<<xxx[index*3+0]<<" y: "<<xxx[index*3+1]<<" z: "<<xxx[index*3+2]<<endl;
			el[col]->GetPotentialAndFieldAt(xxx[index*3+0],xxx[index*3+1],xxx[index*3+2],pot,ex,ey,ez);
			for(elcnt=0,it=subelement.begin();   elcnt<amountOfElectrodes; ++elcnt,++it){
				double vol=((D3electrode*)(*it))->GetVoltage();
		
				charge=dfdnAll[n*elcnt+col]*vol;
				Potential[index]+=pot*charge;
				FieldX[index]+=ex*charge;
				FieldY[index]+=ey*charge;
				FieldZ[index]+=ez*charge;
			}
		}
		
	}

	long ii;
	for(ii=0;ii<xxxnum;ii++) FieldX[ii]/=spaceunit;
	for(ii=0;ii<xxxnum;ii++) FieldY[ii]/=spaceunit;
	for(ii=0;ii<xxxnum;ii++) FieldZ[ii]/=spaceunit;
}

void D3world::calc_slow(double xmin2,double xmax2,int nx2,double ymin2,double ymax2,int ny2,double zmin2,double zmax2,int nz2){
	n=GetAmountOfSubelements();
	amountOfElectrodes=subelement.size();

	freeunusedmem();
	xmin=xmin2;
	xmax=xmax2;
	nx=nx2;
	ymin=ymin2;
	ymax=ymax2;
	ny=ny2;
	zmin=zmin2;
	zmax=zmax2;
	nz=nz2;

	docache=true;
	int num=nx*ny*nz;
	potcache=new double[num*amountOfElectrodes];
	feldxcache=new double[num*amountOfElectrodes];
	feldycache=new double[num*amountOfElectrodes];
	feldzcache=new double[num*amountOfElectrodes];
	
	//potcache=new double[segmentation];
	//feldxcache=new double[segmentation];
	//feldycache=new double[segmentation];
	//feldzcache=new double[segmentation];
	
	
	
	char fname[2000];
	strcpy(fname,cachefilename);
	strcat(fname,"calc");
	if(loadcalc(fname)) return;
	//if(loadcalc(fname))return;

	double *xxx=new double[num*3];
	
	int ix,iy,iz;
	double x,y,z;
	int cnt=0;
	int col;
	double charge,pot,ex,ey,ez;
	for(ix=0;ix<nx;ix++){
		if(nx>1) x=xmin+double(ix)*(xmax-xmin)/double(nx-1);
		else x=xmin;
		for(iy=0;iy<ny;iy++){
			cout <<cnt++<<"/"<<nx*ny<<endl;
			if(ny>1) y=ymin+double(iy)*(ymax-ymin)/double(ny-1);
			else y=ymin;
			for(iz=0;iz<nz;iz++){
				if(nz>1) z=zmin+double(iz)*(zmax-zmin)/double(nz-1);
				else z=zmin;
				int index=ix*nz*ny+iy*nz+iz;
				int elcnt,index2;
				for(elcnt=0,index2=index;elcnt<amountOfElectrodes;++elcnt,index2+=num){
					potcache[index2]=0;
					feldxcache[index2]=0;
					feldycache[index2]=0;
					feldzcache[index2]=0;
				}
				for(col=0;col<n;++col){	
					el[col]->GetPotentialAndFieldAt(x,y,z,pot,ex,ey,ez);
					for(elcnt=0,index2=index;elcnt<amountOfElectrodes;++elcnt,index2+=num){
						charge=dfdnAll[n*elcnt+col];
						potcache[index2]+=pot*charge;
						feldxcache[index2]+=ex*charge;
						feldycache[index2]+=ey*charge;
						feldzcache[index2]+=ez*charge;
					}
				}
				
			}
		}
	}

	long ii;
	for(ii=0;ii<num*amountOfElectrodes;ii++) feldxcache[ii]/=spaceunit;
	for(ii=0;ii<num*amountOfElectrodes;ii++) feldycache[ii]/=spaceunit;
	for(ii=0;ii<num*amountOfElectrodes;ii++) feldzcache[ii]/=spaceunit;


	savecalc(fname);

}

	
void D3world::calc(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ){
	unsigned long count=xxxnum/segmentation;
	if(xxxnum%segmentation!=0) ++count;
	
	unsigned long pos=0;
	for(int i=0;i<count;i++){
		cout <<"*******************************************************\n";
		cout <<"Calculating partition: "<<i<<"/"<<(count-1)<<"\n";
		cout <<"*******************************************************\n";

		if(i<count-1) calc2(segmentation,xxx+pos*3,Potential+pos,FieldX+pos,FieldY+pos,FieldZ+pos);
		else calc2(xxxnum-pos,xxx+pos*3,Potential+pos,FieldX+pos,FieldY+pos,FieldZ+pos);

		pos+=segmentation;

	}
}

void D3world::calc2(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ){
	double *xxx4=new double[xxxnum*4*DIMEN];
	double *potfieldxyz=new double[xxxnum*4];
	memcpy(xxx4,xxx,xxxnum*DIMEN*sizeof(double));
	memcpy(xxx4+xxxnum*DIMEN,xxx,xxxnum*DIMEN*sizeof(double));
	memcpy(xxx4+2*xxxnum*DIMEN,xxx,xxxnum*DIMEN*sizeof(double));
	memcpy(xxx4+3*xxxnum*DIMEN,xxx,xxxnum*DIMEN*sizeof(double));
	//electrodeIndexLimit=new int[amountOfElectrodes+1];
	//el=new PD3element[n];
	//int cnt=0;
	//GetListOfBaseElements(el,electrodeIndexLimit,cnt);
	int row,col;
	int elcnt;
	double a[3],b[3],c[3],d[3],color[3];
	double *dfdn=dfdnAll;
	for(col=0;col<n;col++){//geht alle Flächenelemente durch
		type[col]=NEUMANN;
		f[col]=0;
		if(dynamic_cast<D3triangle *>(el[col])) {
			el[col]->GetTriangle(a,b,c,color);
			shape[col]=TRIANGLE;
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];

			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];
		}
		else if(dynamic_cast<D3rectangle*>(el[col])) {
			el[col]->GetRectangle(a,b,c,d,color);
			shape[col]=QUADRILAT;//TRIANGLE
			x[col*VERTS*DIMEN]=a[0];
			x[col*VERTS*DIMEN+1]=a[1];
			x[col*VERTS*DIMEN+2]=a[2];
			x[col*VERTS*DIMEN+3]=b[0];
			x[col*VERTS*DIMEN+4]=b[1];
			x[col*VERTS*DIMEN+5]=b[2];
			x[col*VERTS*DIMEN+6]=c[0];
			x[col*VERTS*DIMEN+7]=c[1];
			x[col*VERTS*DIMEN+8]=c[2];	
			x[col*VERTS*DIMEN+9]=d[0];
			x[col*VERTS*DIMEN+10]=d[1];
			x[col*VERTS*DIMEN+11]=d[2];
		}
		else{
			int hh=10;
			
		}
	}
	


	 /* This is a Fieldcalculation */
	fljob = 0;
	/* Set up for the fastlap call and save the exact solution for comparison with the 
     computed solution.  Note that recovery of the correct signs for Green's Thm. are 
     obtained by kidding fastlap about the signs on the lhs and rhs vectors. */
	
	

	for(i=0; i<n; i++) {
    
      rhstype[i] = CONSTANT_SOURCE;
	  lhstype[i] = CONSTANT_SOURCE;
     
      rhsvect[i] =0; 
	  PD3elementList::iterator it;
	  for(elcnt=0,it=subelement.begin();elcnt<amountOfElectrodes;elcnt++,++it){
		  rhsvect[i]+=dfdnAll[elcnt*n+i]* ((D3electrode*)(*it))->GetVoltage();
	  }

      rhsindex[i*VERTS] = i;
      lhsindex[i*VERTS] = i; 
	 // Dcentroid(shape[i], &x[i*VERTS*DIMEN], &xcoll[i*DIMEN]);
      /* fprintf(stdout, "Panel:%d    Centroid:%.8g %.8g %.8g\n",i, xcoll[i*DIMEN],xcoll[i*DIMEN+1],xcoll[i*DIMEN+2]); */
    
	}
	nlhs=xxxnum*4;
	
	int *dtype2=new int[xxxnum*4];
	double *xnrm2=new double[xxxnum*4*DIMEN];
	for(i=0;i<xxxnum;i++){
		dtype2[i]=0;
		xnrm2[i*DIMEN]=1;
		xnrm2[i*DIMEN+1]=1;
		xnrm2[i*DIMEN+2]=1;
	}
	for(i=xxxnum;i<xxxnum*2;i++){
		dtype2[i]=1;
		xnrm2[i*DIMEN]=1;
		xnrm2[i*DIMEN+1]=0;
		xnrm2[i*DIMEN+2]=0;
	}
	for(i=xxxnum*2;i<xxxnum*3;i++){
		dtype2[i]=1;
		xnrm2[i*DIMEN]=0;
		xnrm2[i*DIMEN+1]=1;
		xnrm2[i*DIMEN+2]=0;
	}
	for(i=xxxnum*3;i<xxxnum*4;i++){
		dtype2[i]=1;
		xnrm2[i*DIMEN]=0;
		xnrm2[i*DIMEN+1]=0;
		xnrm2[i*DIMEN+2]=1;
	}

	double tol2=tol;
	numit = fastlap(&nlhs,&n,&n,x,shape,dtype2,lhstype,rhstype,lhsindex,rhsindex,potfieldxyz,rhsvect,xxx4,xnrm2,&numLev,&numMom,&maxit,&tol2,&fljob);

	long ii;
	for(ii=xxxnum;ii<4*xxxnum;ii++) potfieldxyz[ii]/=spaceunit;


	memcpy(Potential,potfieldxyz,xxxnum*sizeof(double));
	memcpy(FieldX,potfieldxyz+xxxnum,xxxnum*sizeof(double));
	memcpy(FieldY,potfieldxyz+2*xxxnum,xxxnum*sizeof(double));
	memcpy(FieldZ,potfieldxyz+3*xxxnum,xxxnum*sizeof(double));
	//delete[] electrodeIndexLimit;
	//delete[] el;
	delete[] potfieldxyz;
	delete[] xxx4;
	delete[] dtype2;
	delete[] xnrm2;
}



	


void D3world::draw(){}

void D3world::propagateForwardEuler(double x[3],double v[3],double h,double qDivM){
	int i;
	static double lastx[3]={-1e32,-1e32,-1e32};
	double v12[3],pot;
	static double f[3];
	bool samex=true;
	for(i=0;i<3;i++) samex&=(lastx[i]==x[i]);
	if(!samex) calc(x[0],x[1],x[2],pot,f[0],f[1],f[2]);

	for(i=0;i<3;i++){
		x[i]+=h*v[i];
		v[i]+=h*(-qDivM*f[i]/spaceunit);
	}
}


void D3world::propagateForwardVerlet(double x[3],double v[3],double h,double qDivM,bool onedim){
	int i;
	static double lastx[3]={-1e32,-1e32,-1e32};
	double v12[3],pot;
	static double f[3];
	bool samex=true;
	for(i=0;i<3;i++) samex&=(lastx[i]==x[i]);
	if(!samex) calc(x[0],x[1],x[2],pot,f[0],f[1],f[2]);
	if(onedim){
		f[0]=0;
		f[1]=0;
	}// x''=-qDivM*grad(V)/scale^2. 
	// So if x is given in [mm] scale has to be 1000000.
	// 1000 because grad(V) is by a factor 1000 too small this is already included in the fieldcalculation
	// and another 1000 to get from [m] to [mm].

	for(i=0;i<3;i++){
		v12[i]=v[i]+h/2.*(-qDivM*f[i]/spaceunit);
		x[i]+=h*v12[i];
	}
	calc(x[0],x[1],x[2],pot,f[0],f[1],f[2]);
	if(onedim){
		f[0]=0;
		f[1]=0;
	}
	for(i=0;i<3;i++){
		v[i]=v12[i]+h/2.*(-qDivM*f[i]/spaceunit);
		lastx[i]=x[i];
	}
}


void D3world::propagateForwardVerletRotSymX(double x[3],double v[3],double h,double qDivM){
	int i;
	static double lastx[3]={-1e32,-1e32,-1e32};
	double v12[3],pot;
	static double f[3];
	bool samex=true;
	static double r;
	static double fr,f0;

	for(i=0;i<3;i++) samex&=(lastx[i]==x[i]);

	
	if(!samex) {
		r=sqrt(x[1]*x[1]+x[2]*x[2]);
		calc(x[0],r,0,pot,f[0],fr,f0);
		f[1]=fr*x[1]/r;
		f[2]=fr*x[2]/r;
	}
	
	// x''=-qDivM*grad(V)/scale^2. 
	// So if x is given in [mm] scale has to be 1000000.
	// 1000 because grad(V) is by a factor 1000 too small this is already included in the fieldcalculation
	// and another 1000 to get from [m] to [mm].

	for(i=0;i<3;i++){
		v12[i]=v[i]+h/2.*(-qDivM*f[i]/spaceunit);
		x[i]+=h*v12[i];
	}
	r=sqrt(x[1]*x[1]+x[2]*x[2]);
	calc(x[0],r,0,pot,f[0],fr,f0);
	f[1]=fr*x[1]/r;
	f[2]=fr*x[2]/r;
	for(i=0;i<3;i++){
		v[i]=v12[i]+h/2.*(-qDivM*f[i]/spaceunit);
		lastx[i]=x[i];
	}
}

void D3world::SetScalePostCalc(double xscale2,double yscale2,double zscale2){
	xscale=xscale2;
	yscale=yscale2;
	zscale=zscale2;
}
 

void D3world::propagateForwardVerletRotSymY(double x[3],double v[3],double h,double qDivM){
	int i;
	static double lastx[3]={-1e32,-1e32,-1e32};
	double v12[3],pot;
	static double f[3];
	bool samex=true;
	static double r;
	static double fr,f0;
	for(i=0;i<3;i++) samex&=(lastx[i]==x[i]);
	if(!samex) {
		r=sqrt(x[0]*x[0]+x[2]*x[2]);
		calc(0,x[1],r,pot,f0,f[1],fr);
		f[0]=fr*x[0]/r;
		f[2]=fr*x[2]/r;
	}
	// x''=-qDivM*grad(V)/scale^2. 
	// So if x is given in [mm] scale has to be 1000000.
	// 1000 because grad(V) is by a factor 1000 too small this is already included in the fieldcalculation
	// and another 1000 to get from [m] to [mm].

	for(i=0;i<3;i++){
		v12[i]=v[i]+h/2.*(-qDivM*f[i]/spaceunit);
		x[i]+=h*v12[i];
	}
	r=sqrt(x[0]*x[0]+x[2]*x[2]);
	calc(0,x[1],r,pot,f0,f[1],fr);
	f[0]=fr*x[0]/r;
	f[2]=fr*x[2]/r;
	for(i=0;i<3;i++){
		v[i]=v12[i]+h/2.*(-qDivM*f[i]/spaceunit);
		lastx[i]=x[i];
	}
}

void D3world::propagateForwardVerletRotSymZ(double x[3],double v[3],double h,double qDivM){
	int i;
	static double lastx[3]={-1e32,-1e32,-1e32};
	double v12[3],pot;
	static double f[3];
	bool samex=true;
	static double r;
	static double fr,f0;
	for(i=0;i<3;i++) samex&=(lastx[i]==x[i]);
	if(!samex) {
		r=sqrt(x[0]*x[0]+x[1]*x[1]);
		calc(r,0,x[2],pot,fr,f0,f[2]);
		f[0]=fr*x[0]/r;
		f[1]=fr*x[1]/r;
	}
	// x''=-qDivM*grad(V)/scale^2. 
	// So if x is given in [mm] scale has to be 1000000.
	// 1000 because grad(V) is by a factor 1000 too small this is already included in the fieldcalculation
	// and another 1000 to get from [m] to [mm].

	for(i=0;i<3;i++){
		v12[i]=v[i]+h/2.*(-qDivM*f[i]/spaceunit);
		x[i]+=h*v12[i];
	}
	r=sqrt(x[0]*x[0]+x[1]*x[1]);
	calc(r,0,x[2],pot,fr,f0,f[2]);
	f[0]=fr*x[0]/r;
	f[1]=fr*x[1]/r;
	
	
	for(i=0;i<3;i++){
		v[i]=v12[i]+h/2.*(-qDivM*f[i]/spaceunit);
		lastx[i]=x[i];
	}
}
double D3world::calc_slow(double x,double y,double z){
	n=GetAmountOfSubelements();
	amountOfElectrodes=subelement.size();

	

	double *electrodeVoltage=new double[amountOfElectrodes];
	double *dfdn;
	double *dfdnSuperposition=new double[n];

	int col;



	PD3elementList::iterator i;
	int elcnt=0;

	for(i=subelement.begin(),elcnt=0;i!=subelement.end();++i,++elcnt){
		double voltage=((D3electrode*)(*i))->GetVoltage();
		electrodeVoltage[elcnt]=voltage;
		//for(col=electrodeIndexLimit[elcnt];col<electrodeIndexLimit[elcnt+1];col++) f[col]=voltage;
	}


	dfdn=dfdnAll;
		
	for (col=0; col<n; col++)	{	//to call a Fortran routine from C we
		dfdnSuperposition[col]=dfdn[col]*electrodeVoltage[0];			  
	}

	for(elcnt=1;elcnt<amountOfElectrodes;++elcnt){
		dfdn=dfdnAll+n*elcnt;
		
		for (col=0; col<n; col++)	{	//to call a Fortran routine from C we
			dfdnSuperposition[col]+=dfdn[col]*electrodeVoltage[elcnt];			  
		}
	}
	double pot=0;
	for(col=0;col<n;++col){	
		pot+=el[col]->GetPotentialAt(x,y,z)*dfdnSuperposition[col];//the last term of the sum is negledgible
	}

	delete[] electrodeVoltage;
	delete[] dfdnSuperposition;
	return pot;
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////D3SlowWorld::
D3slowworld::D3slowworld():D3element(0,0,0,0,0,0,false,false,NULL),f(NULL),dfdnAll(NULL),el(NULL),electrodeIndexLimit(NULL){};
D3slowworld::~D3slowworld(){
	delete[] f;
	delete[] dfdnAll;
	delete[] el;
	delete[] electrodeIndexLimit;
};
void D3slowworld::insert(D3electrode *el){
	el->parent=this;
	subelement.push_back(el);//triangle at the left bottom corner
};
void D3slowworld::GetListOfBaseElements(PD3element *el,int *electrodeIndexLimit,int &cnt){
	if (isBaseElement) {
		el[cnt++]=this;
		return;
	}
	PD3elementList::iterator i;
	int electrodecnt;
	for(i=subelement.begin(),electrodecnt=0;i!=subelement.end();++i,electrodecnt++){
		electrodeIndexLimit[electrodecnt]=cnt;
		(*i)->GetListOfBaseElements(el,cnt);
		
	}
	electrodeIndexLimit[amountOfElectrodes]=cnt;

};

void relaxdiag(long int n,long int amountOfElectrodes, double *bAll,double *alpha,double *aAll,double *tmp,double errmax){
	long int elcnt,row,col;
	double *f=new double[n];
	double *corr=new double[n];
	for (row=0; row<n; row++)		
	{				
		f[row]=10;
		for(col=0; col<n; col++) {
			if(col!=row) f[row]+=alpha[row+n*col];
		}
	}	
	
	double err;
	int cntdown;
	double mul=1;
	double oldmaxerr=1e99;
	for(elcnt=0;elcnt<amountOfElectrodes;elcnt++){
		double *b,*a;
		a=aAll+n*elcnt;
		for(col=0; col<n; col++) a[col]=1;
		b=bAll+n*elcnt;
		
	
		cntdown=0;
		err=1e99;
		mul=1;
		oldmaxerr=1e99;
		double maxerr;
		do{
			if(cntdown>0) cntdown--;
			for (row=0; row<n; row++)		
			{				
				tmp[row]=0;
				for(col=0; col<n; col++) {
					if(row==col) tmp[row]+=(alpha[row+n*col]+f[row]*mul)*a[col];
					else tmp[row]+=alpha[row+n*col]*a[col];
				}
			}		
			
			double val;
			maxerr=0;
			for(col=0;col<n;col++){
				val=(b[col]+f[col]*mul*a[col]-tmp[col]);
				corr[col]=val/(alpha[col+n*col]+f[col]*mul);
				err=(val>0?val:-val);
				if(err>maxerr) maxerr=err;
			}
			
			if(maxerr/oldmaxerr>1){
				mul*=1.2;
				cntdown=10;
			}
			else{
				for(col=0;col<n;col++) a[col]+=corr[col];
				if((maxerr/oldmaxerr>0.99) && cntdown==0){
					mul/=1.2;
					cntdown=10;
				}		
			}
				
				
			oldmaxerr=maxerr;
			printf("%e %e\n",maxerr,mul);
		}while(maxerr>errmax);
	}
	delete[] f;

}


void D3slowworld::solve(){
	n=GetAmountOfSubelements();
	amountOfElectrodes=subelement.size();
	delete[] electrodeIndexLimit;
	delete[] f;
	delete[] dfdnAll;
	delete[] el;

	electrodeIndexLimit=new int[amountOfElectrodes+1];
	f=new double[n];
	dfdnAll=new double[n*amountOfElectrodes];
	double *alpha=new double[n*n];
	double *beta=new double[n*n];
	double *dfdn;
	el=new PD3element[n];
	int cnt=0;
	GetListOfBaseElements(el,electrodeIndexLimit,cnt);
	int row,col;
	double x,y,z;
	
	for(col=0;col<n;col++){//geht alle Flächenelemente durch
		for(row=0;row<n;row++){//geht alle Reference Positionen durch
			if(row==col){
				alpha[row+n*col]=el[col]->GetSelfPotential()/(4.*pi);
				beta[row+n*col]=el[col]->GetSelfDoubleLayerPotential()/(4.*pi);			
			}
			else
			{
				el[row]->GetReferencePoint(x,y,z);
				if((col==72) && (row==64)){
					int i;
				//	i++;
				}
				alpha[row+n*col]=el[col]->GetPotentialAt(x,y,z)/(4.*pi);
				
				if(alpha[row+n*col]!=alpha[row+n*col]){
					int ii=0;
					ii++;
				}
				beta[row+n*col]=el[col]->GetDoubleLayerPotentialAt(x,y,z)/(4.*pi);
				if(beta[row+n*col]!=beta[row+n*col]){
					int ii=0;
					ii++;
				}
			}
		}
	}
	int elcnt;
	
	long int *pivot=new long int[n];
	for(elcnt=0;elcnt<amountOfElectrodes;elcnt++){
		PD3elementList::iterator i;
		int elcnt2;
		//for(i=subelement.begin(),elcnt2=0;i!=subelement.end();++i,++elcnt2){ 
		//	if(elcnt==elcnt2) ((D3electrode*)(*i))->SetVoltage(1);
		//	else ((D3electrode*)(*i))->SetVoltage(0);
		//}
	
		for(elcnt2=0;elcnt2<amountOfElectrodes;elcnt2++)
			for(col=electrodeIndexLimit[elcnt2];col<electrodeIndexLimit[elcnt2+1];col++) f[col]=(elcnt2==elcnt)?1:0;
		dfdn=dfdnAll+n*elcnt;
		

		for (row=0; row<n; row++)		// to call a Fortran routine from C we 
		{				// have to used the transform of the matrix 
			dfdn[row]=-0.5*f[row];
			for(col=0; col<n; col++) {
				dfdn[row]+=beta[row+n*col]*f[col];
				if(dfdn[row]!=dfdn[row]){
					int ii;
					ii++;
				}
			}
		}					
		
	}

	long int  c1, c2,ok;
	c1=n;				// n: Amount of Subelements and put all numbers we want to pass 
	c2=amountOfElectrodes;    			// to the routine in variables 
#define MATRIXINVERSE
//#define EXPORTDATA
#ifdef EXPORTDATA
	ofstream ofs("c:\\test.txt",ios_base::binary);
	ofs <<n<<"\n";
	ofs <<amountOfElectrodes<<"\n";
	ofs <<n*n*amountOfElectrodes*sizeof(double)<<"\n";
	ofs.write((char*)alpha,n*n*amountOfElectrodes*sizeof(double));
	ofs <<n*amountOfElectrodes*sizeof(double)<<"\n";
	ofs.write((char*)dfdn,n*amountOfElectrodes*sizeof(double));
	ofs.close();

#endif
#ifdef MATRIXINVERSE
	
	// find solution using LAPACK routine SGESV, all the arguments have to 
	// be pointers and you have to add an underscore to the routine name 
#ifdef WIN32
	DGESV(&c1, &c2, alpha, &c1, pivot, dfdnAll, &c1, &ok);
#else
	dgesv_(&c1, &c2, alpha, &c1, pivot, dfdnAll, &c1, &ok);
#endif
#endif
#ifdef LINEARLEASTSQUARE
	double RCOND,LWORKret;
	long int RANK,INFO,LWORK,NLVL,LIWORK,SMLSIZ=25;
	NLVL=long int(log(double(n))/(SMLSIZ+1))+1;
	LIWORK=3*n*NLVL+11*n;
	LWORK=-1;
	double *S=new double[n];
	long int *IWORK=new long int[LIWORK];
	

	RCOND=0.1;

 

	DGELSD(&c1, &c1, &c2, alpha, &c1, dfdnAll, &c1, S, &RCOND,
            &RANK, &LWORKret, &LWORK, IWORK, &INFO);
	LWORK=LWORKret;
	double *WORK=new double[LWORK];
	DGELSD(&c1, &c1, &c2, alpha, &c1, dfdnAll, &c1, S, &RCOND,
            &RANK,  WORK,&LWORK, IWORK, &INFO);
	
	delete[] S;
	delete[] WORK;
	delete[] IWORK;
#endif

#ifdef RELAXDIAG
	double *aAll=new double[n*amountOfElectrodes];
	double *tmp=new double[n];
	relaxdiag(n,amountOfElectrodes, dfdnAll,alpha,aAll,tmp,1e-1); 
	delete[] tmp;
	delete[] dfdnAll;
	dfdnAll=aAll;
#endif

	delete[] pivot;		
	delete[] alpha;
	delete[] beta;
};


double D3slowworld::calc(double x,double y,double z){
	n=GetAmountOfSubelements();
	amountOfElectrodes=subelement.size();

	
	static double *alpha2=new double[n];
	double *beta2=new double[n];
	double *electrodeVoltage=new double[amountOfElectrodes];
	double *dfdn;
	double *dfdnSuperposition=new double[n];

	int col;
	
	for(col=0;col<n;col++){//geht alle Flächenelemente durch
		
			double dummy=el[col]->GetPotentialAt(x,y,z)/(4.*pi);
			//if(fabs(alpha2[col]-dummy)>1e-4){
			//	printf("%d\n",col);
			//}
			alpha2[col]=dummy;

			beta2[col]=el[col]->GetDoubleLayerPotentialAt(x,y,z)/(4.*pi);
		
			dfdnSuperposition[col]=0;
	}

	PD3elementList::iterator i;
	int elcnt=0;
	for(i=subelement.begin(),elcnt=0;i!=subelement.end();++i,elcnt++){
		double voltage=((D3electrode*)(*i))->GetVoltage();
		electrodeVoltage[elcnt]=voltage;
		for(col=electrodeIndexLimit[elcnt];col<electrodeIndexLimit[elcnt+1];col++) f[col]=voltage;
	}

	dfdnbackup=dfdnAll;
	for(elcnt=0;elcnt<amountOfElectrodes;elcnt++){
		dfdn=dfdnAll+n*elcnt;
		
		for (col=0; col<n; col++){		//to call a Fortran routine from C we 
			dfdnSuperposition[col]+=dfdn[col]*electrodeVoltage[elcnt];	
		}
	}
	double pot=0;
	for(col=0;col<n;col++){	
		pot+=-alpha2[col]*dfdnSuperposition[col]+beta2[col]*f[col];//the last term of the sum is negledgible
	}
//	delete[] alpha2;
	delete[] beta2;
	delete[] electrodeVoltage;
	delete[] dfdnSuperposition;
	return pot;
};

void D3slowworld::draw(){}
