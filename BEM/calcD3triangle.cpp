#include "calcD3triangle.h"
#include <cmath>
#include <cstring>

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

void calcD3triangle::triangleint_pot_exyz(int N,double* x0,double* y0,double* z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double q,double* potf){
	double hs=q*sqrt(sqr(x2*y1 - x3*y1 - x1*y2 + x3*y2 + x1*y3 - x2*y3) + sqr(x2*z1 - x3*z1 - x1*z2 + x3*z2 + x1*z3 - x2*z3) + sqr(y2*z1 - y3*z1 - y1*z2 + y3*z2 + y1*z3 - y2*z3));
	
	for (int k = 0; k < N; ++k) {
		double pot=0, ex=0, ey=0, ez=0;
        for(int i=0;i<n;i++){
            double zeta,eta2,xi2;
                
            eta2=eta[i];
            xi2=xi[i];
            zeta=1-xi2-eta2;

            double xq=x1*zeta+x2*xi2+x3*eta2-x0[k];
            double yq=y1*zeta+y2*xi2+y3*eta2-y0[k];
            double zq=z1*zeta+z2*xi2+z3*eta2-z0[k];

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

		potf[4*k] += pot;
		potf[4*k+1] += ex;
		potf[4*k+2] += ey;
		potf[4*k+3] += ez;
	}
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
	//return calc.triangleint(calcD3triangle::G3D,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3);
	
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

void calcD3triangle::G3Danalytic(int N,double* x0,double* y0,double* z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double q,double* pot){
	//first find out if the point of evaluation is in the plane of the triangle or outside the plane
	//therefore calculate norm of triangle
//undokument if numerical sollution should be enforced
	//return calc.triangleint(calcD3triangle::G3D,xp,yp,zp,x1,y1,z1,x2,y2,z2,x3,y3,z3);

	double px,py,pz,ax,ay,az,bx,by,bz,cx,cy,cz,nx,ny,nz,n;
	double a,b,c,psi,cospsi,sinpsi;

	ax=x3-x1;ay=y3-y1;az=z3-z1;
	bx=x2-x1;by=y2-y1;bz=z2-z1;
	cx=x3-x2;cy=y3-y2;cz=z3-z2;
	a=sqrt(sqr(ax)+sqr(ay)+sqr(az));
	b=sqrt(sqr(bx)+sqr(by)+sqr(bz));
	c=sqrt(sqr(cx)+sqr(cy)+sqr(cz));
	if(abs(a)<1e-13) return;
	if(abs(b)<1e-13) return;
	if(abs(c)<1e-13) return;
	cospsi=(ax*bx+ay*by+az*bz)/(a*b);
    if(1.-abs(cospsi)<1e-13) return;
    sinpsi=sqrt(1.-sqr(cospsi));
	psi=acos_safe(cospsi);
	
	for (int k = 0; k < N; ++k) {
        double xp = x0[k], yp = y0[k], zp = z0[k];
        px=-x1+xp;py=-y1+yp;pz=-z1+zp;
        // b x a
        nx=az*by - ay*bz;
        ny=ax*bz - az*bx;
        nz=ay*bx - ax*by;
        n=sqrt(sqr(nx)+sqr(ny)+sqr(nz));
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
		if(abs(ri1)<1e-13) {
			pot[k] += q*(onplane?GetIntL1(xp,yp,zp,x2,y2,z2,x3,y3,z3):GetIntL1OutOfPlane(xp,yp,zp,x2,y2,z2,x3,y3,z3,rp));
			continue;
		}
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
        pot[k] += q*(n1*Trianglei23+n2*Trianglei13+n3*Trianglei12);
    }
}

#ifdef CALCD3_TEST
#include <iostream>
calcD3triangle calc;
struct D3triangle {
	double x1;double y1;double z1;
	double x2;double y2;double z2;
	double x3;double y3;double z3;
	bool inversenorm;
	D3triangle(double nodes[]) : x1(nodes[0]), y1(nodes[1]), z1(nodes[2]), x2(nodes[3]), y2(nodes[4]), z2(nodes[5]), x3(nodes[6]), y3(nodes[7]), z3(nodes[8]), inversenorm(false) {}
	void GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez){
        calc.triangleint_pot_exyz(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm,pot,ex,ey,ez);
    }
    double GetPotentialAt(double x,double y,double z){
        std::cout << calc.G3DdnAnalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm) << std::endl;
        return calc.G3Danalytic(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,inversenorm);
    };
};

int main() {
    double pot, ex, ey, ez;
    double P[3] = {0., 1., 0.};
    double S[9] = {0., 0., 0., 1., 0., 0., 0., 1., 1.};
    D3triangle tri(S);
    tri.GetPotentialAndFieldAt(P[0],P[1],P[2],pot,ex,ey,ez);
    std::cout << pot << ' ' << ex << ' ' << ey << ' ' << ez << std::endl;
    std::cout << tri.GetPotentialAt(P[0],P[1],P[2]) << std::endl;
    return 0;
}
#endif