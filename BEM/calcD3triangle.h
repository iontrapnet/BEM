#ifdef _WIN32
#undef M_PI
#define M_PI   3.14159265358979323846
#pragma warning(disable : 4800)
#endif

//#define sqr(x) pow((x),2)
#define sqr(x) ((x)*(x))

typedef double trianglefn(double,double,double,double,double,double);

class calcD3triangle{
public:
	calcD3triangle();
	~calcD3triangle();
	double GetIntL1(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3);
	double GetIntL1OutOfPlane(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double rp);
	double GetIntL1divR3(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double rp);
	//triangleint(&testfn,0,0,0, 1,0,0, 0,1,0);

	double triangleint(trianglefn *fn,double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3);
	static double G3D(double x0,double y0,double z0,double x1,double y1,double z1);
	static double dG3Ddx(double x0,double y0,double z0,double x1,double y1,double z1);
	static double dG3Ddy(double x0,double y0,double z0,double x1,double y1,double z1);
	static double dG3Ddz(double x0,double y0,double z0,double x1,double y1,double z1);
	double G3Danalytic(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm);
	void G3Danalytic(int N,double* x0,double* y0,double* z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double q,double* pot);
	double G3DdnAnalytic(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm);
	void triangleint_pot_exyz(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm,double &pot,double &ex,double &ey,double &ez);
	void triangleint_pot_exyz(int N,double* x0,double* y0,double* z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,double q,double* potf);

protected:
	void gauleg(double x1, double x2, double x[], double w[], int n); //from nr c p152 c4-5.pdf
	//Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
	//arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
	//Legendre n-point quadrature formula.

	void gauslegsimplextriangle(double xi[],double eta[],double c[],int n);

	double EPS; //EPS is the relative precision.
	double *xi;
	double *eta;
	double *w;
	int numquad;
	int n;
};