//Todo heap stack problem with D3ImportedElectrode

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


//Euler angles as defined in goldstein page 146. Best understood as:
//phi, theta psi:
//first rotate object around z axis with angle psi
//then rotate in the objects rotated frame around the local x axis with angle theta
//finalay rotat in the objects rotated frame around the local z axis.
//todo put spaceunits into checksum

#ifndef __BEM_H
#define __BEM_H
#include <complex> //modified std::complex with inversed storage order or real and imag in order to comply with fftw_complex
using namespace std;

#define BEMREVISION 2

const int kBlue = 600;
#include "dl_creationadapter.h"
#include "dl_dxf.h"
#include <float.h>
#include <math.h>

#else
typedef double mxArray;
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

#include <vector>
#include <list>
#include <map>

//#define sqr(x) pow((x),2)
#define sqr(x) ((x)*(x))

const double pi=3.1415926535897932384626433832795028841971693993751058209749445923078164062862090;

// solving the matrix equation A*x=b using LAPACK 
 

#define size_ 3				// dimension of matrix 

double testfn(double x,double y,double z);

void q();


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
	double G3DdnAnalytic(double xp,double yp,double zp,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm);
	void triangleint_pot_exyz(double x0,double y0,double z0,double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3,bool inversenorm,double &pot,double &ex,double &ey,double &ez);

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




//x father class of everything

//D3element needs list of sollutions for each 0 ...0 1 0...0 electrode voltage configuration
//new electrode container class
//new world container class
//MatlabDisplay needs other depth than real element depth
class D3element;
typedef D3element *PD3element;
typedef list<PD3element> PD3elementList;
class D3element{
public:
	/// 
	/// \brief  <b>Uses ray-shooting to determine the norm and should only be used for externally imported geometrical figures.</b>
	///
	/// From Point (x0, y0, z0)  Rays are being shot towards the object in order to 
    /// determine in which direction the normal vector of each surface is pointing.
	/// Generally, all normal vectors can point either to the inside of the object or
	/// to the outside of the object.
	/// We should make sure that all normal vectors are pointing in the same direction!
	/// 
	/// A blue color on each outside-surface indicates a correct norm, which means that
	/// all normal Vectors are pointing in the same direction. The ions that fly by the electrodes see
	/// only the outside-surfaces of the electrodes. If an outside-surface has an incorrect norm, the ion will
	/// not be influenced correctly by the surface. <b>This means that all outside-surfaces should be blue!</b>
	///
	/// A red color on at least one outside-surface indicates that at least one normal vector
	/// is pointing in a different direction than the others. 
	/// In this case we need to modify the parameters x0/y0/z0.
	///
	///
	/// <b>Important for the "shooting point" (x0,y0,z0):</b>
	/// It must not be in the same plane as one of the surfaces!!
	/// Furthermore it must not be inside a geometrical figure!!
	/// Otherwise the Ray-Shooting method will not work properly.
	///
    /// \param x0 x-coordinate of the shooting-point
    /// \param y0 y-coordinate of the shooting-point
	/// \param z0 z-coordinate of the shooting-point
	///	                                             
	virtual	void correctNorm(double x0,double y0,double z0);           	                                                         
	virtual bool IntersectWithRay(double x0,double y0,double z0,double xdir,double ydir,double zdir,double &mul);
	virtual double GetSelfPotential();
	virtual double GetPotentialAt(double x,double y,double z);
	virtual double GetSelfDoubleLayerPotential();
	virtual double GetDoubleLayerPotentialAt(double x,double y,double z);
	virtual void GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez);
	virtual double GetArea();
	virtual void GetCenter(double &x,double &y,double &z);
	D3element(	double x_,double y_,double z_,
				double phi_,double theta_,double psi_,
				bool inversenorm_=false,bool refineable_=true,D3element* parent_=NULL,double epsilon=0,string name="");
	~D3element();
	virtual void GetReferencePoint(double &x,double &y, double &z);
	virtual void init(bool ignorefirst=false);
	///
	/// \brief <b>Refines a geometrical figure.</b>
	///
	/// This means that the original geometrical object will be filled with smaller versions of itself.
	/// Works good with cubic objects. 
	///
	/// Does not work good with other objects like cylinders etc..?? //x  \n \n
    ///
	/// The parameter length determines the size of the smaller objects.\n \n
    ///  Example:
	/// Cube with sides (1x1x1) and parameter length 0.1 => 10 small cubes will be fit along each
	/// side, so the cube will then contain 10*10*10=1000 small cubes.
	///
	/// Press ctrl+w to see the lattice model and get an impression of the refinement.
	/// Press ctrl+r to see the original model again
	///
	/// \param length size of the smaller objects that fill up the incident object
	///
	virtual void refine(double length);	
	virtual void refine(double length,int num); //?? unterschied zu oben? -> evtl um sachen zu überspringen? //x
	int GetAmountOfSubelements();
	int PrintAmountOfSubelements();
	void GetListOfBaseElements(PD3element *el,int &cnt);
	virtual void rotate(double phi_,double theta_,double psi_);
	virtual void shift(double xs,double ys,double zs);
	virtual void rotate2(double phi_,double theta_,double psi_,double x,double y,double z,double &x2,double &y2,double &z2);
	virtual void rotateinv(double phi_,double theta_,double psi_,double x,double y,double z,double &x2,double &y2,double &z2);

	//	virtual int Get3DMatlab( mxArray *plhs[],int cnt=-1);
	D3element *parent;
	void Add(PD3element el);
	virtual void GetTriangle(double *A,double *B,double *C,double *COL);
	virtual void GetRectangle(double *A,double *B,double *C,double *D,double *COL);
	virtual void SetNormTowards(double x,double y,double z,bool towards);
	int Color;
	virtual void createNewSubelements(double length);
	virtual string getName();
	virtual void setName(string name2);
	virtual bool GetInversenorm();
	string name;
	list<D3element*> subelement;
protected:
	double epsilon;
	//TGeoHMatrix h;
	virtual void deleteSubelements();
	double x;
	double y;
	double z;
	double phi;
	double theta;
	double psi;
	bool isBaseElement;
	
	bool refineable;
	bool inversenorm;

};



class D3electrode;

class D3electrode:public D3element{
public:
	///
	/// \brief <b>Constructor for an electrode. </b>
	///
	/// The father class is D3element, which means an electrode consists of at least one 3-dimensional element and additional attributes such as voltage etc.	
	///
	D3electrode();	
	//D3electrode(const D3electrode &cp);
	~D3electrode();
	///
	/// \brief <b>Allocates an D3element to an electrode, multiple elements can be allocated to one electrode.</b>
	///
	/// \param el D3element that is being allocated to the electrode
	///
	void insert(PD3element el);	
	///
	/// \brief Sets the voltage for the electrode. Unit: Volt.
	///
	/// \param voltage_ voltage for the electrode, unit: Volt.
	///
	void SetVoltage(double voltage_);
    /// 
	/// \brief Gives the voltage of the electrode. Unit: Volt.
	///
	double GetVoltage();
	///
	/// \brief Eliminates small numerical deviations from perfect symmetry.
	///
	/// If you need an electrode that is symmetric to a certain plane or direction, you have to place it properly at first. \n
	/// <b>The method will now only correct small numerical deviations from perfect symmetry!</b> It cannot provide large-scale translations or rotations!\n \n
    /// 
	/// <b>Additionally to the methode SetSymmetry , this method also provides to set the current electrode in perfect symmetry with the one specified by *syme12</b>\n \n
    ///
	/// xsym=true means symmetric with respect to x=0 plane \n
	/// ysym=true means symmetric with respect to y=0 plane \n
	/// zsym=true means symmetric with respect to z=0 plane \n
	/// diagmirror=true means symmetric with respect of exchange of x and y (if z is set as axis in SymmetrizeCharges())  \n
	/// diagmirror=true means symmetric with respect of exchange of y and z (if x is set as axis in SymmetrizeCharges())  \n
	/// diagmirror=true means symmetric with respect of exchange of x and z (if y is set as axis in SymmetrizeCharges())  \n
	/// symel2 is the other electrode unequal to *this that the symmetry refers to \n
	/// rotsym means that symmetric with respect to rotsym/totalsym*360 degree rotation (set in a loop if you have rotational symmetry)\n
	///
	/// \param xsym  true, if symmetric with respect to x=0 plane
	/// \param ysym  true, if symmetric with respect to y=0 plane
	/// \param zsym  true, if symmetric with respect to z=0 plane
	/// \param diagmirror true, if see above
	/// \param symel2 other electrode unequal to *this that the symmetry refers to
	/// \param rotsym  ?? int ??
	/// \param totalrot ??
	///
	void SetSymmetryWith(bool xsym,bool ysym,bool zsym,bool diagmirror,D3electrode *symel2,int rotsym=0,int totalrot=0);
	///
	/// \brief Eliminates small numerical deviations from perfect symmetry.
	///
	/// If you need an electrode that is symmetric to a certain plane or direction, you have to place it properly at first. \n
	/// <b>The method will now only correct small numerical deviations from perfect symmetry!</b> It cannot provide large-scale translations or rotations! \n \n
    ///
	/// xsym=true means symmetric with respect to x=0 plane \n
	/// ysym=true means symmetric with respect to y=0 plane  \n
	/// zsym=true means symmetric with respect to z=0 plane \n
	/// diagmirror=true means symmetric with respect of exchange of x and y (if z is set as axis in SymmetrizeCharges())  \n
	/// diagmirror=true means symmetric with respect of exchange of y and z (if x is set as axis in SymmetrizeCharges())  \n
	/// diagmirror=true means symmetric with respect of exchange of x and z (if y is set as axis in SymmetrizeCharges())  \n
	/// rotsym means that symmetric with respect to rotsym/totalsym*360 degree rotation (set in a loop if you have rotational symmetry)  \n 
	///
    ///
	/// \param xsym  true, if symmetric with respect to x=0 plane
	/// \param ysym  true, if symmetric with respect to y=0 plane
	/// \param zsym  true, if symmetric with respect to z=0 plane
	/// \param diagmirror true, if see above
	/// \param rotsym  ?? int ??
	/// \param totalrot2 ??
	///
	void SetSymmetry(bool xsym,bool ysym,bool zsym,bool diagmirror,int rotsym=0,int totalrot2=0);
	D3electrode *GetSymmetryWith(int num);
	int GetTotalrot();

	void SetCardinalNumber(int num);
	int GetCardinalNumber();
protected:
	int cardinalnum;
	double voltage;
	D3electrode **symel;
	int totalrot;
};

class D3sortedelement{
public:
	D3sortedelement():done(false){}
	~D3sortedelement(){}
	friend bool SmallerThan(D3sortedelement &a,D3sortedelement &b,double eps){
		if(abs(a.axis-b.axis)<eps) {
			if(a.radius<b.radius) return true;
			else return false;
		}
		else if(a.axis<b.axis) return true;
		else return false;
		
	}
	bool done;
	double axis;//Major sort value
	double radius;//minor sort value
	int index; //für dfdnAll und x
	D3electrode *electrodeptr; //Ptr auf die das Element beinhaltende Elektrode
};

class D3sorter{
private:
	D3sortedelement *a;
	int n;
public:
    void sort(D3sortedelement *arr,int size,double epsi)
    {
        a=arr;
        n=size;
		eps=epsi;
        quicksort(0, n-1);
    }

	void quicksort (int l, int r)
    {
       if(r>l){
		int i=l-1, j=r, tmp;
		if(r-l > 3){ //Median of three
			int m=l+(r-l)/2;
			D3sortedelement &x=a[m];
			if(SmallerThan(x,a[l],eps)) exchange(m,l);
			if(SmallerThan(a[r],a[l],eps)) exchange(r,l);
			else if(SmallerThan(a[m],a[r],eps)) exchange(r,m);
		}

		for(;;){
			while(SmallerThan(a[++i],a[r],eps));
			while(SmallerThan(a[r],a[--j],eps) && j>i);
			if(i>=j) break;
			exchange(i,j);
		}
		exchange(i,r);
		

		quicksort( l, i-1);
		quicksort( i+1, r);

	   }
		
    };

    void exchange(int i, int j)
    {
		D3sortedelement t=a[i];
        a[i]=a[j];
        a[j]=t;
    }
protected:
	double eps;


};



//x
class D3ImportedElectrodes:public DL_CreationAdapter{
public:
    ///
	/// \brief Can be used to import an Electrode from an external .dxf file, that was created by AutoCAD.
	///
	/// The lines below show an example for proper use of the method \n \n
	///
    /// 
	/// D3ImportedElectrodes *impel=new D3ImportedElectrodes(); \n
	///	TString importfilename;  \n
	///	logfile.AbsoluteFileName("./Ansaug_endcap.dxf",importfilename);\n
	///	if(!impel->Import(importfilename)) return 0; \n \n
    ///
	///	endcapl=& (impel->FindElectrode("endcap")); \n
	///	endcapl->correctNorm(0,0,-10); \n
	///	endcapl->refine(0.02); \n \n \n
    ///
	///       so genug??
    ///
    D3ImportedElectrodes();
	~D3ImportedElectrodes();
    ///
    /// \brief  Loads the .dxf file into memory.
    ///
	/// Displays the filename given by  file  , in a command line output. \n
	/// Returns false and gives an error message if the file cannot be opened.  \n \n
	///
	/// so genug??
    ///
    /// \param file name of the .dxf file
    ///
	bool Import(const char* file,bool ignore3DFace_=false,bool ignorePolyline_=false);
	
    virtual void addLayer(const DL_LayerData& data);
    virtual void addPoint(const DL_PointData& data);
    virtual void addLine(const DL_LineData& data);
    virtual void addArc(const DL_ArcData& data);
    virtual void addCircle(const DL_CircleData& data);
    virtual void addPolyline(const DL_PolylineData& data);
    virtual void addVertex(const DL_VertexData& data);
	virtual void add3DFace(const DL_3DFaceData &data);
	virtual void addBlockRecord(const DL_BlockRecordData &data);
	//virtual void add3dFace(const DL_3dFaceData &data);
	//virtual void addBlock(const DL_BlockData &data);
	void endSequence();
	///
	/// \brief Iterates over all electrodes and displays their names in a command line output.
	///
	void ListElectrodes(); 
	D3electrode &FindElectrode(const char* name);
    void printAttributes();
	map<string,D3electrode> electrode;
protected:
	int Model_Space_Handle;
	bool ignore_Model_Space_Handle;
	bool ignore3DFace;
	bool ignorePolyline;
	bool importingPolyline;
	int vertexcnt;
	double x[4]; 
	double y[4];
	double z[4];
};
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


class D3triangle:public D3element{
public:
	virtual bool IntersectWithRay(double x0,double y0,double z0,double xdir,double ydir,double zdir,double &mul);
	static bool IntersectWithRay(double ax,double ay,double az,double bx,double by,double bz,double cx,double cy,double cz,double xdir,double ydir,double zdir,double &mul);

	virtual void SetNormTowards(double x,double y,double z,bool towards);
	virtual void GetReferencePoint(double &x,double &y, double &z);
	virtual double GetArea();
	virtual void GetCenter(double &x,double &y,double &z);
	virtual double GetSelfPotential();
	virtual inline double GetPotentialAt(double x,double y,double z);
	virtual void GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez);
	virtual double GetSelfDoubleLayerPotential();
	virtual inline double GetDoubleLayerPotentialAt(double x,double y,double z);
	///
    /// \brief Creates a triangle from 2 2-dimensional vectors, 1 3-dimensional center vector and 3 Euler-angles
	///
	/// The parameters x,y,z determine the 1. corner of the triangle. The Euler-rotations are also performed with respect to this corner.\n
	/// (x,y,z) can also be interpreted as a center of a coordinate system. \n \n
	/// The 2. corner is determined by (xa_ + x, ya_ + y). The 3. corner is determined by (xb_ + x, yb_ + y).
	/// The ordering does not matter!! \n \n
	///
	/// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles. \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n \n
	/// They can be used to generate z-coordinates for the triangle \n \n
	///
	/// The Euler-rotations occur in the following order: \n \n
	///    
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n
	///

	D3triangle(	double xa_,double ya_,double xb_,double yb_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL,bool stripeit=true);
	///
	/// \brief Creates a triangle from 3 3-dimensional vectors
	///
	///  (x1_,y1_z1_), (x2_,y2_,z2_), (x3_,y3_,z3_) define the corners of the triangle. The ordering does not matter!!
	///

	D3triangle(double x1_,double y1_,double z1_,double x2_,double y2_,double z2_,double x3_,double y3_,double z3_,
			bool inversenorm=false,PD3element parent=NULL,bool stripeit=true);

	~D3triangle();	
	virtual void rotate(double phi_,double theta_,double psi_);
	virtual void shift(double xs,double ys,double zs);
	virtual void GetTriangle(double *A,double *B,double *C,double *COL);
	virtual void createNewSubelements(double length);
	
	double x1;double y1;double z1;
	double x2;double y2;double z2;
	double x3;double y3;double z3;
protected:
	void Stripeit(double length);
	void Newtriangles(double length);
	bool stripeit;
	double xa;double ya;
	double xb;double yb;
};


class D3rectangle:public D3element{
public:
	virtual bool IntersectWithRay(double x0,double y0,double z0,double xdir,double ydir,double zdir,double &mul);
	virtual void SetNormTowards(double x,double y,double z,bool towards);
	virtual void GetReferencePoint(double &x,double &y, double &z);
	virtual double GetArea();
	virtual void GetCenter(double &x,double &y,double &z);
	virtual double GetSelfPotential();
	virtual inline double GetPotentialAt(double x,double y,double z);
	virtual double GetSelfDoubleLayerPotential();
	virtual inline double GetDoubleLayerPotentialAt(double x,double y,double z);
	virtual void GetPotentialAndFieldAt(double x,double y,double z,double &pot,double &ex,double &ey,double &ez);
	virtual void GetRectangle(double *A,double *B,double *C,double *D,double *COL);
	virtual void rotate(double phi_,double theta_,double psi_);
	virtual void shift(double xs,double ys,double zs);
	///
	/// \brief Creates a rectangle from 2 2-dimensional vectors, 1 3-dimensional center vector and 3 Euler angles
    ///
	/// The parameters x,y,z determine the 1. corner of the rectangle. The Euler-rotations are also performed with respect to this corner. \n
	/// (x,y,z) can also be interpreted as a center of a coordinate system.  \n
	/// The 2. corner is determined by (xa_ + x, ya_ + y). \n
	/// The 3. corner is determined by (xb_ + x, yb_ + y).   \n 
	/// The 4. corner is determined by (xa_ + xb_ + x, ya_ + yb_ + y). This represents a vector addition of the vectors for the first 2 corners.  \n \n
	/// 
	/// <b>The ordering of the corners does not matter!!</b>   \n \n
	///
    /// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n
	///  The Center (x,y,z) specifies the center of the Euler-rotations \n \n
	/// <b>They can be used to generate z-coordinates for the rectangle</b>  \n \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n
	///
	D3rectangle(double xa_,double ya_,double xc_,double yc_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL,bool subdivideInSqares=false);
    ///
	/// \brief Creates a rectangle from 3 2-dimensional vectors, 1 center vector and 3 Euler angles
	///
	/// The parameters x,y,z determine the 1. corner of the rectangle. The Euler-rotations are also performed with respect to this corner. \n
	/// (x,y,z) can also be interpreted as a center of a coordinate system.  \n
	/// The 2. corner is determined by (xa_ + x, ya_ + y). \n
	/// The 3. corner is determined by (xb_ + x, yb_ + y).   \n 
	/// The 4. corner is not being calculated automatically, it is specified via (xc_,yc_)! \n \n
	///
	/// <b>The ordering of the corners does not matter!!</b>  \n \n
	///
	/// <b>If the 3 corner-vectors would result in a non-convex object only a triangle will be drawn!</b>  \n \n
	///
    /// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.    \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n \n
	/// <b>They can be used to generate z-coordinates for the rectangle</b>  \n \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n
	///
	D3rectangle(double xa_,double ya_,double xb_,double yb_,double xc_,double yc_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL);
    ///
    /// \brief Creates a rectangle from coordinates for each corner
    /// 
	/// Parameters x1_, y1_, z1_ define the coordinates of the first corner \n
	/// The remaining paramters define the other corners 
	///
	/// \n \n wofür der string???--> nur wegen überladung?
    ///
	                                                                                        
	D3rectangle(char *str,double x1_,double y1_,double z1_,double x2_,double y2_,double z2_,double x3_,double y3_,double z3_,double x4_,double y4_,double z4_,
			bool inversenorm=false,PD3element parent=NULL);
	~D3rectangle(){};
	virtual void createNewSubelements(double length);
	double x1;double y1;double z1;
	double x2;double y2;double z2;
	double x3;double y3;double z3;
	double x4;double y4;double z4;
protected:
	double xa;double ya;
	double xb;double yb;
	double xc;double yc;
};


class D3box:public D3element{
public:
    ///
	/// \brief Creates a 3-dimensional symmetric box
	/// 
	/// Initial orientation of the Box and coordinate system: \n
	/// x-axis: towards the viewer \n
	/// y-axis: upwards \n
	/// z-axis: to the right \n \n
	///
	/// x,y,z: Coordinates of the boxes origin   \n \n
    ///
	/// x1_: defines length of box in x-direction  \n
	/// y1_: defines length of box in y-direction   \n
	/// z1_: defines length of box in z-direction   \n
	/// negative lengths are also possible    \n
	///
	///  x + x1_, y + y1_, z + z1_ define the edges of the box   \n \n
	///
    ///
	/// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n
	///
	D3box(double xl_,double yl_,double zl_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL):
	                xl(xl_),yl(yl_),zl(zl_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent){
						subelement.push_back(new D3rectangle(xl,0,0,yl,0,0,0,0,0,0,!inversenorm,this));//triangle at the left bottom corner
						subelement.push_back(new D3rectangle(xl,0,0,yl,0,0,zl,0,0,0,inversenorm,this));//triangle at the right bottom corner
						subelement.push_back(new D3rectangle(xl,0,0,zl,0,0,0,0,90,0,inversenorm,this));//triangle at the right bottom corner
						subelement.push_back(new D3rectangle(xl,0,0,zl,0,yl,0,0,90,0,!inversenorm,this));//triangle at the right bottom corner
						subelement.push_back(new D3rectangle(yl,0,0,zl,0,0,0,0,90,90,!inversenorm,this));//triangle at the right bottom corner
						subelement.push_back(new D3rectangle(yl,0,0,zl,xl,0,0,0,90,90,inversenorm,this));//triangle at the right bottom corner

					};
	~D3box(){};	

protected:
	double xl;double yl;double zl;
};


class D3icosahedron:public D3element{
public:
	D3icosahedron(double edgelength,double refinement,double x,double y,double z,double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL):r(r),refinement(refinement),
					D3element(x,y,z,phi,theta,psi,inversenorm,false,parent){
						double e=edgelength/2.,t=(1.+sqrt(5.))/2.*edgelength/2.;
						double VV[12][3]={{0,e,t},{e,t,0},{t,0,e},
											 {0,e,-t},{e,-t,0},{t,0,-e},
											 {0,-e,t},{-e,t,0},{-t,0,e},
											 {0,-e,-t},{-e,-t,0},{-t,0,-e}};
						V=&VV;
							
						
						triangle( 0, 1, 2, true);
						triangle(1,2,5);
						triangle(4,6,2,true);

						triangle(2,5,4,true);
						triangle(0,2,6,true);
						triangle(4,5,9,true);

						triangle(4,9,10,true);
						triangle(4,6,10);//??
						triangle(10,6,8);

						triangle(8,6,0);
						triangle(8,7,0,true);
						triangle(0,1,7);

						triangle(1,3,5,true);
						triangle(3,5,9);
						triangle(11,10,9,true);

						triangle(11,10,8);
						triangle(11,8,7);
						triangle(11,7,3);
						triangle(11,3,9);

						triangle(1 ,3,7);

					
					};//x  --> eventuell fehler drin??, änderung von edgelength und refinement ändert nichts
	~D3icosahedron(){};	
	void triangle(int a,int b,int c,bool inversenorm=false){
		subelement.push_back(new D3triangle((*V)[a][0],(*V)[a][1],(*V)[a][2],
											(*V)[b][0],(*V)[b][1],(*V)[b][2],
											(*V)[c][0],(*V)[c][1],(*V)[c][2],
											inversenorm,this,false));
	}
protected:
	double r;
	double refinement;
	double (*V)[12][3];
};


class D3sphere:public D3element{
public:
	/// 
    /// \brief Creates a sphere with radius r.
	///
	/// The sphere is centered at (x,y,z). \n \n
    ///
	/// The parameter refinement gives the order of magnitude for the small triangles which assemble the sphere. \n
	/// Press ctrl+w to see the lattice model and get an impression of the refinement.\n
	/// Press ctrl+r to see the original model again. \n
	///
	/// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n
	///
	D3sphere(double r,double refinement,double x,double y,double z,double phi_,double theta_,double psi_,bool inversenorm=false,PD3element parent=NULL):r(r),refinement(refinement),
					D3element(x,y,z,phi,theta,psi,inversenorm,false,parent){
						double edgelength=4./sqrt(10.+2.*sqrt(5.));//radius 1
						D3icosahedron *ico=new D3icosahedron(edgelength,0,0,0,0,0,0,inversenorm,NULL);
						ico->refine(refinement/r);
						int n=ico->GetAmountOfSubelements();
						PD3element *el=new PD3element[n];
						int cnt=0;
						ico->GetListOfBaseElements(el,cnt);

						double a[3],b[3],c[3],anorm,bnorm,cnorm,color[3];
						int col;
						for(col=0;col<n;col++){//geht alle Flaechenelemente durch
							if(dynamic_cast<D3triangle *>(el[col])) {
								el[col]->GetTriangle(a,b,c,color);
								anorm=sqrt(sqr(a[0])+sqr(a[1])+sqr(a[2]));
								bnorm=sqrt(sqr(b[0])+sqr(b[1])+sqr(b[2]));
								cnorm=sqrt(sqr(c[0])+sqr(c[1])+sqr(c[2]));
								subelement.push_back(new D3triangle(x+a[0]/anorm*r,y+a[1]/anorm*r,z+a[2]/anorm*r,
																	x+b[0]/bnorm*r,y+b[1]/bnorm*r,z+b[2]/bnorm*r,
																	x+c[0]/cnorm*r,y+c[1]/cnorm*r,z+c[2]/cnorm*r,
																	inversenorm,this,false));
							}
						}
						delete[] el;
						delete ico;
					
					};


	~D3sphere(){};	
protected:
	double r;
	double refinement;
};

class D3disk:public D3element{
public:
	///
	/// \brief Creates a disk with radius r_.
	///
	/// The disk is centered at (x,y,z). \n  \n 
    ///
	/// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n \n \n
	///
	/// For phi_ = theta_ = psi_ = 0, the disk lies in the x,y-plane and there is no extension in z-direction \n \n
	///
	/// The parameter sectors_ provides an option for refinement. It is set to 32 by initialization but can also be set manually. \n
	///	A value of sectors_ = 8 creates a disk with 8 edges. sectors_ is not limited to powers of two or even numbers!! 5 or 11 for example work as well.
	///

	D3disk(double r_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,int sectors_=32,bool subdivide=true,bool inversenorm=false,PD3element parent=NULL):
					r(r_),sectors(sectors_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent){
						int i;
						double step=2*pi/double(sectors);
						double dr=r*step;
						int subdivr=(r/dr);
						if(subdivr<1||!subdivide){
							for(i=0;i<sectors;i++){
								subelement.push_back(new D3triangle(r*cos(step*double(i)),r*sin(step*double(i)),r*cos(step*double(i+1)),r*sin(step*double(i+1)),0,0,0,0,0,0,inversenorm,this));//triangle at the left bottom corner
							}
						}
						else{
							for(i=0;i<sectors;i++){
								subelement.push_back(new D3triangle(r/double(subdivr)*cos(step*double(i)),r/double(subdivr)*sin(step*double(i)),r/double(subdivr)*cos(step*double(i+1)),r/double(subdivr)*sin(step*double(i+1)),0,0,0,0,0,0,inversenorm,this));//triangle at the left bottom corner
								int ring;
								for(ring=1;ring<subdivr;ring++){
									double xa=r/double(subdivr)*double(ring)*cos(step*double(i));
									double ya=r/double(subdivr)*double(ring)*sin(step*double(i));
									double xd=r/double(subdivr)*double(ring)*cos(step*double(i+1));
									double yd=r/double(subdivr)*double(ring)*sin(step*double(i+1));
									double xb=r/double(subdivr)*double(ring+1)*cos(step*double(i));
									double yb=r/double(subdivr)*double(ring+1)*sin(step*double(i));
									double xc=r/double(subdivr)*double(ring+1)*cos(step*double(i+1));
									double yc=r/double(subdivr)*double(ring+1)*sin(step*double(i+1));
									subelement.push_back(new D3triangle(xb-xa,yb-ya,xc-xa,yc-ya,xa,ya,0,0,0,0,inversenorm,this));
									subelement.push_back(new D3triangle(xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,0,0,0,0,inversenorm,this));
								}						
							}
							

						}
					};
					
	~D3disk(){};	
protected:
	double r;
	double sectors;

};

class D3cylinder:public D3element{
public:
	///
	/// \brief Creates a cylinder with radius r_ and length l_.
	///
	/// The bottom of the cylinder is centered at (x,y,z). \n \n
    ///
    /// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n \n \n
	///
	/// For phi_ = theta_ = psi_ = 0, the bottom of the cylinder lies in the x,y-plane and the z-extension correspond to the length l_.
	///
	/// The parameter sectors_ provides an option for refinement. It is set to 12 by initialization but can also be set manually.
	///	A value of sectors_ = 8 creates a cylinder with 8 edges on the bottom for example. sectors_ is not limited to powers of two or even numbers!! 5 or 11 works as well.
	///
	D3cylinder(double r_,double l_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,int sectors_=12,bool subdivide=true,bool inversenorm=false,PD3element parent=NULL):
					r(r_),l(l_),sectors(sectors_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent){
						int i;
						double step=2*pi/double(sectors);
						double dr=r*sqrt(2.*(1.-cos(step)));
						int subdivr=(r/dr);
						if(subdivr<=1||!subdivide){
							for(i=0;i<sectors;i++){
								subelement.push_back(new D3triangle(r*cos(step*double(i)),r*sin(step*double(i)),r*cos(step*double(i+1)),r*sin(step*double(i+1)),0,0,0,0,0,0,!inversenorm,this));//triangle at the left bottom corner
								subelement.push_back(new D3triangle(r*cos(step*double(i)),r*sin(step*double(i)),r*cos(step*double(i+1)),r*sin(step*double(i+1)),0,0,l,0,0,0,inversenorm,this));//triangle at the left bottom corner
								double xb=r*cos(step*double(i));
								double yb=r*sin(step*double(i));
								subelement.push_back(new D3rectangle(dr,0,0,l,xb,yb,0,0,90,90.+180./pi*step*(double(i)+0.5),inversenorm,this,true));

							}
							
						}
						else{
							for(i=0;i<sectors;i++){
								//center circle of front faces
								subelement.push_back(new D3triangle(r/double(subdivr)*cos(step*double(i)),r/double(subdivr)*sin(step*double(i)),r/double(subdivr)*cos(step*double(i+1)),r/double(subdivr)*sin(step*double(i+1)),0,0,0,0,0,0,!inversenorm,this));//triangle at the left bottom corner
								subelement.push_back(new D3triangle(r/double(subdivr)*cos(step*double(i)),r/double(subdivr)*sin(step*double(i)),r/double(subdivr)*cos(step*double(i+1)),r/double(subdivr)*sin(step*double(i+1)),0,0,l,0,0,0,inversenorm,this));//triangle at the left bottom corner
								int ring;
								double xb,yb;
								for(ring=1;ring<subdivr;ring++){
									double xa=r/double(subdivr)*double(ring)*cos(step*double(i));
									double ya=r/double(subdivr)*double(ring)*sin(step*double(i));
									double xd=r/double(subdivr)*double(ring)*cos(step*double(i+1));
									double yd=r/double(subdivr)*double(ring)*sin(step*double(i+1));
									xb=r/double(subdivr)*double(ring+1)*cos(step*double(i));
									yb=r/double(subdivr)*double(ring+1)*sin(step*double(i));
									double xc=r/double(subdivr)*double(ring+1)*cos(step*double(i+1));
									double yc=r/double(subdivr)*double(ring+1)*sin(step*double(i+1));

									//sectors of front faces as rectangles
									subelement.push_back(new D3rectangle(xb-xa,yb-ya,xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,0,0,0,0,!inversenorm,this));
									subelement.push_back(new D3rectangle(xb-xa,yb-ya,xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,l,0,0,0,inversenorm,this));
									//sectors of front faces as triangles (are documented out)
									//	subelement.push_back(new D3triangle(xb-xa,yb-ya,xc-xa,yc-ya,xa,ya,0,0,0,0,!inversenorm,this));
									//	subelement.push_back(new D3triangle(xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,0,0,0,0,!inversenorm,this));

									//	subelement.push_back(new D3triangle(xb-xa,yb-ya,xc-xa,yc-ya,xa,ya,l,0,0,0,inversenorm,this));
									//	subelement.push_back(new D3triangle(xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,l,0,0,0,inversenorm,this));				
								}			
								//outer side element
								subelement.push_back(new D3rectangle(dr,0,0,l,xb,yb,0,0,90,90.+180./pi*step*(double(i)+0.5),inversenorm,this,true));
							}
							

						}
					};				
	
	~D3cylinder(){};	
protected:
	double r;
	double l;
	double sectors;

};


class D3tube:public D3element{
public:
    ///
	/// \brief  Creates a tube with outer radius r_, inner radius r2_ and length l_.
	///
	/// r2_ should be smaller than r_ for the correct norm. \n \n
	///
	/// The bottom of the tube is centered at (x,y,z). \n \n
    ///
    /// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n \n \n
	///
	/// For phi_ = theta_ = psi_ = 0, the bottom of the tube lies in the x,y-plane and the z-extension correspond to the length l_. \n \n
	///
	/// The parameter sectors_ provides an option for refinement. It is set to 12 by initialization but can also be set manually. \n
	///	A value of sectors_ = 8 creates a cylinder with 8 edges on the bottom for example. sectors_ is not limited to powers of two or even numbers!! 5 or 11 works as well.
	///
	D3tube(double r_,double r2_,double l_,
				double x,double y,double z,
				double phi_,double theta_,double psi_,int sectors_=12,bool subdivide=true,bool inversenorm=false,PD3element parent=NULL):
					r(r_),r2(r2_),l(l_),sectors(sectors_),
					D3element(x,y,z,phi_,theta_,psi_,inversenorm,false,parent){
						if(!subdivide){
							cout <<"subdivide false not yet implemented"<<endl;					
						}
						subdivide=true;
						int i;
						double step=2*pi/double(sectors);
						double dr=r*sqrt(2.*(1.-cos(step)));
						int subdivr=((r-r2)/dr);
						if(subdivr<=1) subdivr=2;
						if(subdivr<=1||!subdivide){
							for(i=0;i<sectors;i++){
								subelement.push_back(new D3triangle(r*cos(step*double(i)),r*sin(step*double(i)),r*cos(step*double(i+1)),r*sin(step*double(i+1)),0,0,0,0,0,0,!inversenorm,this));//triangle at the left bottom corner
								subelement.push_back(new D3triangle(r*cos(step*double(i)),r*sin(step*double(i)),r*cos(step*double(i+1)),r*sin(step*double(i+1)),0,0,l,0,0,0,inversenorm,this));//triangle at the left bottom corner
								double xb=r*cos(step*double(i));
								double yb=r*sin(step*double(i));
								subelement.push_back(new D3rectangle(dr,0,0,l,xb,yb,0,0,90,90.+180./pi*step*(double(i)+0.5),inversenorm,this,true));

							}
							
						}
						else{
							for(i=0;i<sectors;i++){
								//center circle of front faces
								//subelement.push_back(new D3triangle(r/double(subdivr)*cos(step*double(i)),r/double(subdivr)*sin(step*double(i)),r/double(subdivr)*cos(step*double(i+1)),r/double(subdivr)*sin(step*double(i+1)),0,0,0,0,0,0,!inversenorm,this));//triangle at the left bottom corner
								//subelement.push_back(new D3triangle(r/double(subdivr)*cos(step*double(i)),r/double(subdivr)*sin(step*double(i)),r/double(subdivr)*cos(step*double(i+1)),r/double(subdivr)*sin(step*double(i+1)),0,0,l,0,0,0,inversenorm,this));//triangle at the left bottom corner
								int ring;
								double xb,yb;
								for(ring=0;ring<subdivr;ring++){
									double rr=r2+(r-r2)/double(subdivr)*double(ring);
									double rr2=r2+(r-r2)/double(subdivr)*double(ring+1);
									double xa=rr*cos(step*double(i));
									double ya=rr*sin(step*double(i));
									double xd=rr*cos(step*double(i+1));
									double yd=rr*sin(step*double(i+1));
									xb=rr2*cos(step*double(i));
									yb=rr2*sin(step*double(i));
									double xc=rr2*cos(step*double(i+1));
									double yc=rr2*sin(step*double(i+1));

									//sectors of front faces as rectangles
									subelement.push_back(new D3rectangle(xb-xa,yb-ya,xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,0,0,0,0,!inversenorm,this));
									subelement.push_back(new D3rectangle(xb-xa,yb-ya,xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,l,0,0,0,inversenorm,this));
									//sectors of front faces as triangles (are documented out)
									//	subelement.push_back(new D3triangle(xb-xa,yb-ya,xc-xa,yc-ya,xa,ya,0,0,0,0,!inversenorm,this));
									//	subelement.push_back(new D3triangle(xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,0,0,0,0,!inversenorm,this));

									//	subelement.push_back(new D3triangle(xb-xa,yb-ya,xc-xa,yc-ya,xa,ya,l,0,0,0,inversenorm,this));
									//	subelement.push_back(new D3triangle(xc-xa,yc-ya,xd-xa,yd-ya,xa,ya,l,0,0,0,inversenorm,this));				
								}			
								//outer side element
								subelement.push_back(new D3rectangle(dr,0,0,l,xb,yb,0,0,90,90.+180./pi*step*(double(i)+0.5),inversenorm,this,true));
								//inner side element
							
								xb=r2*cos(step*double(i));
								yb=r2*sin(step*double(i));
								double xb2=r2*cos(step*double(i+1));
								double yb2=r2*sin(step*double(i+1));
								double dr2=sqrt((xb-xb2)*(xb-xb2)+(yb-yb2)*(yb-yb2));
								subelement.push_back(new D3rectangle(dr2,0,0,l,xb,yb,0,0,90,90.+180./pi*step*(double(i)+0.5),!inversenorm,this,true));
							}
							

						}
					};	
	~D3tube(){};	
protected:
	double r;
	double r2;
	double l;
	double sectors;

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Electron charge     
#define eee 1.602e-19
//Mass of Calcium
#define mmm (40.078*1.66e-27)

class D3world:public D3element{
	
public:
	bool rangeerror;
	D3world(unsigned long num1,unsigned long num2,unsigned long num3,unsigned long num4,unsigned long d,unsigned long m,unsigned long y); //x
	/// aus der cpp:
	/// num is registration key, d is day, m is month, y is year
	///
	///   ???  

	D3world(const char *_cachefilename,double tol=0.001,int maxit=32,int numMom=2,int numLev=4,double spaceunit=0.001,int segmentation=1000);
	//D3world(const char *_cachefilename,double tol=0.001,int maxit=32,int numMom=2,int numLev=4,double spaceunit=0.001,int segmentation=100000);
	~D3world();
	void insert(D3electrode *el);//x
	void GetListOfBaseElements(PD3element *el,int *electrodeIndexLimit,int &cnt);
	void SymmetrizeCharges(int axis=1,double epsilon=0.00001,bool ignoremirror=false);//x=1;y=2;z=3 call it after solve!//x
	void solve();//x  ??
	/// 
	/// \brief Computes matrices that are required for the calculations of the potentials
	///
	double calc(double x,double y,double z);//x
	void calc(double xxx,double y,double z,double &pot,double &feldx,double &feldy,double &feldz);//x

	void calc(double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz);///< Caching data for access with calc(x,y,z) or calc(x,y,z,&pot,&feldx,&feldy,&feldz). //x
	void calc(int xxxnum,double *xxx,double *Potential);///< Calculating only p			otentials. //x
	void calc(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ);///<Calculating potentials and fields. //x
	double calc_slow(double x,double y,double z); //x
	void calc_slow(double xmin,double xmax,int nx,double ymin,double ymax,int ny,double zmin,double zmax,int nz); //x
	void calc_slow(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ);///<Calculating potentials and fields. //x
	void SetScalePostCalc(double xscale2,double yscale2,double zscale2); 
    bool IsEqualSurfaceElement(int amountOfVertices,double eps,double *x,double *X);
	void draw(); //x
	void Dcentroid(int shape,double * pc,double* xcout);
	void save(char *fname);
	bool load(char *fname);//true if saved data is compatible
	int cut(int n,int max);
	void savecalc(char *fname);
	bool loadcalc(char *fname);//true if saved data is compatible
	void AssignColors();
	unsigned long update_adler32double(unsigned long old, double *buf, unsigned long len,int ignorebytes=3);
	void RefreshChecksum();
	void exportGeometry(const char *fname);//x
	/// Propagates particle at position x[3] one timestep h forward. 
	/// \param x[3] contains x,y,z position at current timestep. After call contains the position at next timestep.
	/// \param v[3] contains x,y,z velocity at current timestep.
	/// \param h is timestep in seconds.
	/// \param qDivM is charge of particle divided by mass.
	void propagateForwardVerlet(double x[3],double v[3],double h,double qDivM=eee/mmm,bool onedim=false);//x
	void propagateForwardEuler(double x[3],double v[3],double h,double qDivM=eee/mmm);//x
	void propagateForwardVerletRotSymX(double x[3],double v[3],double h,double qDivM=eee/mmm);//x
	void propagateForwardVerletRotSymY(double x[3],double v[3],double h,double qDivM=eee/mmm);//x
	void propagateForwardVerletRotSymZ(double x[3],double v[3],double h,double qDivM=eee/mmm);//x
	void calc2(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ);///<Calculating potentials and fields.  //x
	void calc_slow2(int xxxnum,double *xxx,double *Potential,double *FieldX,double *FieldY,double *FieldZ);///<Calculating potentials and fields. //x

protected:
	double xscale;
	double yscale;
	double zscale;
	int totalrot;
	void exch(double &a,double &b);
	void RotMirrorSurfaceElement(int axis,int sym,double *xx);
	bool refreshchecksum;
	int currentel;
	int segmentation;
	double spaceunit;
	char cachefilename[2000];
//caching stuff
	bool docache;
	double xmin;
	double xmax;
	int nx;
	double ymin;
	double ymax;
	int ny;
	double zmin;
	double zmax;
	int nz;
	double* potcache;
	double* feldxcache;
	double* feldycache;
	double* feldzcache;
//end of caching stuff
  int size;
  int nlhs;
  int nrhs;
  int numMom;
  int numLev;
  int i;
  int j;
  char *shapechar;
  double *x;
  double *poten;
  double *dbydnpotenAll;
  double *dbydbpoten;
  double *xcoll;
  double *xnrm;
  double *lhsvect;
  double *rhsvect;
  int *shape;
  int *type;
  int *dtype;
  int *rhstype;
  int *lhstype;
  int *rhsindex;
  int *lhsindex;
  int job;
  int fljob; 
  double error;
  double max_diri;
  double ave_diri;
  double max_neum;
  double ave_neum;
  double cnt_diri;
  double cnt_neum; 
  /* Set the tolerance and max iterations for GMRES called by fastlap. */
  double tol;
  int maxit;
  int numit;
  unsigned long checksum;
  unsigned long checksumcalc;
  
  double *f;
	double *dfdnAll;
	int *electrodeIndexLimit;
	PD3element *el;
	D3electrode **electrodes;
	int n;
	int amountOfElectrodes;
	unsigned long update_adler32(unsigned long old, unsigned char *buf, unsigned long len);

};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class D3slowworld:public D3element{
public:
	D3slowworld();
	~D3slowworld();
	void insert(D3electrode *el);
	void GetListOfBaseElements(PD3element *el,int *electrodeIndexLimit,int &cnt);
	void solve();
	double calc(double x,double y,double z);
	void draw();
	
protected:
	double *f;
	double *dfdnAll;
	int *electrodeIndexLimit;
	PD3element *el;
	int n;
	int amountOfElectrodes;
};

extern  calcD3triangle calc;


class D3thicktriangle:public D3element{
public:
	///
	/// \brief Creates a thick triangle from 2 2-dimensional vectors, 1 3-dimensional center vector and 3 Euler-angles
	///
	/// The parameter htriangle determines the thickness of the triangle. For phi_ = theta_= psi_ = 0 the thickness corresponds to the extension in z-direction. \n \n
	///
	/// The parameters x,y,z determine the 1. corner of the triangle. The Euler-rotations are also performed with respect to this corner. \n
	/// (x,y,z) can also be interpreted as a center of a coordinate system. \n
	/// The 2. corner is determined by (xa_ + x, ya_ + y). The 3. corner is determined by (xb_ + x, yb_ + y). \n \n
	/// The ordering does not matter!!
	///
    /// The angles phi_, theta_ and  psi_ can be declared in degrees and represent Euler-angles.   \n \n
	/// The Center (x,y,z) specifies the center of the Euler-rotations \n
	///
	/// The Euler-rotations occur in the following order:  \n   \n 
	/// 
	/// phi_:..........Rotation around z-axis...................................x -> x'......y  -> y'.......z -> z \n
	/// theta_:......Rotation around the rotated x-axis (x'):.........x'-> x'......y' -> y''......z -> z' \n
	/// psi_:..........Rotation around z'-axis:.................................x'-> x''......y''-> y'''.....z'-> z' \n \n \n
    ///
	D3thicktriangle(double xa,double ya,double xb,double yb,double htriangle,
				double x,double y,double z,
				double phi,double theta,double psi,bool inversenorm=false,PD3element parent=NULL):
					xa(xa),ya(ya),xb(xb),yb(yb),htriangle(htriangle),
					D3element(x,y,z,phi,theta,psi,inversenorm,false,parent){
						subelement.push_back(new D3triangle(xa,ya,xb,yb,0,0,0,0,0,0,!inversenorm,this));
						subelement.push_back(new D3triangle(xa,ya,xb,yb,0,0,htriangle,0,0,0,inversenorm,this));

					double l1=sqrt(sqr(xa)+sqr(ya));
						double w1=atan(ya/xa);
						subelement.push_back(new D3rectangle(l1,0,0,htriangle,0,0,0,0,90,w1/pi*180,!inversenorm,this));
						
						double l2=sqrt(sqr(xb)+sqr(yb));
						double w2=atan(yb/xb);
						subelement.push_back(new D3rectangle(l2,0,0,htriangle,0,0,0,0,90,w2/pi*180,inversenorm,this));
				
						double l3=sqrt(sqr(xb-xa)+sqr(yb-ya));
						double w3=atan((xb-xa)/(yb-ya));

						subelement.push_back(new D3rectangle(l3,0,0,htriangle,xa,ya,0,0,90,90-w3/pi*180,inversenorm,this));
				
					

	};

	~D3thicktriangle(){};	
protected:
	double xa;double ya;double xb;double yb;double htriangle;
};