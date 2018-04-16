#ifndef GMALL_H
#define GMALL_H

#ifndef GMHEADER_H
#define GMHEADER_H
#endif

#ifndef M_PI
#define M_PI	((2)*(acos(0.0)))
#endif

#define GM_STARTYEAR 	1900
#define RAD2DEG(rad)	((rad)*(180.0L/M_PI))
#define DEG2RAD(deg)	((deg)*(M_PI/180.0L))
#define ATanH(x)	(0.5 * log((1 + x) / (1 - x)))
#define MU_0		4*M_PI / 10000000
#define R_e		6.371 * 1000000

#define TRUE		((int)1)
#define FALSE		((int)0)


typedef struct {
			int Day;
			int Month;
			int Year;
			double DecimalYear;
			int DayNumber;
			} GMtype_Date;

typedef struct {
			double lambda;// longitude
			double phi; // geodetic latitude
			double HeightAboveEllipsoid; // height above the ellipsoid (HaE)
			} GMtype_CoordGeodetic;

typedef struct {
			double r; //radial distance
			double phig; //latitude
			double lambda; //longitude
			} GMtype_CoordSpherical;
typedef struct {
			double lambda;
			double phi;
			} GMtype_CoordDipole;

typedef struct {
			double x;
			double y;
			double z;
			} GMtype_CoordCartesian;

typedef struct {
			double a; /*semi-major axis of the ellipsoid*/
			double b; /*semi-minor axis of the ellipsoid*/
			double fla; /* flattening */
			double epssq; /*first eccentricity squared */
			double eps; /* first eccentricity */
			double re;/* mean radius of  ellipsoid*/
			} GMtype_Ellipsoid;

typedef struct {
			int degree;
			double coef[35];
			} GMtype_Polynomial;

typedef struct {
			double h1;
			double g1;
			double g0;
			} GMtype_Model;

typedef struct {
			double M; /*Magnetic Moment*/
			double phi;
			double lambda;
			} GMtype_Pole;
			
typedef struct{
			int rows;
			int columns;
			double element[35][35];/*rows then columns*/
			} GMtype_Matrix;

typedef struct{
			int size;
			double element[35];
			} GMtype_Data;


//GM Cord functions
void GM_CartesianToSpherical(GMtype_CoordCartesian CoordCartesian, GMtype_CoordSpherical *CoordSpherical);
void GM_CORD (GMtype_CoordGeodetic location, GMtype_Date *date, GMtype_Ellipsoid Ellip, GMtype_Data g0d, GMtype_Data g1d, GMtype_Data h1d, GMtype_CoordDipole *CoordDipole);
int GM_DateToYear(GMtype_Date *Date);
void GM_EarthCartToDipoleCartCD(GMtype_Pole Pole, GMtype_CoordCartesian EarthCoord, GMtype_CoordCartesian *DipoleCoords);
void GM_GeodeticToSpherical(GMtype_Ellipsoid Ellip, GMtype_CoordGeodetic CoordGeodetic, GMtype_CoordSpherical *CoordSpherical);
void GM_GetUserInput(GMtype_CoordGeodetic *location, GMtype_Date *date);
void GM_PoleLocation(GMtype_Model Model, GMtype_Pole *Pole);
void GM_PrintUserData(GMtype_CoordGeodetic location, GMtype_Date date, GMtype_CoordDipole DipLocation);
void GM_ScanIGRF(GMtype_Data *G0, GMtype_Data *G1, GMtype_Data *H1);
void GM_SetEllipsoid(GMtype_Ellipsoid *Ellip);
void GM_SphericalToCartesian(GMtype_CoordSpherical CoordSpherical, GMtype_CoordCartesian *CoordCartesian);
void GM_TimeAdjustCoefs(GMtype_Date Date, GMtype_Data g0d, GMtype_Data g1d, GMtype_Data h1d, GMtype_Model *Model);

//Math Functions
double GM_DotProduct(GMtype_Data VectorA, GMtype_Data VectorB);
double GM_LinearInterpolation(double x1, double x2, double y1, double y2, double x);
void GM_LUDecomposition(GMtype_Matrix A, GMtype_Matrix *L, GMtype_Matrix *U, GMtype_Matrix *P);
void GM_LUSolve(GMtype_Matrix L, GMtype_Matrix U, GMtype_Matrix P, GMtype_Matrix *x, GMtype_Matrix b);
double GM_MatDet(GMtype_Matrix Matrix);
void GM_MatMultiply(GMtype_Matrix MatrixA, GMtype_Matrix MatrixB, GMtype_Matrix *MatrixC);
void GM_MatInverse(GMtype_Matrix Matrix, GMtype_Matrix *InvertedMatrix);
void GM_MatTranspose(GMtype_Matrix Matrix, GMtype_Matrix *TMatrix);
double GM_Mean(GMtype_Data Data);
void GM_Median(GMtype_Data Data, double *upper, double *lower);
void GM_PolyFit(GMtype_Data DataX, GMtype_Data DataY, GMtype_Polynomial *Polynomial);
double GM_Pow(double x, int y);
void GM_PrintMatrix(GMtype_Matrix X);
double GM_SolvePolynomial(GMtype_Polynomial Polynomial, double x);
void GM_Sort(GMtype_Data *Data);
double GM_StandardDeviation(GMtype_Data Data);
void GM_Swap(double *x, double *y);
void GM_SwapRows(GMtype_Matrix *Matrix, int Row1, int Row2);

#endif
