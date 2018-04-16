#include "lane/GM_SubLibrary.hpp"



void GM_CartesianToSpherical(GMtype_CoordCartesian CoordCartesian, GMtype_CoordSpherical *CoordSpherical)
{
    /*This function converts a point from Cartesian coordinates into spherical coordinates*/
    double X, Y, Z;

    X = CoordCartesian.x;
    Y = CoordCartesian.y;
    Z = CoordCartesian.z;

    CoordSpherical->r = sqrt(X * X + Y * Y + Z * Z);
    CoordSpherical->phig = RAD2DEG(asin(Z / (CoordSpherical->r)));
    CoordSpherical->lambda = RAD2DEG(atan2(Y, X));
} /*GM_CartesianToSpherical*/

void GM_CORD (GMtype_CoordGeodetic location, GMtype_Date *date, GMtype_Ellipsoid Ellip, GMtype_Data g0d, GMtype_Data g1d, GMtype_Data h1d, GMtype_CoordDipole *CoordDipole)
{
    /*This function can be used in a wrapper file to call all of the necessary functions to convert from geocentric coordinates to geomagnetic coordinates for a given time*/
    GMtype_CoordSpherical CoordSpherical, DipoleSpherical;
    GMtype_CoordCartesian CoordCartesian, DipoleCartesian;
    GMtype_Model Model;
    GMtype_Pole Pole;
    GM_DateToYear(date);
    GM_GeodeticToSpherical(Ellip, location, &CoordSpherical);
    GM_SphericalToCartesian(CoordSpherical,  &CoordCartesian);
    GM_TimeAdjustCoefs(*date, g0d, g1d, h1d, &Model);
    GM_PoleLocation(Model, &Pole);
    GM_EarthCartToDipoleCartCD(Pole, CoordCartesian, &DipoleCartesian);
    GM_CartesianToSpherical(DipoleCartesian, &DipoleSpherical);
    CoordDipole->phi = DipoleSpherical.phig;
    CoordDipole->lambda = DipoleSpherical.lambda;
} /*GM_CORD*/

int GM_DateToYear(GMtype_Date *CalendarDate)
{
    /*This function converts a Date into Decimal years.  It was taken from the WMM SubLibrary on the NGDC website*/
    int temp = 0; /*Total number of days */
    int MonthDays[13];
    int ExtraDay = 0;
    int i;
    if((CalendarDate->Year%4 == 0 && CalendarDate->Year%100 != 0) || CalendarDate->Year%400 == 0)
        ExtraDay = 1;
    MonthDays[0] = 0;
    MonthDays[1] = 31;
    MonthDays[2] = 28 + ExtraDay;
    MonthDays[3] = 31;
    MonthDays[4] = 30;
    MonthDays[5] = 31;
    MonthDays[6] = 30;
    MonthDays[7] = 31;
    MonthDays[8] = 31;
    MonthDays[9] = 30;
    MonthDays[10] = 31;
    MonthDays[11] = 30;
    MonthDays[12] = 31;

    for(i = 1; i <= CalendarDate->Month; i++)
        temp+=MonthDays[i-1];
    temp+=CalendarDate->Day;
    CalendarDate->DayNumber = temp;
    CalendarDate->DecimalYear = CalendarDate->Year + (temp-1)/(365.0 + ExtraDay);
    return TRUE;
}  /*GM_DateToYear*/

void GM_EarthCartToDipoleCartCD(GMtype_Pole Pole, GMtype_CoordCartesian EarthCoord, GMtype_CoordCartesian *DipoleCoords)
{
    /*This function converts from Earth centered cartesian coordinates to dipole centered cartesian coordinates*/
    double X, Y, Z, CosPhi, SinPhi, CosLambda, SinLambda;
    CosPhi = cos(DEG2RAD(Pole.phi));
    SinPhi = sin(DEG2RAD(Pole.phi));
    CosLambda = cos(DEG2RAD(Pole.lambda));
    SinLambda = sin(DEG2RAD(Pole.lambda));
    X = EarthCoord.x;
    Y = EarthCoord.y;
    Z = EarthCoord.z;

    /*These equations are taken from a document by Wallace H. Campbell*/
    DipoleCoords->x = X * CosPhi * CosLambda + Y * CosPhi * SinLambda - Z * SinPhi;
    DipoleCoords->y = -X * SinLambda + Y * CosLambda;
    DipoleCoords->z = X * SinPhi * CosLambda + Y * SinPhi * SinLambda + Z * CosPhi;
} /*GM_EarthCartToDipoleCartCD*/


void GM_GeodeticToSpherical(GMtype_Ellipsoid Ellip, GMtype_CoordGeodetic CoordGeodetic, GMtype_CoordSpherical *CoordSpherical)
{
    double CosLat, SinLat, rc, xp, zp; /*all local variables */
    /*
    ** Convert geodetic coordinates, (defined by the WGS-84
    ** reference ellipsoid), to Earth Centered Earth Fixed Cartesian
    ** coordinates, and then to spherical coordinates.
    */

    CosLat = cos(DEG2RAD(CoordGeodetic.phi));
    SinLat = sin(DEG2RAD(CoordGeodetic.phi));

    /* compute the local radius of curvature on the WGS-84 reference ellipsoid */

    rc = Ellip.a / sqrt(1.0 - Ellip.epssq * SinLat * SinLat);

    /* compute ECEF Cartesian coordinates of specified point (for longitude=0) */

    xp = (rc + CoordGeodetic.HeightAboveEllipsoid) * CosLat;
    zp = (rc*(1.0 - Ellip.epssq) + CoordGeodetic.HeightAboveEllipsoid) * SinLat;

    /* compute spherical radius and angle lambda and phi of specified point */

    CoordSpherical->r = sqrt(xp * xp + zp * zp);
    CoordSpherical->phig = RAD2DEG(asin(zp / CoordSpherical->r));     /* geocentric latitude */
    CoordSpherical->lambda = CoordGeodetic.lambda;                   /* longitude */
} /*GM_GeodeticToSpherical*/

void GM_GetUserInput(GMtype_CoordGeodetic *location, GMtype_Date *date)
{
    /*This function is used in the GMCORD program to get the users input data through the console*/
    char buffer[20];
    int flag; /*This is a loop control that allows the program to check the validity of the users input*/
    flag = 1;
    while(flag == 1)
    {
        printf("ENTER the geographic decimal latitude (+ North; - South)? ");
        fgets(buffer, 20, stdin);
        sscanf(buffer, "%lf", &location->phi);
        flag = 0;
        if(location->phi <-90 || location->phi > 90)
        {
            printf("Input must be in decimal degrees, between -90 and 90\n");
            flag = 1;
        }
        strcpy(buffer, "");
    }
    flag = 1;
    while(flag == 1)
    {
        printf("ENTER the geographic longitude (+ East; - West)? ");
        fgets(buffer, 20, stdin);
        sscanf(buffer, "%lf", &location->lambda);
        flag = 0;
        if(location->lambda < -180 || location->lambda > 360)
        {
            printf("Input must be in decimal degrees, between -180 and 360\n");
            flag = 1;
        }
        strcpy(buffer, "");
    }
    flag = 1;
    while(flag == 1)
    {
        printf("ENTER the analysis year (1900 to 2015)? ");
        fgets(buffer, 20, stdin);
        sscanf(buffer, "%d", &date->Year);
        flag = 0;
        if(date->Year < 1900 || date->Year > 2015)
        {
            printf("Input must be in integer years between 1900 and 2015\n");
            flag = 1;
        }
        strcpy(buffer, "");
    }
    flag = 1;
    while(flag == 1)
    {
        printf("Month number (1 to 12)? ");
        fgets(buffer, 20, stdin);
        sscanf(buffer, "%d", &date->Month);
        flag = 0;
        if(date->Month < 1 || date->Month > 12)
        {
            printf("Input must be an integer from 1 to 12.\n");
            flag = 1;
        }
        strcpy(buffer, "");
    }
    flag = 1;
    while(flag == 1)
    {
        printf("Day in month (1 to 31)? ");
        fgets(buffer, 20, stdin);
        sscanf(buffer, "%d", &date->Day);
        flag = 0;
        if(date->Day < 1 || date->Day > 31)
        {
            printf("Input must be an integer from 1 to 31.\n");
            flag = 1;
        }
        strcpy(buffer, "");
    }
} /*GM_GetUserInput*/

void GM_PoleLocation(GMtype_Model Model, GMtype_Pole *Pole)
{
    /*This function finds the location of the north magnetic pole in spherical coordinates.  The equations are
    **from Wallace H. Campbell's Introduction to Geomagnetic Fields*/

    Pole->phi = RAD2DEG(-atan(sqrt(Model.h1 * Model.h1 + Model.g1 * Model.g1)/Model.g0));
    Pole->lambda = RAD2DEG(atan(Model.h1/Model.g1));
} /*GM_PoleLocation*/

void GM_PrintUserData(GMtype_CoordGeodetic location, GMtype_Date date, GMtype_CoordDipole DipLocation)
{
    /*This function prints the users data to the console*/
    printf("Centered Geomagnetic Dipole location at %lf, %lf:\n", location.phi, location.lambda);
    if(DipLocation.phi >= 0)
        printf("\t+%.2lf deg geomag lat. ", DipLocation.phi);
    else
        printf("\t%.2lf deg geomag lat. ", DipLocation.phi);
    if(DipLocation.lambda >= 0)
        printf("and +%.2lf deg E. geomag long.\n", DipLocation.lambda);
    else
        printf("and %.2lf deg E. geomag long.\n", DipLocation.lambda);
} /*GM_PrintUserData*/

void GM_ScanIGRF(GMtype_Data *G0, GMtype_Data *G1, GMtype_Data *H1)
{
    /*This function scans inputs G0, G1, and H1 of the IGRF table into 3 data arrays*/
    int i;
    double temp;
    char buffer[200];
    FILE *IGRF;
    IGRF = fopen("IGRF.tab", "r");
    G0->size = 25;
    G1->size = 25;
    H1->size = 25;
    for( i = 0; i < 4; i++)
    {
        fgets(buffer, 200, IGRF);
    }
    fscanf(IGRF, "g 1 0 %lf ", &G0->element[0]);
    for(i = 1; i <= 22; i++)
    {
        fscanf(IGRF ,"%lf ", &G0->element[i]);
    }
    fscanf(IGRF ,"%lf\n", &temp);
    G0->element[23] = temp * 5 + G0->element[22];
    G0->element[24] = G0->element[23] + 5 * temp;
    fscanf(IGRF, "g 1 1 %lf ", &G1->element[0]);
    for(i = 1; i <= 22; i++)
    {
        fscanf( IGRF, "%lf ", &G1->element[i]);
    }
    fscanf(IGRF, "%lf\n", &temp);
    G1->element[23] = temp * 5 + G1->element[22];
    G1->element[24] = temp * 5 + G1->element[23];
    fscanf(IGRF, "h 1 1 %lf ", &H1->element[0]);
    for(i = 1; i <= 22; i++)
    {
        fscanf( IGRF, "%lf ", &H1->element[i]);
    }
    fscanf(IGRF, "%lf\n", &temp);
    H1->element[23] = temp * 5 + H1->element[22];
    H1->element[24] = temp * 5 + H1->element[23];
} /*GM_ScanIGRF*/

void GM_SetEllipsoid(GMtype_Ellipsoid *Ellip)
{
    /*This function sets the WGS84 reference ellipsoid to its default values*/
    Ellip->a	=			6378.137; /*semi-major axis of the ellipsoid in */
    Ellip->b	=			6356.7523142;/*semi-minor axis of the ellipsoid in */
    Ellip->fla	=			1/298.257223563;/* flattening */
    Ellip->eps	=			sqrt(1- ( Ellip->b *	Ellip->b) / (Ellip->a * Ellip->a ));  /*first eccentricity */
    Ellip->epssq	=			(Ellip->eps * Ellip->eps);   /*first eccentricity squared */
    Ellip->re	=			6371.2;/* Earth's radius */
} /*GM_SetEllipsoid*/

void GM_SphericalToCartesian(GMtype_CoordSpherical CoordSpherical, GMtype_CoordCartesian *CoordCartesian)
{
    /*This function converts spherical coordinates into Cartesian coordinates*/
    double CosPhi = cos(DEG2RAD(CoordSpherical.phig));
    double SinPhi = sin(DEG2RAD(CoordSpherical.phig));
    double CosLambda = cos(DEG2RAD(CoordSpherical.lambda));
    double SinLambda = sin(DEG2RAD(CoordSpherical.lambda));

    CoordCartesian->x = CoordSpherical.r * CosPhi * CosLambda;
    CoordCartesian->y = CoordSpherical.r * CosPhi * SinLambda;
    CoordCartesian->z = CoordSpherical.r * SinPhi;
} /*GM_SphericalToCartesian*/

void GM_TimeAdjustCoefs(GMtype_Date Date, GMtype_Data g0d, GMtype_Data g1d, GMtype_Data h1d, GMtype_Model *Model)
{
    /*This function calls GM_LinearInterpolation for the coefficients to estimate the value of the
    **coefficient for the given date*/
    int index;
    double x;
    index = (Date.Year - GM_STARTYEAR) / 5;
    x = (Date.DecimalYear - GM_STARTYEAR) / 5;
    Model->g0 = GM_LinearInterpolation(index, index+1, g0d.element[index], g0d.element[index+1], x);
    Model->g1 = GM_LinearInterpolation(index, index+1, g1d.element[index], g1d.element[index+1], x);
    Model->h1 = GM_LinearInterpolation(index, index+1, h1d.element[index], h1d.element[index+1], x);
} /*GM_TimeAdjustCoefs*/

//--------------------------------------------------------------------------------------------
//General math functions
//--------------------------------------------------------------------------------------------

double GM_DotProduct(GMtype_Data VectorA, GMtype_Data VectorB)
{
    /*This function takes the dot product of the elements stored in VectorA and VectorB*/
    double Total = 0;
    int i = 0;

    /*SUM(A_i * B_i)*/
    for(i = 0; i < VectorA.size; i++)
    {
        Total += VectorA.element[i] * VectorB.element[i];
    }
    return Total;
} /*GM_DotProduct*/

double GM_LinearInterpolation(double x1, double x2, double y1, double y2, double x)
{
    /*This function takes a linear interpolation between two given points for x*/
    double weight, y;
    weight  = (x - x1) / (x2 - x1);
    y = y1 * (1 - weight) + y2 * weight;
    return y;
}/*GM_LinearInterpolation*/

void GM_LUDecomposition(GMtype_Matrix A, GMtype_Matrix *L, GMtype_Matrix *U, GMtype_Matrix *P)
{
    /*This function LU decomposes a given matrix A and stores it in two new matricies L and U.
    **It uses the Gaussian Elimination with partial pivoting algorithm from Golub & Van Loan, Matrix Computations, **Algorithm 3.4.1.*/
    int i, j, pivot, k, N;
    double max;

    N = A.rows;
    L->rows = A.rows;
    L->columns = A.columns;
    U->rows = A.rows;
    U->columns = A.columns;
    for(i = 0; i < N - 1; i++)
    {
        max  = fabs(A.element[i][i]);
        pivot = i;

        for(j = i+1; j < N; j++)
        {
            if(fabs(A.element[j][i]) > max)
            {
                max = fabs(A.element[j][i]);
                pivot = j;
            }
        }

        if (pivot != i)
        {
            GM_SwapRows (&A, i, pivot);
            GM_SwapRows (P, i, pivot);
        }

        if(A.element[i][i] != 0.0)
        {
            for (j = i + 1; j < N; j++)
            {
                A.element[j][i] = A.element[j][i] / A.element[i][i];
                for (k = i + 1; k < N; k++)
                {
                    A.element[j][k] = A.element[j][k] - A.element[j][i] * A.element[i][k];
                }
            }
        }
    }
    for(i = 0; i < N; i++)
    {
        L->element[i][i] = 1;
        for(j = 0; j < N; j++)
        {
            if(j < i)
            {
                L->element[i][j] = A.element[i][j];
            }
            else
            {
                U->element[i][j] = A.element[i][j];
            }
        }
    }
} /*GM_LUDecomposition*/

void GM_LUSolve(GMtype_Matrix L, GMtype_Matrix U, GMtype_Matrix P, GMtype_Matrix *x, GMtype_Matrix b)
{
    /*This function solves an LU decomposed matrix equation LUx = b for x*/
    GMtype_Matrix y;
    int N, i, j;
    /*Apply permutation*/
    GM_MatMultiply(P, b, &b);
    N = U.rows;
    x->rows = b.rows;
    x->columns = 1;
    /*y: = Ux, thus Ly = b*/
    /*solve for y by forward substitution*/
    y.rows = b.rows;
    y.columns = 1;
    y.element[0][0] = b.element[0][0] / L.element[0][0];
    for(i = 1; i < N; i++)
    {
        y.element[i][0] = b.element[i][0] / L.element[i][i];
        for(j = 0; j < i; j++)
        {
            y.element[i][0] += -L.element[i][j] * y.element[j][0] / L.element[i][i];
        }
    }
    /*Ux = y*/
    /*Solve for x by backward substitution*/

    x->element[N-1][0] = y.element[N-1][0] / U.element[N-1][N-1];
    for(i = N - 2; i >= 0; i--)
    {
        x->element[i][0] = y.element[i][0] / U.element[i][i];
        for(j = i+1; j < N; j++)
        {
            x->element[i][0] += -U.element[i][j] * x->element[j][0] / U.element[i][i];
        }
    }
} /*GM_LUSolve*/

double GM_MatDet(GMtype_Matrix Matrix)
{
    /*This is a recursive function that calculates the determinant of the matrix*/
    double value = 0, SubDet;
    GMtype_Matrix TempMatrix;
    int i, j, k, l, m;

    if(Matrix.rows != Matrix.columns)
    {
        return 0;
    }
    if(Matrix.rows == 1)
    {
        value = Matrix.element[0][0];
        return value;
    }
    else if(Matrix.rows == 2)
    {
        value = Matrix.element[0][0] * Matrix.element[1][1] - Matrix.element[0][1] * Matrix.element[1][0];
        return value;
    }
    else if(Matrix.rows > 2)
    {
        TempMatrix.rows = Matrix.rows - 1;
        TempMatrix.columns = Matrix.columns - 1;
        for(m=0; m < Matrix.columns; m++)
        {
            /*This loop finds the subdeterminant of the matrix, multiplies it by a the element
            **in a column and sums or subtracts it to/from the total*/
            k = 0;
            for(i=0; i < Matrix.rows; i++)
            {
                l = 0;
                for(j=0; j < Matrix.columns; j++)
                {
                    if(i!=0 && j!=m)
                    {
                        TempMatrix.element[k][l] = Matrix.element[i][j];
                        l++;
                    }
                }
                if(i!=0) k++;
            }
            SubDet = GM_MatDet(TempMatrix);
            if (m%2 == 1) SubDet = SubDet * -1;
            value += SubDet * Matrix.element[0][m];
        }
        return value;
    }
    else
    {
        return 0;
    }
} /*GM_MatDet*/

void GM_MatInverse(GMtype_Matrix Matrix, GMtype_Matrix *InvertedMatrix)
{
    /*This function calculates the inverse of an input matrix.  Error starts to become noticeable around 6x6 matricies.
    **The function starts to take very large amounts of time beyond 8x8 matricies.  LU decomposition is probably a
    **better approach for large matricies*/
    int x, y, m, n, i, j;
    double MatDet;
    GMtype_Matrix Temp;

    MatDet = GM_MatDet(Matrix);
    Temp.rows = Matrix.rows - 1;
    Temp.columns = Matrix.columns - 1;
    InvertedMatrix->rows = Matrix.rows;
    InvertedMatrix->columns = Matrix.columns;
    for(m = 0; m < Matrix.rows; m++)
    {
        for(n = 0; n < Matrix.columns; n++)
        {
            x = 0;
            for(i = 0; i < Matrix.rows; i++)
            {
                y = 0;
                for(j = 0; j < Matrix.columns; j++)
                {
                    if(i != m && j != n)
                    {
                        Temp.element[x][y] = Matrix.element[i][j];
                        y++;
                    }
                }
                if (i!=m) x++;
            }
            InvertedMatrix->element[n][m] = GM_Pow(-1.0, n+m) * GM_MatDet(Temp) / MatDet;
        }
    }
} /*GM_MatInverse*/

void GM_MatMultiply(GMtype_Matrix MatrixA, GMtype_Matrix MatrixB, GMtype_Matrix *MatrixC)
{
    /*This function multiplies two input matricies A and B, and stores their product in C*/
    int i, j, k;
    GMtype_Matrix TMatrixB;
    GMtype_Data Row, Column;


    GM_MatTranspose(MatrixB, &TMatrixB);
    MatrixC->rows = MatrixA.rows;
    MatrixC->columns = MatrixB.columns;
    for(i = 0; i < MatrixA.rows; i++)
    {
        for(j = 0; j < MatrixB.columns; j++)
        {
            Row.size = MatrixA.columns;
            for(k = 0; k < MatrixA.columns; k++)
            {
                Row.element[k] = MatrixA.element[i][k];
            }
            Column.size = TMatrixB.columns;
            for(k = 0; k < TMatrixB.columns; k++)
            {
                Column.element[k] = TMatrixB.element[j][k];
            }
            MatrixC->element[i][j] = GM_DotProduct(Row, Column);
        }
    }
} /*GM_MatMultiply*/

void GM_MatTranspose(GMtype_Matrix Matrix, GMtype_Matrix *TMatrix)
{
    /*This function takes the transpose of the input matrix*/
    int i, j;
    TMatrix->rows = Matrix.columns;
    TMatrix->columns = Matrix.rows;
    for(i = 0; i < Matrix.rows; i++)
    {
        for(j = 0; j < Matrix.columns; j++)
        {
            TMatrix->element[j][i] = Matrix.element[i][j];
        }
    }
} /*GM_MatTranspose*/

double GM_Mean(GMtype_Data Data)
{
    /*This function average the elements of the input, Data*/
    int i;
    double total = 0, average;
    for(i = 0; i <= Data.size; i++)
    {
        total+=Data.element[i];
    }
    average = total / Data.size;
    return average;
} /*GM_Mean*/

void GM_Median(GMtype_Data Data, double *upper, double *lower)
{
    /*This function finds the median of the input, Data*/
    GM_Sort(&Data);
    if(Data.size%2 == 0) /*If the number of elements is even then there is an upper and lower bound for the median*/
    {
        *upper = Data.element[Data.size/2];
        *lower = Data.element[Data.size/2 - 1];
    }
    else /*If the number of elements is odd then the upper and lower bound are equivalent*/
    {
        *upper = Data.element[Data.size/2];
        *lower = Data.element[Data.size/2];
    }
} /*GM_Meadian*/

void GM_PolyFit(GMtype_Data DataX, GMtype_Data DataY, GMtype_Polynomial *Polynomial)
{
    /*This finds a best polynomial fit of the input degree for the input data set.  Polynomial must have
    **a value in degree, otherwise the function will not know what degree to fit to.*/
    GMtype_Matrix X, TX, XSquared, b, Final, Y, L, U, P;
    int i, j;
    for(i = 0; i < DataX.size; i++)
    {
        for(j = 0; j <= Polynomial->degree; j++)
        {
            X.element[i][j] = GM_Pow(DataX.element[i],j);
        }
    }
    X.rows = DataX.size;
    X.columns = Polynomial->degree + 1;
    for(i = 0; i<= DataY.size; i++)
    {
        Y.element[i][0] = DataY.element[i];
    }
    Y.rows = DataY.size;
    Y.columns = 1;
    GM_MatTranspose(X, &TX);
    GM_MatMultiply(TX, X, &XSquared);
    P.rows = XSquared.rows;
    P.columns = XSquared.columns;
    for(i = 0; i < XSquared.rows; i++)
    {
        for(j = 0; j < XSquared.columns; j++)
        {
            if(i == j) P.element[i][j] = 1;
            else P.element[i][j] = 0;
        }
    }
    GM_LUDecomposition(XSquared, &L, &U, &P);
    GM_MatMultiply(TX, Y, &b);
    GM_LUSolve(L, U, P, &Final, b);
    for(i = 0; i <= Polynomial->degree; i++){
        Polynomial->coef[i] = Final.element[i][0];
    }
} /*GM_PolyFit*/

double GM_Pow(double x, int y)
{
    /*This function raises x to the y power, only works for positive exponents*/
    double value = 1;
    int i;
    if (x == 0)
        return 0;
    if(y<0)
    {
        for(i = 0; i < -1 * y; i++)
        {
            value = value / x;
        }
    }
    else
    {
        for(i = 0; i < y; i++)
        {
            value = value * x;
        }
    }

    return value;
} /*GM_Pow*/

void GM_PrintMatrix(GMtype_Matrix X)
{
    /*This function is used to print a matrix to the console*/
    int i, j;
    printf("\nRows = %d\nColumns = %d\n", X.rows, X.columns);
    for(i = 0; i < X.rows; i++)
    {
        printf("\n");
        for(j = 0; j < X.columns; j++)
        {
            printf("%lf   ", X.element[i][j]);
        }
    }
} /*GM_PrintMatrix*/

double GM_SolvePolynomial(GMtype_Polynomial Polynomial, double x)
{
    /*This function solves a polynomial for a given x*/
    int i;
    double y = 0;
    for(i = 0; i <= Polynomial.degree; i++)
    {
        if(x != 0)
            y += Polynomial.coef[i] * GM_Pow(x, i);
    }
    return y;
} /*GM_SolvePolynomial*/

void GM_Sort(GMtype_Data *Data)
{
    /*This function sorts the array of elements in the input Data*/
    double lowest;
    int i, m, LowestIndex;
    lowest = Data->element[0];
    for(m = 0; m < Data->size; m++)
    {
        lowest = Data->element[m];
        LowestIndex = m;
        for(i = m; i < Data->size; i++)
        {
            if(lowest > Data->element[i])
            {
                lowest = Data->element[i];
                LowestIndex = i;
            }
        }
        GM_Swap(&Data->element[m], &Data->element[LowestIndex]);
    }
} /*GM_Sort*/

double GM_StandardDeviation(GMtype_Data Data)
{
    /*This function takes the standard deviation of the input Data*/
    double stdev, mean, variance;
    int i;
    mean = GM_Mean(Data);
    variance = 0;
    for(i = 0; i < Data.size; i++)
    {
        variance += (Data.element[i] - mean) * (Data.element[i] - mean) / (Data.size-1);
    }
    stdev = sqrt(variance);
    return stdev;
} /*GM_StandardDeviation*/

void GM_Swap(double *x, double *y)
{
    /*This function swaps two numbers*/
    double temp;
    temp = *x;
    *x = *y;
    *y = temp;
} /*GM_Swap*/

void GM_SwapRows(GMtype_Matrix *Matrix, int Row1, int Row2)
{
    int i;
    double a;
    for(i = 0; i < Matrix->columns; i++)
    {
        a = Matrix->element[Row1][i];
        Matrix->element[Row1][i] = Matrix->element[Row2][i];
        Matrix->element[Row2][i] = a;
    }
} /*GM_SwapRows*/

