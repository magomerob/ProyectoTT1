/**
 * @file tests.cpp
 * @author Marcos Gómez Robres
 * @brief Archivo que contiene todos los test
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * 
 */
#define _USE_MATH_DEFINES
#include "../include/JPL_Eph.h"
#include "../include/matrix.h"
#include "../include/cheb3d.h"
#include "../include/eccAnom.h"
#include "../include/accelPointMass.h"
#include "../include/timediff.h"
#include "../include/meanObliquity.h"
#include "../include/frac.h"
#include "../include/mjday.h"
#include "../include/mjday_tdb.h"
#include "../include/position.h"
#include "../include/sign.h"
#include "../include/r_x.h"
#include "../include/r_y.h"
#include "../include/r_z.h"
#include "../include/timediff.h"
#include "../include/azelpa.h"
#include "../include/iers.h"
#include "../include/legendre.h"
#include "../include/timeUpdate.h"
#include "../include/gmst.h"
#include "../include/nutAngles.h"
#include "../include/precMatrix.h"
#include "../include/global.h"
#include "../include/poleMatrix.h"
#include "../include/nutMatrix.h"
#include "../include/LTC.h"
#include "../include/eqnEquinox.h"
#include "../include/accelHarmonic.h"
#include "../include/g_accelHarmonic.h"
#include "../include/varEqn.h"
#include "../include/ghaMatrix.h"
#include "../include/gast.h"
#include "../include/accel.h"
#include "../include/measUpdate.h"
#include "../include/DEInteg.h"
#include <cstdio>	
#include <cmath>
#include <tuple>
int tests_run = 0;

using namespace std;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


int m_equals(Matrix A, Matrix B, double e = 1e-10) {
	if (A.fil != B.fil || A.col != B.col){
		printf("Matrix sizes do not match: %i %i != %i %i\n", A.fil, A.col, B.fil, B.col);
		return 0;
	}else
		for(int i = 1; i <= A.fil; i++)
			for(int j = 1; j <= A.col; j++)
				if(fabs(A(i,j)-B(i,j)) > e) {
					printf("%2.20lf %2.20lf pos %i %i\n",A(i,j),B(i,j),i,j);
					return 0;
				}
	
	return 1;
}

int test_sum()
{
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = 5.0;
	B(1,2) = 6.0;
	B(2,1) = 7.0;
	B(2,2) = 8.0;

	C = A + B;

	Matrix D(2,2);
	D(1,1) = 6.0;
	D(1,2) = 8.0;
	D(2,1) = 10.0;
	D(2,2) = 12.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_sub()
{
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = 5.0;
	B(1,2) = 6.0;
	B(2,1) = 7.0;
	B(2,2) = 8.0;

	C = A - B;

	Matrix D(2,2);
	D(1,1) = -4.0;
	D(1,2) = -4.0;
	D(2,1) = -4.0;	
	D(2,2) = -4.0;
	_assert(m_equals(C,D));

	return 0;
}

int test_mult()
{
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = 5.0;
	B(1,2) = 6.0;
	B(2,1) = 7.0;
	B(2,2) = 8.0;

	C = A * B;

	Matrix D(2,2);
	D(1,1) = 19.0;
	D(1,2) = 22.0;
	D(2,1) = 43.0;	
	D(2,2) = 50.0;
	_assert(m_equals(C,D));

	return 0;
}

int test_div()
{
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = 5.0;
	B(1,2) = 6.0;
	B(2,1) = 7.0;
	B(2,2) = 8.0;

	C = A / B;

	Matrix D(2,2);
	D(1,1) = 3;
	D(1,2) = -2;
	D(2,1) = 2;	
	D(2,2) = -1;

	_assert(m_equals(C,D));

	return 0;
}

int test_set()
{
	Matrix A(2,2);
	Matrix B(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 1.0;
	A(2,1) = 1.0;
	A(2,2) = 1.0;

	B(1,1) = 1.0;
	B(1,2) = 2.0;
	B(2,1) = 3.0;
	B(2,2) = 4.0;

	A = B;

	_assert(m_equals(A,B));

	return 0;
}

int test_zeros()
{
	Matrix A = zeros(2,2);
	Matrix B(2,2);

	B(1,1) = 0.0;
	B(1,2) = 0.0;
	B(2,1) = 0.0;
	B(2,2) = 0.0;

	_assert(m_equals(A,B));

	return 0;
}

int test_eye()
{
	Matrix A = eye(2,2);
	Matrix B(2,2);

	B(1,1) = 1.0;
	B(1,2) = 0.0;
	B(2,1) = 0.0;
	B(2,2) = 1.0;

	_assert(m_equals(A,B));

	return 0;
}

int test_transpose()
{
	Matrix A(2,2);
	Matrix B(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = 1.0;
	B(1,2) = 3.0;
	B(2,1) = 2.0;
	B(2,2) = 4.0;

	A = transpose(A);

	_assert(m_equals(A,B));

	return 0;
}

int test_inverse()
{
	Matrix A(2,2);
	Matrix B(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1,1) = -2.0;
	B(1,2) = 1.0;
	B(2,1) = 1.5;
	B(2,2) = -0.5;

	A = inv(A);

	_assert(m_equals(A,B));

	return 0;
}

int test_sum_const()
{
	Matrix A(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	C = A + 1.0;

	Matrix D(2,2);
	D(1,1) = 2.0;
	D(1,2) = 3.0;
	D(2,1) = 4.0;
	D(2,2) = 5.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_sub_const()
{
	Matrix A(2,2);
	Matrix B(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	C = A - 1.0;

	Matrix D(2,2);
	D(1,1) = 0.0;
	D(1,2) = 1.0;
	D(2,1) = 2.0;	
	D(2,2) = 3.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_mult_const()
{
	Matrix A(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	C = A * 2.0;

	Matrix D(2,2);
	D(1,1) = 2.0;
	D(1,2) = 4.0;
	D(2,1) = 6.0;	
	D(2,2) = 8.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_div_const()
{
	Matrix A(2,2);
	Matrix C(2,2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	C = A / 2.0;

	Matrix D(2,2);
	D(1,1) = 0.5;
	D(1,2) = 1.0;
	D(2,1) = 1.5;	
	D(2,2) = 2.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_vector(){
	Matrix A(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;

	_assert(A(2) == 2.0);

	return 0;
}

int test_zeros_vector(){
	Matrix A = zeros(3);
	Matrix B(3);

	B(1) = 0.0;
	B(2) = 0.0;
	B(3) = 0.0;

	_assert(m_equals(A,B));

	return 0;
}

int test_norm(){
	Matrix A(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;

	double vnorm = norm(A);

	_assert(fabs(vnorm - sqrt(14)) < 1e-10);

	return 0;
}

int test_dot(){
	Matrix A(3);
	Matrix B(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;

	B(1) = 4.0;
	B(2) = 5.0;
	B(3) = 6.0;

	double vdot = dot(A,B);

	_assert(fabs(vdot - 32) < 1e-10);

	return 0;
}

int test_cross(){
	Matrix A(3);
	Matrix B(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;

	B(1) = 4.0;
	B(2) = 5.0;
	B(3) = 6.0;

	Matrix C = cross(A,B);

	Matrix D(3);
	D(1) = -3.0;
	D(2) = 6.0;
	D(3) = -3.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_extract_vector()
{
	Matrix A(6);
	Matrix B(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;
	A(4) = 4.0;
	A(5) = 5.0;
	A(6) = 6.0;

	B = extract_vector(A, 2, 5);
	
	Matrix D(4);
	D(1) = 2.0;
	D(2) = 3.0;
	D(3) = 4.0;
	D(4) = 5.0;

	_assert(m_equals(B,D));

	return 0;
}

int test_union_vector()
{
	Matrix A(3);
	Matrix B(3);
	
	A(1) = 1.0;
	A(2) = 2.0;
	A(3) = 3.0;

	B(1) = 4.0;
	B(2) = 5.0;
	B(3) = 6.0;

	Matrix C = union_vector(A,B);

	Matrix D(6);
	D(1) = 1.0;
	D(2) = 2.0;
	D(3) = 3.0;
	D(4) = 4.0;
	D(5) = 5.0;
	D(6) = 6.0;

	_assert(m_equals(C,D));

	return 0;
}

int test_extract_row()
{
	Matrix A(2,2);
	Matrix B(2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B = A.extract_row(1);

	Matrix D(2);
	D(1) = 1.0;
	D(2) = 2.0;

	_assert(m_equals(B,D));

	return 0;
}

int test_extract_column()
{
	Matrix A(2,2);
	Matrix B(2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B = A.extract_column(1);

	Matrix D(2);
	D(1) = 1.0;
	D(2) = 3.0;

	_assert(m_equals(B,D));

	return 0;
}

int test_assign_row()
{
	Matrix A(2,2);
	Matrix B(2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1) = 5.0;
	B(2) = 6.0;

	A.assign_row(B,1);

	Matrix D(2,2);
	D(1,1) = 5.0;
	D(1,2) = 6.0;
	D(2,1) = 3.0;	
	D(2,2) = 4.0;

	_assert(m_equals(A,D));

	return 0;
}

int test_assign_column()
{
	Matrix A(2,2);
	Matrix B(2);
	
	A(1,1) = 1.0;
	A(1,2) = 2.0;
	A(2,1) = 3.0;
	A(2,2) = 4.0;

	B(1) = 5.0;
	B(2) = 6.0;

	A.assign_column(B,1);

	Matrix D(2,2);
	D(1,1) = 5.0;
	D(1,2) = 2.0;
	D(2,1) = 6.0;	
	D(2,2) = 4.0;

	_assert(m_equals(A,D));

	return 0;
}

int test_accelPointMass()
{
	Matrix u(3);
	Matrix v(3);
	u(1)=0;u(2)=1;u(3)=2;
	v(1)=4;v(2)=3;v(3)=2;
	double GM = 2;
	
	Matrix R(3);
	R(1)=0.0382164189132187;
	R(2)=0.00630163440991609;
	R(3)=-0.0256131500933865;
	
	Matrix B = AccelPointMass(u,v,GM);

	_assert(m_equals(R,B));
	
	return 0;
}

int test_cheb3d()
{
	Matrix Cx(5);
	Matrix Cy(5);
	Matrix Cz(5);
	Cx(1)=1;Cx(2)=2;Cx(3)=3;Cx(4)=4;Cx(5)=5;
	Cy(1)=-1;Cy(2)=-2;Cy(3)=-3;Cy(4)=-4;Cy(5)=-5;
	Cz(1)=0;Cz(2)=1;Cz(3)=2;Cz(4)=3;Cz(5)=4;
	double t = 2;
	double N = 5;
	double Ta = 1;
	double Tb = 3;
	
	Matrix R(3);
	R(1)=3;R(2)=-3;R(3)=2;
	Matrix A = Cheb3D(t,N,Ta,Tb,Cx,Cy,Cz);
	_assert(m_equals(R,A, 1e-10));
	
	return 0;
}

int test_eccAnom()
{
	double M = 2;
	double e = 0.16;
	double E = eccAnom(M,e);
	_assert(fabs(E - 2.13518640114338) < 1e-10);
	return 0;
}

int test_frac(){
	double x=5.2352345;
	
	_assert(abs(0.2352345 - frac(x)) < 1e-10);
	
	return 0;
}

int test_mean_obliquity()
{
	double Mjd_TT=13;
		
	_assert(abs(0.409412989421167 - MOblq(Mjd_TT)) < 1e-10);
	
	return 0;
}

int test_mjd()
{
	_assert(abs(52561.08333333349 - Mjday(2002,10,14,2,00,0)) < 1e-10);
	
	return 0;
}

int test_mjday_tdb()
{
	_assert(abs(58340.4456185067 - Mjd_TDB(58340.445618518)) < 1e-10);
	
	return 0;
}

int test_position(){
	
	Matrix R(3);
	R(1)=-264748.079560846;
	R(2)=4802726.74400269;
	R(3)=-4175137.27246505;
	
	Matrix B = Position(42.46656962435925, -2.423409704452636, 384);
	
	_assert(m_equals(R,B,1e-6));
	
	return 0;
}

int test_sign()
{
	
	_assert(abs(-5-sign(5,-1))<1e-10);
	return 0;

}

int test_rx()
{
	Matrix A(3,3);
	Matrix B(3,3);
	A(1,1) = 1.0;A(1,2) = 0;A(1,3) = 0;
	A(2,1) = 0;A(2,2) = 0;A(2,3) = 1.0;
	A(3,1) = 0;A(3,2) = -1.0;A(3,3) = 0;

	B = r_x(M_PI/2);

	_assert(m_equals(A,B));

	return 0;
}

int test_ry()
{
	Matrix A(3,3);
	Matrix B(3,3);
	A(1,1) = 0;A(1,2) = 0;A(1,3) = -1.0;
	A(2,1) = 0;A(2,2) = 1.0;A(2,3) = 0;
	A(3,1) = 1.0;A(3,2) = 0;A(3,3) = 0;

	B = r_y(M_PI/2);

	_assert(m_equals(A,B));

	return 0;
}

int test_rz()
{
	Matrix A(3,3);
	Matrix B(3,3);
	A(1,1) = 0;A(1,2) = 1.0;A(1,3) = 0;
	A(2,1) = -1.0;A(2,2) = 0;A(2,3) = 0;
	A(3,1) = 0;A(3,2) = 0;A(3,3) = 1.0;

	B = r_z(M_PI/2);

	_assert(m_equals(A,B));

	return 0;
}

int test_timediff()
{
	double t1 = 6.067604166666651e04;
	double t2 = 7.0e04;
	auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(t1,t2);
	_assert(fabs(UT1_TAI - -9.323958333333489e+03) < 1e-10);
	_assert(fabs(UTC_GPS - -69981) < 1e-10);
	_assert(fabs(UT1_GPS - -9.304958333333489e+03) < 1e-10);
	_assert(fabs(TT_UTC - 7.003218399999999e+04) < 1e-10);
	_assert(fabs(GPS_UTC - 69981) < 1e-10);
	return 0;
}

int test_azelpa()
{
	Matrix r(3);
	r(1)=1; r(2)=2; r(3)=3;
	double a = 0.463647609000806;
	double b = 0.930274014115472;
	tuple<double, double, Matrix, Matrix> result = azelpa(r);

	_assert(fabs(get<0>(result) - a) < 1e-10);
	_assert(fabs(get<1>(result) - b) < 1e-10);
	return 0;
}

int test_iers()
{
	tuple<double, double, double, double, double, double, double, double, double> result = iers(eopdata, 40340, 'l');
	
	_assert(fabs(get<0>(result) - 1.24267442741996e-07) < 1e-10);
	_assert(fabs(get<1>(result) - 1.72124370831681e-06) < 1e-10);
	_assert(fabs(get<2>(result) - -0.0025048) < 1e-10);
	_assert(fabs(get<3>(result) -  0.003476) < 1e-10);
	_assert(fabs(get<4>(result) - 2.81565241577985e-07) < 1e-10);
	_assert(fabs(get<5>(result) - 1.7133315490411e-08) < 1e-10);
	_assert(fabs(get<6>(result) -  0) < 1e-10);
	_assert(fabs(get<7>(result) - 0) < 1e-10);
	_assert(fabs(get<8>(result) - 7) < 1e-10);
	return 0;
}

int test_legendre()
{
	auto [A,B] = legendre(2, 2, 1);

	Matrix C(3,3);
	Matrix D(3,3);
    C(1,1)=1.00000000000000000000;C(1,2)=0.00000000000000000000;C(1,3)=0.00000000000000000000;
    C(2,1)=1.45747045027529753547;C(2,2)=0.93583099453595375294;C(2,3)=0.00000000000000000000;
    C(3,1)=1.25691629764617052167;C(3,2)=1.76084674147066455596;C(3,3)=0.56531333344858780698;

    D(1,1)=0.00000000000000000000;D(1,2)=0.00000000000000000000;D(1,3)=0.00000000000000000000;
    D(2,1)=0.93583099453595375294;D(2,2)=-1.45747045027529753547;D(2,3)=0.00000000000000000000;
    D(3,1)=3.04987602056929052452;D(3,2)=-1.61172970742831900282;D(3,3)=-1.76084674147066455596;

	_assert(m_equals(A,C));
	_assert(m_equals(B,D));
	return 0;
}

int test_nutangles()
{
	tuple<double, double> result= NutAngles(34000);
	_assert(fabs(3.51485038567364e-05-get<0>(result))<1e-10);
	_assert(fabs(3.75448367402296e-05-get<1>(result))<1e-10);
	return 0;
}

int test_timeupdate()
{
	Matrix r=zeros(3,3);
	r(1,1)=1;r(1,2)=1;r(1,3)=1;
	r(2,1)=1;r(2,2)=1;r(2,3)=1;
	r(3,1)=1;r(3,2)=1;r(3,3)=1;
	Matrix v(3);
	v(1)=1; v(2)=2; v(3)=3;
	Matrix result = TimeUpdate(r, v, 1.0);
	Matrix expected(1);
	expected(1)=37;

	_assert(m_equals(result, transpose(expected), 1e-10));
	return 0;
}

int test_gmst()
{
	_assert(abs(3.86479-gmst(2460803))<1e-5);
	return 0;

}

int test_precMatrix()
{
	Matrix m = precMatrix(2450803,2460803);
	Matrix m2 = zeros(3,3);

	
	m2(1,1)=0.999976404960704;m2(1,2)=-0.00636469924721845;m2(1,3)=-0.00258459384809681;
	m2(2,1)=0.00636469924622329;m2(2,2)=0.999979745062791;m2(2,3)=-8.22556332843431e-06;
	m2(3,1)=0.00258459385054743;m2(3,2)=-8.22479327083077e-06;m2(3,3)=0.999996659897912;
	
	_assert(m_equals(m,m2));
	return 0;
}

int test_jpl()
{
	auto [a, b, c, d, e, f, g, h, i, j, k] = JPL_Eph(50000);
	Matrix A(3);
	A(1)=-106013268114.407;
	A(2)=-10704774110.6778;
	A(3)=-5733414115.03038;
	_assert(m_equals(a,A, 1e-2));
	Matrix B(3);
	B(1)=-215428191214.142;
	B(2)=-113786449444.223;
	B(3)=-46005254439.599;
	_assert(m_equals(b,B, 1e-2));
	Matrix C(3);
	C(1)=142969439950.568;
	C(2)=39449684228.9347;
	C(3)=17113648323.9793;
	_assert(m_equals(c,C, 1e-2));
	Matrix D(3);
	D(1)=-198514069852.68;
	D(2)=-234901827954.399;
	D(3)=-105286846732.194;
	_assert(m_equals(d,D, 1e-2));
	Matrix E(3);
	E(1)=-268216727765.582;
	E(2)=-758862052731.16;
	E(3)=-322431019777.682;
	_assert(m_equals(e,E, 1e-2));
	Matrix F(3);
	F(1)=1278811202046.14;
	F(2)=-195324072722.162;
	F(3)=-142638833473.58;
	_assert(m_equals(f,F, 1e-2));
	Matrix G(3);
	G(1)=1307661864082.46;
	G(2)=-2385778983942.11;
	G(3)=-1065262865019.22;
	_assert(m_equals(g,G, 1e-2));
	Matrix H(3);
	H(1)=1743981369482.03;
	H(2)=-3815582390654.74;
	H(3)=-1609685101213.35;
	_assert(m_equals(h,H, 1e-2));
	Matrix I(3);
	I(1)=-2300231562026.75;
	I(2)=-3910278366850.89;
	I(3)=-575107437973.647;
	_assert(m_equals(i,I, 1e-2));
	Matrix J(3);
	J(1)=332092361.309373;
	J(2)=192658201.912627;
	J(3)=79678751.0600327;
	_assert(m_equals(j,J, 1e-2));
	Matrix K(3);
	K(1)=-143417093159.564;
	K(2)=-38421114445.7267;
	K(3)=-16658394353.0936;
	_assert(m_equals(k,K, 1e-2));

	return 0;
}

int test_poleMatrix()
{
	Matrix A(3,3);
	Matrix B(3,3);
	A(1,1) = 0;A(1,2) = 1.0;A(1,3) = 0;
	A(2,1) = 0;A(2,2) = 0;A(2,3) = -1.0;
	A(3,1) = -1.0;A(3,2) = 0;A(3,3) = 0;

	B = poleMatrix(M_PI/2,M_PI/2);

	_assert(m_equals(A,B));

	return 0;
}

int test_nutMatrix()
{
	Matrix A  = nutMatrix(50000);
	Matrix B(3,3);
	B(1,1)=0.999999999436159;B(1,2)=-3.08098181393545e-05;B(1,3)=-1.3358042226427e-05;
	B(2,1)=3.08103201452149e-05;B(2,2)=0.99999999881915;B(2,3)=3.75822191258646e-05;
	B(3,1)=1.33568843093166e-05;B(3,2)=-3.75826306702187e-05;B(3,3)=0.99999999920457;	

	_assert(m_equals(A,B,1e-10));
	return 0;
}

int test_LTC()
{
	Matrix A(3,3);
	Matrix B(3,3);
	A(1,1) = -1.22464679914735e-16;A(1,2) = -1;A(1,3) = 0;
	A(2,1) = 1.22464679914735e-16;A(2,2) = -1.49975978266186e-32;A(2,3) = -1;
	A(3,1) = 1.0;A(3,2) = -1.22464679914735e-16;A(3,3) = 1.22464679914735e-16;

	B = LTC(M_PI,M_PI);

	_assert(m_equals(A,B));

	return 0;
}

int test_eqnEquinox()
{
	double a = eqnEquinox(50000);
	_assert(fabs(a - 3.08098181451451e-05) < 1e-10);

	return 0;
}

int test_accelHarmonic()
{
	double a = -4.783104375562398e+09;
	Matrix r(3,3);
	r(1,1)=1; r(1,2)=2; r(1,3)=3;
	r(2,1)=4; r(2,2)=5; r(2,3)=6;
	r(3,1)=7; r(3,2)=8; r(3,3)=9;
	Matrix E(3);
	E(1)=10; E(2)=10; E(3)=10;
	Matrix expected(3);
	expected(1)=-3.851815228199548e+10; expected(2)=-4.592548925930230e+10; expected(3)=-5.333282623660913e+10;
	Matrix result = accelHarmonic(transpose(E),r, 1, 10);
	_assert(m_equals(transpose(expected), result, 1e-10));
	return 0;
}

int test_g_accelHarmonic()
{
	Matrix r(3);
	Matrix E(3,3);
	r(1)=5;r(2)=5;r(3)=5;
	E(1,1)=1;E(1,2)=4;E(1,3)=7;
	E(2,1)=2;E(2,2)=5;E(2,3)=8;
	E(3,1)=3;E(3,2)=6;E(3,3)=9;
	
	Matrix B = g_accelHarmonic(transpose(r),E,1,10);

	Matrix A(3,3);
	A(1,1)=4528002282.97842;A(1,2)=11211752958.2118;A(1,3)=17948017765.0721;
	A(2,1)=11190806609.3086;A(2,2)=26978280555.6971;A(2,3)=42887485497.8079;
	A(3,1)=17853610935.6388;A(3,2)=42744808153.1825;A(3,3)=67826953230.5435;

	_assert(m_equals(A,B, 1e-4));
	
	return 0;
}

int test_gast()
{
	double Mjd_TT=50000;
	double exp = 0.316428666579827;
	_assert(abs(exp - gast(Mjd_TT)) < 1e-10);
	return 0;
}

int test_ghaMatrix()
{
	double Mjd_TT=50000;
	Matrix exp(3,3);
	exp(1,1)=0.950352784296835;exp(1,2)=0.311174525594968;exp(1,3)=0;
	exp(2,1)=-0.311174525594968;exp(2,2)=0.950352784296835;exp(2,3)=0;
	exp(3,1)=0;exp(3,2)=0;exp(3,3)=1;
	
	_assert(m_equals(exp, GHAMatrix(Mjd_TT), 1e-10));

	return 0;
}

int test_varEqn()
{
	double x = 5.38970808087706;
	Matrix yPhi(42, 1);
	yPhi(1,1)=7101800.90695315;
	yPhi(2,1)=1293997.58115302;
	yPhi(3,1)=10114.014948955;
	yPhi(4,1)=573.068082065557;
	yPhi(5,1)=-3085.15736953138;
	yPhi(6,1)=-6736.03068347156;
	yPhi(7,1)=1.0000293469741;
	yPhi(8,1)=8.22733917593032e-06;
	yPhi(9,1)=2.17104932968693e-07;
	yPhi(10,1)=1.08925458231315e-05;
	yPhi(11,1)=3.04673932160225e-06;
	yPhi(12,1)=6.63504292706821e-08;
	yPhi(13,1)=8.22733944423959e-06;
	yPhi(14,1)=0.999986101965304;
	yPhi(15,1)=3.99927483270551e-08;
	yPhi(16,1)=3.04673960163327e-06;
	yPhi(17,1)=-5.1596062466179e-06;
	yPhi(18,1)=1.22075292404534e-08;
	yPhi(19,1)=2.17105640392839e-07;
	yPhi(20,1)=3.9992870847826e-08;
	yPhi(21,1)=0.999984551298692;
	yPhi(22,1)=6.63510875632706e-08;
	yPhi(23,1)=1.22076480274715e-08;
	yPhi(24,1)=-5.73276287738792e-06;
	yPhi(25,1)=5.38976081674752;
	yPhi(26,1)=1.47507305174403e-05;
	yPhi(27,1)=3.21241787851554e-07;
	yPhi(28,1)=1.00002936035846;
	yPhi(29,1)=8.19365458482084e-06;
	yPhi(30,1)=1.40504658112974e-07;
	yPhi(31,1)=1.47507306419397e-05;
	yPhi(32,1)=5.38968310056198;
	yPhi(33,1)=5.90697768748029e-08;
	yPhi(34,1)=8.19365482653896e-06;
	yPhi(35,1)=0.9999860891763;
	yPhi(36,1)=2.58022974647481e-08;
	yPhi(37,1)=3.21242427100724e-07;
	yPhi(38,1)=5.90698876854246e-08;
	yPhi(39,1)=5.38968032557769;
	yPhi(40,1)=1.4050537070756e-07;
	yPhi(41,1)=2.58024285760964e-08;
	yPhi(42,1)=0.999984550703337;

	Matrix expected(42, 1);
	expected(1,1)=573.068082065557;
	expected(2,1)=-3085.15736953138;
	expected(3,1)=-6736.03068347156;
	expected(4,1)=-7.53489822593659;
	expected(5,1)=-1.37294429126638;
	expected(6,1)=-0.0107597986473575;
	expected(7,1)=1.08925458231315e-05;
	expected(8,1)=3.04673932160225e-06;
	expected(9,1)=6.63504292706821e-08;
	expected(10,1)=2.02239897508587e-06;
	expected(11,1)=5.61811901849645e-07;
	expected(12,1)=4.39846387071934e-09;
	expected(13,1)=3.04673960163327e-06;
	expected(14,1)=-5.1596062466179e-06;
	expected(15,1)=1.22075292404534e-08;
	expected(16,1)=5.61812134084449e-07;
	expected(17,1)=-9.58613689243416e-07;
	expected(18,1)=8.05616500343474e-10;
	expected(19,1)=6.63510875632706e-08;
	expected(20,1)=1.22076480274715e-08;
	expected(21,1)=-5.73276287738792e-06;
	expected(22,1)=4.39895597958216e-09;
	expected(23,1)=8.0570607835305e-10;
	expected(24,1)=-1.06368693580442e-06;
	expected(25,1)=1.00002936035846;
	expected(26,1)=8.19365458482084e-06;
	expected(27,1)=1.40504658112974e-07;
	expected(28,1)=1.08999102436198e-05;
	expected(29,1)=3.02797128053784e-06;
	expected(30,1)=2.37068516291712e-08;
	expected(31,1)=8.19365482653896e-06;
	expected(32,1)=0.9999860891763;
	expected(33,1)=2.58022974647481e-08;
	expected(34,1)=3.02797160153579e-06;
	expected(35,1)=-5.16671243316801e-06;
	expected(36,1)=4.34211426867344e-09;
	expected(37,1)=1.4050537070756e-07;
	expected(38,1)=2.58024285760964e-08;
	expected(39,1)=0.999984550703337;
	expected(40,1)=2.37075280907946e-08;
	expected(41,1)=4.34223837651307e-09;
	expected(42,1)=-5.73302112206999e-06;

	Matrix result = varEqn(x, yPhi);

	_assert(m_equals(result,expected,1e-8));
	return 0;
}

int test_accel()
{
	double x =-543.476874884521;
	Matrix y(6);

	y(1)=5720694.2260585;
	y(2)=2687728.41425142;
	y(3)=3483000.08675422;
	y(4)=4371.83136151615;
	y(5)=-1905.47309296258;
	y(6)=-5698.58341612187;

	Matrix acc = accel(x,transpose(y));

	Matrix exp(6);

	exp(1)=4371.83136151615;
	exp(2)=-1905.47309296258;
	exp(3)=-5698.58341612187;
	exp(4)=-6.0654420426171;
	exp(5)=-2.84977703178268;
	exp(6)=-3.70232534578346;

	_assert(m_equals(acc,transpose(exp), 1e-6));

	return 0;
}

int test_measUpdate()
{
	Matrix ix(6);
	ix(1)=7081709.15235124;
	ix(2)=1384037.58387016;
	ix(3)=208036.562137862;
	ix(4)=794.215435663787;
	ix(5)=-3043.46028560081;
	ix(6)=-6732.63769538921;

	double z = 0.166638800992663;

	double g = 0.166498010440566;

	double s =  0.000242600766027212;

	Matrix G(6);
	G(1)=3.38224487276275e-07;
	G(2)=6.92583388658317e-08;
	G(3)=2.08370308916034e-07;
	G(4)=0;
	G(5)=0;
	G(6)=0;

	Matrix iP(6,6);
	iP(1,1)=16158.3302232928;	iP(1,2)=-6010.47704001827;	iP(1,3)=8909.2403033278;	iP(1,4)=52.3928185697041;	iP(1,5)=-15.4569900201582;	iP(1,6)=25.529764531455;
	iP(2,1)=-6010.47704001827;	iP(2,2)=22876.1753800742;	iP(2,3)=-778.334152884529;	iP(2,4)=-6.03444081531537;	iP(2,5)=60.3316340659165;	iP(2,6)=-26.6862579576286;
	iP(3,1)=8909.2403033278;	iP(3,2)=-778.334152884511;	iP(3,3)=7060.63817759347;	iP(3,4)=30.1600339164459;	iP(3,5)=-2.79055554931059;	iP(3,6)=18.8555622209732;
	iP(4,1)=52.3928185697041;	iP(4,2)=-6.03444081531539;	iP(4,3)=30.1600339164459;	iP(4,4)=0.205777460607875;	iP(4,5)=-0.0250235926760925;iP(4,6)=0.0818687411589902;
	iP(5,1)=-15.4569900201582;	iP(5,2)=60.3316340659165;	iP(5,3)=-2.79055554931064;	iP(5,4)=-0.0250235926760924;iP(5,5)=0.168178931856111;	iP(5,6)=-0.0825059738208085;
	iP(6,1)=25.529764531455;	iP(6,2)=-26.6862579576284;	iP(6,3)=18.8555622209731;	iP(6,4)=0.0818687411589899;	iP(6,5)=-0.0825059738208081;iP(6,6)=0.102560261447662;

	int n = 6;
	auto [K,x,P] = measUpdate(transpose(ix), z, g, s, G, iP, n);

	Matrix expK(6);
	expK(1)=111247.169980744;
	expK(2)=-9838.74201513926;
	expK(3)=71379.5898347798;
	expK(4)=379.997602081885;
	expK(5)=-26.2750278118874;
	expK(6)=172.631075925264;

	_assert(m_equals(K,transpose(expK), 1e-6));

	Matrix expX(6);
	expX(1)=7081724.81490172;
	expX(2)=1384036.19866824;
	expX(3)=208046.611709723;
	expX(4)=794.268935735979;
	expX(5)=-3043.46398487648;
	expX(6)=-6732.61339056473;

	_assert(m_equals(x,transpose(expX), 1e-6));

	Matrix expP(6,6);
	expP(1,1)=15390.1364249244;	expP(1,2)=-5942.53768838276;expP(1,3)=8416.34364581595;	expP(1,4)=49.7688255156362;		expP(1,5)=-15.2755533739896;	expP(1,6)=24.3376971561284;
	expP(2,1)=-5942.53768838276;expP(2,2)=22870.1667979814;	expP(2,3)=-734.742185068871;expP(2,4)=-5.80237387183302;	expP(2,5)=60.3155877399376;		expP(2,6)=-26.5808310656678;
	expP(3,1)=8416.34364581594;	expP(3,2)=-734.742185068853;expP(3,3)=6744.38059324125;	expP(3,4)=28.4763996454698;		expP(3,5)=-2.67414024286502;	expP(3,6)=18.0906953005745;
	expP(4,1)=49.7688255156362;	expP(4,2)=-5.80237387183305;expP(4,3)=28.4763996454698;	expP(4,4)=0.196814436482425;	expP(4,5)=-0.0244038421648142;	expP(4,6)=0.0777968825739801;
	expP(5,1)=-15.2755533739896;expP(5,2)=60.3155877399376;	expP(5,3)=-2.67414024286508;expP(5,4)=-0.0244038421648141;	expP(5,5)=0.168136079054327;	expP(5,6)=-0.0822244241558079;
	expP(6,1)=24.3376971561283;	expP(6,2)=-26.5808310656676;expP(6,3)=18.0906953005745;	expP(6,4)=0.0777968825739798;	expP(6,5)=-0.0822244241558074;	expP(6,6)=0.100710435752198;
	
	_assert(m_equals(P,expP));

	return 0;
}

int test_DEInteg()
{
	AuxParam.Mjd_UTC=4.974611128472211e+04;
	Matrix Y0_apr(6, 1);
	Y0_apr(1, 1) = 6.221397628578685e+06;
	Y0_apr(2, 1) = 2.867713779657379e+06;
	Y0_apr(3, 1) = 3.006155985099489e+06;
	Y0_apr(4, 1) = 4.645047251618060e+03;
	Y0_apr(5, 1) = -2.752215915882042e+03;
	Y0_apr(6, 1) = -7.507999409870306e+03;
	Matrix result = DEInteg(accel, 0.0, -1.349999919533730e+02, 1e-13, 1e-6, 6, Y0_apr);

	Matrix expected(6, 1);
	expected(1, 1) = 5.542555937228607e+06;
	expected(2, 1) = 3.213514867349196e+06;
	expected(3, 1) = 3.990892975876853e+06;
	expected(4, 1) = 5.394068421663513e+03;
	expected(5, 1) = -2.365213378823415e+03;
	expected(6, 1) = -7.061845542002954e+03;
	// cout << result << endl;
	_assert(m_equals(result, expected, 1e-6));
	return 0;
}

int all_tests()
{

	eop19620101(21413);
    DE430Coeff();
    auxparam();
    GGM03S();
    GEOS3();

	_verify(test_sum);
	_verify(test_sub);
	_verify(test_mult);
	_verify(test_div);
	_verify(test_set);
	_verify(test_zeros);
	_verify(test_eye);
	_verify(test_transpose);
	_verify(test_inverse);
	_verify(test_sum_const);
	_verify(test_sub_const);
	_verify(test_mult_const);
	_verify(test_div_const);
	_verify(test_vector);
	_verify(test_zeros_vector);
	_verify(test_norm);
	_verify(test_dot);
	_verify(test_cross);
	_verify(test_extract_vector);
	_verify(test_union_vector);
	_verify(test_extract_row);
	_verify(test_extract_column);
	_verify(test_assign_row);
	_verify(test_assign_column);
	_verify(test_accelPointMass);
	_verify(test_cheb3d);
	_verify(test_eccAnom);
	_verify(test_frac);
	_verify(test_mean_obliquity);
	_verify(test_mjd);
	_verify(test_mjday_tdb);
	_verify(test_position);
	_verify(test_sign);
	_verify(test_rx);
	_verify(test_ry);
	_verify(test_rz);
	_verify(test_timediff);
	_verify(test_azelpa);
	_verify(test_iers);
	_verify(test_legendre);
	_verify(test_nutangles);
	_verify(test_timeupdate);
	_verify(test_gmst);
	_verify(test_precMatrix);
	_verify(test_jpl);
	_verify(test_poleMatrix);
	_verify(test_nutMatrix);
	_verify(test_LTC);
	_verify(test_eqnEquinox);
	_verify(test_accelHarmonic);
	_verify(test_g_accelHarmonic);
	_verify(test_gast);
	_verify(test_ghaMatrix);
	_verify(test_varEqn);
	_verify(test_accel);
	_verify(test_measUpdate);
	_verify(test_DEInteg);
	return 0;
}

int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);
    return 0;
}