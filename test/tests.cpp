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
#include <cstdio>	
#include <cmath>
#include <tuple>
int tests_run = 0;

using namespace std;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


int m_equals(Matrix A, Matrix B, double e = 1e-10) {
	if (A.fil != B.fil || A.col != B.col)
		return 0;
	else
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
	B(1,2) = 1.0;
	B(2,1) = 1.0;
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
	TimeDifferences dt = timediff(t1,t2);
	_assert(fabs(dt.UT1_TAI - -9.323958333333489e+03) < 1e-10);
	_assert(fabs(dt.UTC_GPS - -69981) < 1e-10);
	_assert(fabs(dt.UT1_GPS - -9.304958333333489e+03) < 1e-10);
	_assert(fabs(dt.TT_UTC - 7.003218399999999e+04) < 1e-10);
	_assert(fabs(dt.GPS_UTC - 69981) < 1e-10);
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

test_legendre()
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
	Matrix r=eye(3,3);
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
	Matrix r(3);
	Matrix E(3,3);
	r(1)=5;r(2)=5;r(3)=5;
	E(1,1)=1;E(1,2)=4;E(1,3)=7;
	E(2,1)=2;E(2,2)=5;E(2,3)=8;
	E(3,1)=3;E(3,2)=6;E(3,3)=9;
	
	Matrix B = transpose(accelHarmonic(transpose(r),E,1,10))/1e+10;

	Matrix A(3);
	A(1)=-08.39013884498123;
	A(2)=-20.18877159573608;
	A(3)=-31.98740434649093;

	_assert(m_equals(A,B));
	
	return 0;
}

int all_tests()
{

	eop19620101(21413);
	DE430Coeff();
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