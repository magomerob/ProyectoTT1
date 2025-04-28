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
#include <cstdio>
#include <cmath>
#include <tuple>
int tests_run = 0;

using namespace std;

#define FAIL() printf("/nfailure in %s() line %d/n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


int m_equals(Matrix A, Matrix B, double e = 1e-10) {
	if (A.fil != B.fil || A.col != B.col)
		return 0;
	else
		for(int i = 1; i <= A.fil; i++)
			for(int j = 1; j <= A.col; j++)
				if(fabs(A(i,j)-B(i,j)) > e) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
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


int all_tests()
{
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