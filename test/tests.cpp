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

#include <cstdio>
#include <cmath>
#include <tuple>
int tests_run = 0;

using namespace std;

#define FAIL() printf("/nfailure in %s() line %d/n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)


int all_tests()
{
	/*
    _verify(m_sum_01);
    _verify(m_sub_01);
	_verify(m_mul_01);
	_verify(m_div_01);
	_verify(m_asig_01);
    _verify(m_zeros_01);
	_verify(v_zeros_01);
	_verify(m_eye_01);
	_verify(m_transpose_01);
	_verify(m_inv_01);
	_verify(m_determinante_01);
	_verify(m_submatriz_01);
	_verify(m_sumdouble_01);
	_verify(m_subdouble_01);
	_verify(m_muldouble_01);
	_verify(m_divdouble_01);
	_verify(m_extractvec_01);
	_verify(m_unionvec_01);
	_verify(m_extractrow_01);
	_verify(m_extractcol_01);
	_verify(m_assignrow_01);
	_verify(m_assigncol_01);
	_verify(v_norm_01);
	_verify(v_dot_01);
	_verify(v_cross_01);
	
	_verify(r_x_01);
	_verify(r_y_01);
	_verify(r_z_01);
	_verify(accelpointmass_01);
	_verify(cheb3d_01);
	_verify(eccanom_01);
	_verify(frac_01);
	_verify(meanobliquity_01);
	_verify(mjday_01);
	_verify(mjday_02);
	_verify(mjday_tdb_01);
	_verify(position_01);
	_verify(timediff_01);
	*/
    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED/n");

    printf("Tests run: %d/n", tests_run);
    return 0;
}