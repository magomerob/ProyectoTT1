// $Source$
/**
 * @file g_accelHarmonic.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de g_accelHarmonic.
 * @version 0.1
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/accelHarmonic.h"

Matrix g_accelHarmonic(Matrix r, Matrix u, int n, int m)
{
    double d = 1.0;

    Matrix G = zeros(3,3);
    Matrix dr = zeros(3,1);
    Matrix da;
    for(size_t i=1; i<=3; i++){
        dr = zeros(3,1);
        dr(i,1) = d;
        
        Matrix v1 = r+dr/2;
        Matrix v2 = r-dr/2;
        da = accelHarmonic (v1, u, n, m ) -
             accelHarmonic (v2, u, n, m );

        G.assign_column(transpose((da)/d), i);

    }
    return G;
}