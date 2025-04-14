// $Source$
/**
 * @file r_y.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación del cálculo de la matriz de rotación en el eje y.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/r_y.h"

Matrix r_y(double angle)
{
    double C = cos(angle);
    double S = sin(angle);
    Matrix rotmat = zeros(3,3);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;
    return rotmat;
}