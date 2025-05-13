// $Source$
/**
 * @file nutMatrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de nutMatrix.
 * @version 0.1
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/nutMatrix.h"

Matrix nutMatrix(double Mjd_TT)
{
    double eps = MOblq(Mjd_TT);
    auto [dpsi, deps] = NutAngles(Mjd_TT);

    Matrix nutMat = r_x(-eps-deps)*r_z(-dpsi)*r_x(eps);
    return nutMat;
}