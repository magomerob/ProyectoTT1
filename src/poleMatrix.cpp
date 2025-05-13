// $Source$
/**
 * @file poleMatrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de poleMatrix.
 * @version 0.1
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/poleMatrix.h"

Matrix poleMatrix(double xp,double yp)
{
    Matrix PoleMat = r_y(-xp) * r_x(-yp);
    return PoleMat;
}