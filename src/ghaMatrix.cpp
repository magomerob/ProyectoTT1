// $Source$
/**
 * @file ghaMatrix.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de ghaMatrix.
 * @version 0.1
 * @date 2025-05-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/ghaMatrix.h"

Matrix GHAMatrix(double Mjd_UT1)
{
    return r_z(gast(Mjd_UT1));
}