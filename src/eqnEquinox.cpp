// $Source$
/**
 * @file eqnEquinox.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de eqnEquinox.
 * @version 0.1
 * @date 2025-05-13
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/eqnEquinox.h"

double eqnEquinox(double Mjd_TT)
{
    auto [dpsi, deps] = NutAngles(Mjd_TT);
    return dpsi * cos(MOblq(Mjd_TT));
}