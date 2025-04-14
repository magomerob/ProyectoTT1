// $Source$
/**
 * @file frac.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación el calculo de la parte fraccionaria.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/frac.h"

double frac(double x)
{
    return x-floor(x);
}