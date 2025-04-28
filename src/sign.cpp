// $Source$
/**
 * @file sign.cpp
 * @author Marcos Gómez Robres
 * @brief Implementación de la función sign.
 * @version 0.1
 * @date 2025-04-14
 * 
 * @copyright Copyright (c) 2025
 * @bug No known bugs
 */

#include "../include/sign.h"

double sign(double a, double b)
{
    if(b>=0){
        return fabs(a);
    }
    return -1*fabs(a);
}